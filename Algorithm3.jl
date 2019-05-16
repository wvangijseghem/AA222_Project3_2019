using LinearAlgebra
using Distributions
using PyPlot
using Random

function RouteTimes(r, properties, locations)
    loiter = properties[1]; v = properties[2]
    delivery_time = 0
    arrival_time = 0
    schedule_list = Any[]
    for k = 2:length(r)
        i = r[k]
        j = r[k-1]
        if (i != 0) || (j != 0)
            delivery_time = arrival_time
            edge_time = loiter + distance(i,j, locations)/v
            arrival_time = arrival_time + edge_time
            if (i == 0) && (j != 0)
                push!(schedule_list,(delivery_time,arrival_time))
                arrival_time = 0
            end
        end
    end
    return schedule_list
end

function distance(i,j, locations)
    if i != 0 && j != 0
        loc_1= locations[i,:]
        loc_2= locations[j,:]
        return norm(loc_2 - loc_1)
    elseif j == 0
        return norm(locations[i,:])
    else
        return norm(locations[j,:])
    end
end

function cost(r, objective, constraints, properties, locations, payload)
    budget = constraints[1]; T = constraints[2]; penalty = properties[8]
    energy_cost = EnergyCost(r, properties, locations,payload)
    drone_cost, delivery_time = dcdt(r,energy_cost,objective, properties, locations, constraints)
    c = energy_cost + drone_cost
    if c > budget
        c = c + penalty*(c-budget)
        delivery_time = delivery_time + penalty*(c-budget)
    elseif delivery_time > T
        c = c + penalty*(t-T)
        delivery_time = delivery_time + penalty*(t-T)
    end
    return c, delivery_time
end

function EnergyCost(r, properties, locations,payload)
    loiter = properties[1]; v = properties[2]; Total_Weight = properties[3]; E_price = properties[6]
    density = properties[7]; penalty = properties[8]
    t = 0
    weight = 0
    c = 0
    weight_time = 0
    for k = length(r)-1:-1:1
        i = r[k+1]
        j = r[k]
        if i != 0 || j != 0
            edge_time = loiter + distance(i,j, locations)/v
            t += edge_time
            weight_time += weight*edge_time
            if j != 0
                weight += payload[j]
            else
                weight = weight
            end
            if j == 0 && i != 0
                E = BatteryEnergy(t, weight_time, properties)
                if E > 0
                    c = c + E * E_price
                else
                    c = c - penalty*(E*E_price)
                end
                if weight + E/density > Total_Weight
                    c = c + penalty*(weight + E/density - Total_Weight)
                end
                t = 0
                weight_time = 0
                weight = 0
            end
        end
    end
    return c
end

function BatteryEnergy(t, weight_time, properties)
    α = properties[8]; β = properties[9]; density = properties[7]
    E = (α*weight_time + β*t)/(1 - α/density*t)
    return E
end


function routeschedule(u, n)
  # u: contains pairs of delivery and arrival times. (Array of tuples)
  # n: number of drones. (Integer)

  k = zeros((Int(n), 2)) # Initialize the drone timing vector (Rows: Drones; Column1: Delivery Time; Column2: Arrival Time)
  for i = 1 : length(u) # Go through every route in the route timing vector 'u'
    index = argmin( k[:, 2] ) # Get the drone with the smallest arrival time
    k[index, 1] = k[index, 2] + u[i][1]
    k[index, 2] = k[index, 2] + u[i][2]
  end
  return maximum( k[:, 1] ) # Return the maximum delivery time
end

function dcdt(r, energy_cost, objective, properties, locations, constraints)
    # r: Solution Vector. (Array)
    # energy_cost: Cost of energy for solution r. (Float)
    # objective: A boolean variable that is 1 if the goal is to minimize
    # drone cost, or 0 to minimize overall delivery time.

    # Get the vector 'u' that contains pairs of delivery and arrival times
    # 'u' is an array of tuples.
    u = RouteTimes(r, properties, locations)

    # Perform a binary search to find the best number of drones 'n'
    # while minimizing cost
    if (objective == 1)
        # Set the lower and upper bounds rho and psi.
        # properties[5] is the maximum number of drones
        rho = 1
        psi = properties[5]

        while rho <= (psi - 1)
            n = floor(rho + (psi - rho) / 2)
            if (routeschedule(u, n) <= constraints[2])
                psi = n
            else
                rho = n + 1
            end
        end
        n = rho
        if (routeschedule(u, n) > constraints[2])
            n = psi
        end

    else
        # Purchase as many drones as possible to minimize time, given a budget
        # (constraint[1]) and a cost (properties[4]) of a single drone
        n = (constraints[1] - energy_cost) / properties[4]
        # Make sure number of drones is positive
        if (n < 1)
            n = 1
        end
    end

    # Find the cost of drones and overall delivery time
    drone_cost = n * properties[4]
    overall_delivery_time = routeschedule(u, n)

    return drone_cost, overall_delivery_time
end

function manipulations(r, i, j, type)
    s = copy(r)
    if (type == 1)
        swap = s[i]
        s[i] = s[j]
        s[j] = swap
    elseif (type == 2)
        s = vec(s)
        deleted = s[i]
        deleteat!(s, i)
        splice!(s, j - 1, [s[j - 1], deleted])
    else
        s = vec(s)
        s = reverse(s, i, j)
    end
    return s
end

function sa(sa_properties, objective, N, constraints, properties, locations, payload)
    initial_temp = sa_properties[1]; final_temp = sa_properties[2]
    reduction_factor = sa_properties[3]; number_adjustments = sa_properties[4]
    num_locations = N
    a = collect(1 : num_locations)
    b = fill(0, (N-1, 1))
    r = [0; a; b; 0]

    current_temp = initial_temp
    lst = []

    while current_temp > final_temp
        current_temp = reduction_factor * current_temp
        for k = 1 : number_adjustments
            i = rand(2 : (length(r) - 1))
            j = rand(2 : (length(r) - 1))
            type  = rand(1 : 3)
            r_prime = manipulations(r, i, j, type)
            x = rand()
            if objective == 0
                ind = 2
            else
                ind = 1
            end
            c = cost(r, objective, constraints, properties, locations, payload)
            c_prime = cost(r_prime, objective, constraints, properties, locations, payload)

            if exp(c[ind] - c_prime[ind]) / current_temp >= x
                r = r_prime
                push!(lst,c_prime[ind])
            end
        end
    end
    print(r)
    return r, lst
end

N = 30
Random.seed!(1234)
locations = rand(-100:100,N,2)
delivery_payload = rand(Uniform(0,3),N)
tau = 60 #1
v = 6 #2
Total_Weight = 3 #3
Drone_Cost = 500 #4
n_drones = 10 #5
E_price = 0.1 #6
density = 650 #7
α = 0.217 #8
β = 0.185 #9
penalty = 10000000 #10
objective = 1
constraints = [1500, 10*60]
properties = [tau, v, Total_Weight, Drone_Cost, n_drones, E_price, density, α, β, penalty]
sa_properties = [1, 0.001, 0.9, 1000]


r = [0,1,2,3,0,0,0,8,9,10,0,4,5,6,7,0]
energy_cost = EnergyCost(r, properties, locations,delivery_payload)
c,l = dcdt(r, energy_cost, objective, properties, locations, constraints)
c, delivery_time = cost(r, objective, constraints, properties, locations, delivery_payload)


r, lst = sa(sa_properties, objective, N, constraints, properties, locations, delivery_payload)
figure1 = figure()
ax = figure1.add_subplot(1,1,1)
ax.plot(lst)
title("Convergence cost SA")
ylabel("Cost")
xlabel("Iterations")
gcf()
show()
savefig("SA.png")
