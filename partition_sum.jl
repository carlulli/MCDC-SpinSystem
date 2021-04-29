L = 3 # dimensions of the 2 by 2 matrix (grid)
# this here is my way of finding all possible states  
# ! the way its currently coded , you have to put manually 2^L times (-1,1) in here 
a = collect(Iterators.product((-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)))

#for b in a 
 #   println(a)
#end

InteractionStrength = 0.001 # the interaction strength between the spins
Temperature_factor = 0.01   # the temperature of the system
magField = 10  # the magnetic field
magnetization_density_vector = zeros(convert(Int,(end_temp - start_temp)/temp_step)+1)

# summing over different tempertures 
# these are not the total temperaures though, as these still get multiplied by a factor 
start_temp = 1
temp_step = 1
end_temp = 3000

for s in start_temp:temp_step:end_temp

    state = zeros(L,L)
    energy_vector = zeros(2^(L*L))
    partition_sum = 0
    state_magnetization = zeros(2^(L*L))
    Temperature = s* Temperature_factor
    #magnetization_density_vector = zeros(300)

    # here I create and sum over all possible configurations. 
    @inbounds for k in 1:(2^(L*L))
        tupel = a[k]
        start = 1
        energy = 0
        # first I translate the numbers in one tupel into one configuration
        @inbounds for i in 1:L, j in 1:L
            state[i,j] = tupel[start]
            start += 1 
        end
        #println(state)
        #println(sum(state))
        energy = sum(state)*magField   # the sum over all spins times the magnetic field constitues one part of the energy
        state_magnetization[k] = sum(state) # the sum is also the magnetization 
        
        # interaction energy with periodic boundary condition

        @inbounds for l in 1:L, n in 1:L
            energy += -InteractionStrength*state[l,n]*(state[l%L+1,n]+state[l,n%L+1]) # 
        end
        #println(energy)
        energy_vector[k] = energy
        partition_sum += exp(- (energy / (Temperature)))
        #println(partition_sum)
    end

    # function for probability of being in state r


    # summing over all states multiplying each states magnetization with the probability of being in that state
    total_magnetization = 0
    @inbounds for x in 1:(2^(L*L))
        total_magnetization += state_magnetization[x] * (1/partition_sum) * exp(energy_vector[x]/(Temperature))
    end
    total_magnetization
    #println(total_magnetization)
    magnetization_density = total_magnetization/(L*L)
    println(magnetization_density)
    magnetization_density_vector[s] = magnetization_density
    
    # here I want to test if the total probability turns out to be 1
    #= 
    prob = 0
    @inbounds for x in 1:(2^(L*L))
        prob += (1/partition_sum) * exp(energy_vector[x]/(Temperature))
    end
    prob =#
    
    #println("The probability is")
    #println(prob)

end 
x = start_temp:temp_step:end_temp
using Plots
magnetization_density_vector
plot( (Temperature_factor * x),magnetization_density_vector)
#= plot( (Temperature_factor * x).^-1,magnetization_density_vector)
xlims!(0,1) =#