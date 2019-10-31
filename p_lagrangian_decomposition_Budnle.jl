using SharedArrays, Distributed, PyPlot

include("models_generation.jl")

#intialization of the parameters for the bundle method
Max_iterations = 20

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

# lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
vector_of_lambda_lagrangian = Array{Any}(undef, Max_iterations, (number_of_scenarios - 1) * 2)
[ vector_of_lambda_lagrangian[1, i] = 10.5 .+ 60.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:(number_of_scenarios - 1) * 2 ]

# lagrangian relaxation variables for the y-non anticipativity conditions
#vector_of_mu_lagrangian = Array{Any}(undef, 1, number_of_scenarios-1)
#[ vector_of_mu_lagrangian[i] = 0.5 .+ 1.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:number_of_scenarios-1 ]


# dual function at the lagragian multiplers'vector at correspondent iteration
lagrangian_bundle_phi = Array{Float64}(undef, 1, Max_iterations)

# vector that contains decision variables written in a column (x and y in this case)
# (each row represnets the components of the correspondent lagrangian multiplier)
z_decomposition = SharedArray{Float64}(number_of_continuos_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

# values at each ietration of the variable z uder minimization in the objective function of the cutting plane method
z_bundle = Array{Float64}(undef, 1, Max_iterations)

# the center of gravity at each iteration
M_bundle = Array{Any}(undef, Max_iterations, (number_of_scenarios - 1) * 2)

# dual function at the center of gravity at correspondent iteration
M_bundle_phi = Array{Float64}(undef, 1, Max_iterations)

# subgradient at iteration k
subgradient_bundle = Array{Any}(undef, Max_iterations, (number_of_scenarios - 1) * 2)


#bundle_LB = Array{Float64}(undef, 1, Max_iterations)

#array of the models for every scenario with lambda inital
models = subproblem
#[models[m] =  SW_PC_RNMDT_DEC_generation(m, lambda_bundle[1,:,:,:,:]) for m in Mn]

UB_decomposition = SharedArray{Float64}(1,1)
LB_decomposition = SharedArray{Float64}(1,1)

bundle_subproblem = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0))
@variables bundle_subproblem begin
    z
    mu[1:number_of_continuos_decision_variables, 1:(number_of_scenarios - 1) * 2]
end

m = 0.3
d = 0.2
mu0 = 3 #for the initialization of the M_bundle
ii = 1 #strating counter
initial_time = time()
while (ii <= Max_iterations) #& ((ii>6) ? (norm(lambda_bundle_phi[ii - 5:ii-1].- lambda_bundle_phi[ii-6:ii-2])^2 >= eps_stop) : true)

    UB_decomposition[1] = 0
     for s in 1 : number_of_scenarios
    #objective_update
    @objective( models[s], Max,
            sum( objective_Qs[s][i, i] * models[s][:w_RNMDT][i, i]
                for i = 1 : number_of_continuos_decision_variables )
            + 2 * sum(objective_Qs[s][i, j] * models[s][:w_RNMDT][i, j]
                for i = 1 : number_of_continuos_decision_variables,
                    j = i+1 : number_of_continuos_decision_variables)
            + sum( ( models[s][:x] .* objective_fs[s][1, :] )
                .+ ( models[s][:y] .* objective_fs[s][2, :] ) )
                + objective_fs[s][3, 1]
            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[ii, 1 : number_of_scenarios - 1], s ) .* models[s][:x] )

            +  sum( f_lambda_lagrangian(vector_of_lambda_lagrangian[ii, number_of_scenarios - 1 + 1 :  (number_of_scenarios - 1) * 2], s ) .* models[s][:y] )
                )

    status = optimize!(models[s])
    obj_value = objective_value(models[s])
    UB_decomposition[1] = UB_decomposition[1] + obj_value

    z_decomposition[ 1:number_of_continuos_decision_variables, s] = value.(models[s][:x])
    z_decomposition[ number_of_continuos_decision_variables + 1 : number_of_continuos_decision_variables + number_of_integer_decision_variables ] = value.(models[s][:y])
    end

        lagrangian_bundle_phi[ii] = UB_decomposition[1]

        [subgradient_bundle[ii,  s - 1] =  z_decomposition[ 1 : number_of_continuos_decision_variables, 1] - z_decomposition[ 1 : number_of_continuos_decision_variables, s] for  s in 2: number_of_scenarios ]
        [subgradient_bundle[ii,  end - (number_of_scenarios - s)] =  z_decomposition[ 1 + number_of_continuos_decision_variables : number_of_continuos_decision_variables + number_of_integer_decision_variables, 1] - z_decomposition[ 1 + number_of_continuos_decision_variables : number_of_continuos_decision_variables + number_of_integer_decision_variables, s] for  s in 2: number_of_scenarios ]

        if ii == 1
            M_bundle[ii, :] = vector_of_lambda_lagrangian[ii, :]
            M_bundle_phi[ii] = lagrangian_bundle_phi[ii]
        elseif lagrangian_bundle_phi[ii] - M_bundle_phi[ii-1] <= m * ( z_bundle[ii-1] - d * (ii-1) * sum( sum( vector_of_lambda_lagrangian[ii-1, s] .- M_bundle[ii-1, s] ) .^ 2  for s = 1 : (number_of_scenarios - 1) * 2 ) - lagrangian_bundle_phi[ii-1] )
            M_bundle[ii, :] = vector_of_lambda_lagrangian[ii, :]
            M_bundle_phi[ii] = lagrangian_bundle_phi[ii]
        else
            M_bundle[ii, :] = M_bundle[ii-1, :]
            M_bundle_phi[ii] = M_bundle_phi[ii-1]
        end

        @objective(bundle_subproblem, Min, z + d * ii * sum( sum( ( bundle_subproblem[:mu][:, s] .- M_bundle[ii, s] ).^2 ) for  s in 1 : (number_of_scenarios - 1) * 2) )
        @constraint(bundle_subproblem, z >= lagrangian_bundle_phi[ii] + sum( sum( subgradient_bundle[ii, s] .* ( bundle_subproblem[:mu][s, :] .- vector_of_lambda_lagrangian[ii, s] ) ) for s = 1 : (number_of_scenarios - 1) * 2 ) )
        status = optimize!(bundle_subproblem)
            if ii < Max_iterations
                [ vector_of_lambda_lagrangian[ii+1, s] = value.(bundle_subproblem[:mu][:, s]) for s = 1 : (number_of_scenarios - 1) * 2 ]
                z_bundle[ii] = value.(bundle_subproblem[:z])
            end
        global ii += 1
    end
    final_time = time()-initial_time

    figure("Budnle method")
    plot(Array(1:ii-1),lagrangian_bundle_phi[1: ii-1])
    scatter(Array(1:ii-1), lagrangian_bundle_phi[1: ii-1])
    title("Bundle method insipred lagrangian decomposition, $number_of_scenarios scenarios")
    xlabel("iteration")
    ylabel("Upper bound")
    savefig("Bundle method1")
