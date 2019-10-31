using SharedArrays, Distributed, PyPlot

include("Models_generation.jl")

#intialization of the parameters for the bundle method
number_of_iterations = 20

# auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

# lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
[ vector_of_lambda_lagrangian[1, i] = -10.5 .+ 20.0 .* rand(1,
    number_of_continuos_decision_variables + number_of_integer_decision_variables)
        for i = 1 : number_of_scenarios - 1 ]

# indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
x_indices = 1 : number_of_continuos_decision_variables

# indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
y_indices = number_of_continuos_decision_variables + 1 : number_of_continuos_decision_variables + number_of_integer_decision_variables

# dual function at the lagragian multiplers' vector at correspondent iteration
dual_objective_value_at_lagrangian = Array{Float64}(undef, 1, number_of_iterations)

# vector that contains decision variables written in a column (x and y in this case)
# (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
decision_variables_values_for_each_scenario = SharedArray{Float64}(number_of_continuos_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

# values at each ietration of the variable z uder minimization in the objective function of the cutting plane method
relaxed_dual_objective_value = Array{Float64}(undef, 1, number_of_iterations)

# the center of gravity at each iteration
center_of_gravity = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

# dual function at the center of gravity at correspondent iteration
dual_fiunction_value_at_the_center_of_gravity = Array{Float64}(undef, 1, number_of_iterations)

# subgradient vector at each iteration
subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

#array of the subproblems for every scenario with lambda inital
subproblems = subproblem

# upper bound for the original problem obtained by solving relaxed subproblems
# and summing up the values of the objective functions
ub_of_original_problem = SharedArray{Float64}(1,1)

cutting_plane_subproblem = Model(with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
@variables cutting_plane_subproblem begin
    z
    lagrangian_multipliers_representing_variable[ 1 : number_of_continuos_decision_variables + number_of_integer_decision_variables,
        1 : number_of_scenarios - 1]
end

m = 0.7
d = 0.02

iteration = 1 #strating counter

initial_time = time()
while (iteration <= number_of_iterations) #& ((iteration>6) ? (norm(lambda_bundle_phi[iteration - 5:iteration-1].- lambda_bundle_phi[iteration-6:iteration-2])^2 >= eps_stop) : true)

    ub_of_original_problem[1] = 0
    for s in 1 : number_of_scenarios
        #objective_update
        @objective( subproblems[s], Max,
                sum( objective_Qs[s][i, i] * subproblems[s][:w_RNMDT][i, i]
                    for i = 1 : number_of_continuos_decision_variables )
                + 2 * sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                    for i = 1 : number_of_continuos_decision_variables,
                    j = i+1 : number_of_continuos_decision_variables)
                + sum( ( subproblems[s][:x] .* objective_fs[s][1, :] )
                    .+ ( subproblems[s][:y] .* objective_fs[s][2, :] ) )
                    + objective_fs[s][3, 1]
                +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[x_indices] .* subproblems[s][:x] )

                +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[y_indices] .* subproblems[s][:y] )
                )

    status = optimize!(subproblems[s])
    obj_value = objective_value(subproblems[s])
    ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value

    decision_variables_values_for_each_scenario[ x_indices, s ] = value.(subproblems[s][:x])
    decision_variables_values_for_each_scenario[ y_indices, s ] = value.(subproblems[s][:y])
    end

        dual_objective_value_at_lagrangian[iteration] = ub_of_original_problem[1]

        [ subgradient_vector[iteration,  s - 1] =  decision_variables_values_for_each_scenario[ :, 1] - decision_variables_values_for_each_scenario[ :, s] for  s in 2: number_of_scenarios ]

        if iteration == 1
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_fiunction_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]

        elseif abs( dual_objective_value_at_lagrangian[iteration] - dual_fiunction_value_at_the_center_of_gravity[iteration-1]) >= m * abs( relaxed_dual_objective_value[iteration-1] - d * (iteration-1) * sum( sum( (vector_of_lambda_lagrangian[iteration-1, s] .- center_of_gravity[iteration-1, s]) .^ 2 )  for s = 1 : number_of_scenarios - 1 ) - dual_objective_value_at_lagrangian[iteration-1] )
            #print("left part $(dual_objective_value_at_lagrangian[iteration] - dual_fiunction_value_at_the_center_of_gravity[iteration-1]), right part $( m * ( relaxed_dual_objective_value[iteration-1] - d * (iteration-1) * sum( sum( (vector_of_lambda_lagrangian[iteration-1, s] .- center_of_gravity[iteration-1, s]) .^ 2 )  for s = 1 : number_of_scenarios - 1 ) - dual_objective_value_at_lagrangian[iteration-1] ))  \n ")
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_fiunction_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]

        else
            center_of_gravity[iteration, :] = center_of_gravity[iteration-1, :]
            dual_fiunction_value_at_the_center_of_gravity[iteration] = dual_fiunction_value_at_the_center_of_gravity[iteration-1]

        end

        @objective(cutting_plane_subproblem, Min, z + d * iteration * sum( sum( (cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- center_of_gravity[iteration, s] ).^2 ) for  s in 1 : number_of_scenarios - 1 ) )

        @constraint(cutting_plane_subproblem, z >= dual_objective_value_at_lagrangian[iteration] + sum( sum( subgradient_vector[iteration, s] .* ( cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- vector_of_lambda_lagrangian[iteration, s] ) ) for s = 1 : number_of_scenarios - 1) )

        status = optimize!(cutting_plane_subproblem)
            if iteration < number_of_iterations
                [ vector_of_lambda_lagrangian[iteration+1, s] = value.(cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s]) for s = 1 : number_of_scenarios - 1 ]
                relaxed_dual_objective_value[iteration] = value.(cutting_plane_subproblem[:z])
            end
        global iteration += 1
    end
    final_time = time()-initial_time

    figure("Budnle method")
    plot(Array(1:iteration-1),dual_objective_value_at_lagrangian[1: iteration-1])
    scatter(Array(1:iteration-1), dual_objective_value_at_lagrangian[1: iteration-1])
    title("Bundle method insipred lagrangian decomposition, $number_of_scenarios scenarios")
    xlabel("iteration")
    ylabel("Upper bound")
    savefig("Bundle method3")
