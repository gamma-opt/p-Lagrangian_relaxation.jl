using SharedArrays, Distributed

include("Models_generation.jl")


# the function containing the implementation of the Lagrangian decomposition
# with the bundle method multipliers update applied to the mixed integer based
# relaxation of the original problem

function Lagrangian_decomposition_bundle(p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # limits for the initial values of the Lagrangian multipliers
    lagrangian_multipliers_max = 2
    lagrangian_multipliers_min = 0

    # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
    f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
    [ vector_of_lambda_lagrangian[1, i] = lagrangian_multipliers_min .+ (lagrangian_multipliers_max - lagrangian_multipliers_min)  .* rand(1,
        number_of_continuous_decision_variables + number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
    x_indices = 1 : number_of_continuous_decision_variables

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
    y_indices = number_of_continuous_decision_variables + 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables

    # dual function at the lagragian multiplers' vector at correspondent iteration
    dual_objective_value_at_lagrangian = Array{Float64}(undef, 1, number_of_iterations)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    decision_variables_values_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

    # values at each iteration of the variable z uder minimization in the objective function of the cutting plane method
    relaxed_dual_objective_value = Array{Float64}(undef, 1, number_of_iterations)

    # the center of gravity at each iteration
    center_of_gravity = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    # dual function at the center of gravity at correspondent iteration
    dual_fiunction_value_at_the_center_of_gravity = Array{Float64}(undef, 1, number_of_iterations)

    # subgradient vector at each iteration
    subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    #array of the subproblems for every scenario with lambda inital
    subproblems = LD_RNDMT_problem_generation(p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # upper bound for the original problem obtained by solving relaxed subproblems
    # and summing up the values of the objective functions
    ub_of_original_problem = SharedArray{Float64}(1,1)

    cutting_plane_subproblem = Model(with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
    @variables cutting_plane_subproblem begin
        z
        lagrangian_multipliers_representing_variable[ 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables,
            1 : number_of_scenarios - 1]
        end

    # the values for the parameters of the Bundle method
    m = 0.2
    d = 0.01

    iteration = 1 #strating counter

    eps_stop = 1e-2 # stopping criteria

    initial_time = time()
    while (iteration <= number_of_iterations) & ((iteration>6) ? (norm(dual_objective_value_at_lagrangian[iteration - 5:iteration-1].- dual_objective_value_at_lagrangian[iteration-6:iteration-2])^2 >= eps_stop) : true)

        ub_of_original_problem[1] = 0

        @sync @distributed for s in 1 : number_of_scenarios

            #objective_update
            @objective( subproblems[s], Max,
                (1/number_of_scenarios) *
                    (
                    sum( objective_Qs[s][i, i] * subproblems[s][:w_RNMDT][i, i]
                        for i = 1 : number_of_continuous_decision_variables )
                    + 2 * sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                        for i = 1 : number_of_continuous_decision_variables,
                            j = i+1 : number_of_continuous_decision_variables)
                    + sum( subproblems[s][:x][i] * objective_fs[s][1, i]  for i = 1:number_of_continuous_decision_variables)
                    + sum( subproblems[s][:y][j] * objective_fs[s][2, j]  for j = 1:number_of_integer_decision_variables)
                    )
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

        iteration += 1

    end

    final_time = time()-initial_time

    return [dual_objective_value_at_lagrangian[1 : ( (iteration > number_of_iterations) ? number_of_iterations : iteration )], final_time, dual_fiunction_value_at_the_center_of_gravity]

end

function dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # limits for the initial values of the Lagrangian multipliers
    lagrangian_multipliers_max = 10
    lagrangian_multipliers_min = 0

    # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
    f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
    [ vector_of_lambda_lagrangian[1, i] = lagrangian_multipliers_min .+ (lagrangian_multipliers_max - lagrangian_multipliers_min)  .* rand(1,
        number_of_continuous_decision_variables + number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
    x_indices = 1 : number_of_continuous_decision_variables

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
    y_indices = number_of_continuous_decision_variables + 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables

    # dual function at the lagragian multiplers' vector at correspondent iteration
    dual_objective_value_at_lagrangian = Array{Float64}(undef, 1, number_of_iterations)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    decision_variables_values_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

    # the vector that contains quadraticity representing variables wij for each of the scenarios
    RNMDT_quadraticity_variables_w_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables, number_of_continuous_decision_variables, number_of_scenarios)

    # values at each iteration of the variable z uder minimization in the objective function of the cutting plane method
    relaxed_dual_objective_value = Array{Float64}(undef, 1, number_of_iterations)

    # the center of gravity at each iteration
    center_of_gravity = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    # dual function at the center of gravity at correspondent iteration
    dual_fiunction_value_at_the_center_of_gravity = Array{Float64}(undef, 1, number_of_iterations)

    # subgradient vector at each iteration
    subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    #array of the subproblems for every scenario with lambda inital
    subproblems = dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # upper bound for the original problem obtained by solving relaxed subproblems
    # and summing up the values of the objective functions
    ub_of_original_problem = SharedArray{Float64}(1,1)

    cutting_plane_subproblem = Model(with_optimizer(Gurobi.Optimizer, OutputFlag = 0))
    @variables cutting_plane_subproblem begin
        z
        lagrangian_multipliers_representing_variable[ 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables,
            1 : number_of_scenarios - 1]
        end

    # the values for the parameters of the Bundle method
    m = 0.9
    d = 0.05
    N_to_change_bundle_parameters = 1

    iteration = 0 #strating counter

    eps_stop = 1e-1 # stopping criteria

    initial_time = time()
    while (iteration < number_of_iterations) & ((iteration>6) ? (norm(dual_objective_value_at_lagrangian[iteration - 5:iteration-1].- dual_objective_value_at_lagrangian[iteration-6:iteration-2])^2 >= eps_stop) : true)

        iteration += 1
        print("iteration = $iteration \n ")

        ub_of_original_problem[1] = 0

        @sync @distributed for s in 1 : number_of_scenarios

            #objective_update
            @objective( subproblems[s], Max,
                (1/number_of_scenarios) *
                    (
                    sum( objective_Qs[s][i, i] * subproblems[s][:w_RNMDT][i, i]
                        for i = 1 : number_of_continuous_decision_variables )
                    + 2 * sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                        for i = 1 : number_of_continuous_decision_variables,
                            j = i+1 : number_of_continuous_decision_variables)
                    + sum( subproblems[s][:x][i] * objective_fs[s][1, i]  for i = 1:number_of_continuous_decision_variables)
                    + sum( subproblems[s][:y][j] * objective_fs[s][2, j]  for j = 1:number_of_integer_decision_variables)
                    )
                    +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[x_indices] .* subproblems[s][:x] )

                    +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[y_indices] .* subproblems[s][:y] )
            )

            status = optimize!(subproblems[s])
            obj_value = objective_value(subproblems[s])
            ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value

            decision_variables_values_for_each_scenario[ x_indices, s ] = value.(subproblems[s][:x])
            decision_variables_values_for_each_scenario[ y_indices, s ] = value.(subproblems[s][:y])
            RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(subproblems[s][:w_RNMDT])
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
        if mod(iteration, N_to_change_bundle_parameters) == 0
            m = m * 0.8
            d = d * 0.8
        end

    end

    final_time = time()-initial_time

    #return [dual_objective_value_at_lagrangian[1 : ( (iteration > number_of_iterations) ? number_of_iterations : iteration )], final_time, dual_fiunction_value_at_the_center_of_gravity]
    return dual_objective_value_at_lagrangian[ iteration ], decision_variables_values_for_each_scenario[x_indices, :], decision_variables_values_for_each_scenario[y_indices, :], RNMDT_quadraticity_variables_w_for_each_scenario

end

#result = p_Lagrangian_decomposition_bundle(-4, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

#figure("Budnle method")
#plot(Array(1:iteration-1),dual_objective_value_at_lagrangian[1: iteration-1])
#scatter(Array(1:iteration-1), dual_objective_value_at_lagrangian[1: iteration-1])
#title("Bundle method insipred lagrangian decomposition, $number_of_scenarios scenarios")
#xlabel("iteration")
#ylabel("Upper bound")
#savefig("Bundle method3")
