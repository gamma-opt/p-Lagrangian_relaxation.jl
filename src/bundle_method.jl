using Suppressor

# the function contains the implementation of the Lagrangian decomposition
# with the bundle method multipliers update applied to the mixed integer based
# relaxation of the original problem

function auxiliary_check(vector_of_lambda_lagrangian, current_teration, iterations_to_be_considered, number_of_scenarios, tolerance)
    result = zeros(number_of_scenarios-1, 1)
    result_tolerance = tolerance .* ones(number_of_scenarios-1, 1)
    for i = 1 : iterations_to_be_considered
        result = result .+ norm.(vector_of_lambda_lagrangian[current_teration - i + 1 , :] .- vector_of_lambda_lagrangian[current_teration - i, :])
    end

    return sum(result .> result_tolerance)/ (number_of_scenarios-1)
end

function dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, number_of_iterations, center_of_gravity_inital_value, dynamic_precision_algorithm_iteration, parallelised, seed)

    number_of_the_serious_steps = 0

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
    f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

    # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
    #[ vector_of_lambda_lagrangian[1, i] = lagrangian_multipliers_min .+ (lagrangian_multipliers_max - lagrangian_multipliers_min)  .* rand(1,
        #number_of_integer_decision_variables)
            #for i = 1 : number_of_scenarios - 1 ]
    [ vector_of_lambda_lagrangian[1, i] = center_of_gravity_inital_value[i] for i = 1 : number_of_scenarios - 1 ]

    # dual function at the lagragian multiplers' vector at correspondent iteration
    dual_objective_value_at_lagrangian = Array{Float64}(undef, 1, number_of_iterations)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    integer_decision_variables_values_for_each_scenario = Array{Float64}(undef, number_of_integer_decision_variables, number_of_scenarios)

    # vector that contains decision variables written in a column (x and y in this case)
    # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
    continuous_decision_variables_values_for_each_scenario = Array{Float64}(undef, number_of_continuous_decision_variables, number_of_scenarios)

    # the vector that contains quadraticity representing variables wij for each of the scenarios
    RNMDT_quadraticity_variables_w_for_each_scenario = Array{Float64}(undef, number_of_continuous_decision_variables, number_of_continuous_decision_variables, number_of_scenarios)

    # values at each iteration of the variable z uder minimization in the objective function of the cutting plane method
    relaxed_dual_objective_value = Array{Float64}(undef, 1, number_of_iterations)

    # the center of gravity at each iteration
    center_of_gravity = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    # dual function at the center of gravity at correspondent iteration
    dual_function_value_at_the_center_of_gravity = Array{Float64}(undef, 1, number_of_iterations)

    # subgradient vector at each iteration
    subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

    #array of the subproblems for every scenario with lambda inital
    subproblems = dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    # upper bound for the original problem obtained by solving relaxed subproblems
    # and summing up the values of the objective functions
    ub_of_original_problem = Array{Float64}(undef,1,1)

    cutting_plane_subproblem = Model(with_optimizer(Gurobi.Optimizer, Threads = 1, Method = 4)) #, LogFile = loglink_par_bundle * "$(number_of_scenarios)_scenarios_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_$(seed)_seed_$(Dates.today())_bundle_LD+RNDMT_par_logfile.txt" ))
    @variables cutting_plane_subproblem begin
        z
        lagrangian_multipliers_representing_variable[ 1 : number_of_integer_decision_variables,
            1 : number_of_scenarios - 1]
    end

    # the values for the parameters of the Bundle method
    m = 0.7
    d = 10

    iteration = 0 #strating counter

    eps_stop = 100

    number_of_iteration_for_checking = 3

    initial_time = time()
    while (iteration < number_of_iterations) & ((iteration > number_of_iteration_for_checking + 1 ) ? ( norm(dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking : iteration-1] .- dual_objective_value_at_lagrangian[iteration - number_of_iteration_for_checking - 1 : iteration - 2]) >= eps_stop) : true)
        iteration += 1
        ub_of_original_problem[1] = 0

        if parallelised  == "parallelised"

                @suppress @sync Threads.@threads for s in 1 : number_of_scenarios
                    #objective_update
                    @objective( subproblems[s], Max,
                        round( (1/number_of_scenarios), digits = 3) *
                        ( sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                            for i = 1 : number_of_continuous_decision_variables,
                                j = 1 : number_of_continuous_decision_variables)
                        + sum( subproblems[s][:x][i] * objective_c[i]  for i = 1:number_of_integer_decision_variables)
                        + sum( subproblems[s][:y][j] * objective_fs[s][j]  for j = 1:number_of_continuous_decision_variables)
                        )
                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* subproblems[s][:x] )
                    )
                    status = optimize!(subproblems[s])
                    obj_value = objective_value(subproblems[s])

                    ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                    integer_decision_variables_values_for_each_scenario[ :, s ] = value.(subproblems[s][:x])
                    continuous_decision_variables_values_for_each_scenario[ :, s ] = value.(subproblems[s][:y])
                    RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(subproblems[s][:w_RNMDT])

            end

        else
            @suppress for s in 1 : number_of_scenarios
                #objective_update
                @objective( subproblems[s], Max,
                    round( (1/number_of_scenarios), digits = 3) *
                    ( sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                        for i = 1 : number_of_continuous_decision_variables,
                            j = 1 : number_of_continuous_decision_variables)
                    + sum( subproblems[s][:x][i] * objective_c[i]  for i = 1:number_of_integer_decision_variables)
                    + sum( subproblems[s][:y][j] * objective_fs[s][j]  for j = 1:number_of_continuous_decision_variables)
                    )
                    +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s ) .* subproblems[s][:x] )
                )
                status = optimize!(subproblems[s])
                obj_value = objective_value(subproblems[s])

                ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value
                integer_decision_variables_values_for_each_scenario[ :, s ] = value.(subproblems[s][:x])
                continuous_decision_variables_values_for_each_scenario[ :, s ] = value.(subproblems[s][:y])
                RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(subproblems[s][:w_RNMDT])

            end
        end

        dual_objective_value_at_lagrangian[iteration] = ub_of_original_problem[1]

        [ subgradient_vector[iteration,  s - 1] =  integer_decision_variables_values_for_each_scenario[ :, 1] - integer_decision_variables_values_for_each_scenario[ :, s] for  s in 2: number_of_scenarios ]

        if iteration == 1
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
    elseif  dual_function_value_at_the_center_of_gravity[iteration-1] - dual_objective_value_at_lagrangian[iteration] >= m * ( dual_function_value_at_the_center_of_gravity[iteration-1] -  ( relaxed_dual_objective_value[iteration-1] + d * sum( sum( (vector_of_lambda_lagrangian[iteration, s] .- center_of_gravity[iteration-1, s]) .^ 2 )  for s = 1 : number_of_scenarios - 1 ) ) )
            center_of_gravity[iteration, :] = vector_of_lambda_lagrangian[iteration, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_objective_value_at_lagrangian[iteration]
            number_of_the_serious_steps = number_of_the_serious_steps + 1
        else
            center_of_gravity[iteration, :] = center_of_gravity[iteration-1, :]
            dual_function_value_at_the_center_of_gravity[iteration] = dual_function_value_at_the_center_of_gravity[iteration-1]

        end

        @objective(cutting_plane_subproblem, Min, cutting_plane_subproblem[:z] + d * sum( sum( (cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- center_of_gravity[iteration, s] ).^2 ) for  s in 1 : number_of_scenarios - 1 ) )

        @constraint(cutting_plane_subproblem, cutting_plane_subproblem[:z] >= dual_function_value_at_the_center_of_gravity[iteration] + sum( sum( subgradient_vector[iteration, s] .* ( cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s] .- vector_of_lambda_lagrangian[iteration, s] ) ) for s = 1 : number_of_scenarios - 1) )

        status = optimize!(cutting_plane_subproblem)

            if iteration < number_of_iterations
                [ vector_of_lambda_lagrangian[iteration+1, s] = value.(cutting_plane_subproblem[:lagrangian_multipliers_representing_variable][:, s]) for s = 1 : number_of_scenarios - 1 ]
                relaxed_dual_objective_value[iteration] = value.(cutting_plane_subproblem[:z])
            end

    end

    final_time = time()-initial_time

    return dual_objective_value_at_lagrangian[ iteration], integer_decision_variables_values_for_each_scenario[:, :], continuous_decision_variables_values_for_each_scenario[:, :], RNMDT_quadraticity_variables_w_for_each_scenario, vector_of_lambda_lagrangian[ 1 : iteration, :], number_of_the_serious_steps, center_of_gravity[iteration, :]

end
