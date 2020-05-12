
function dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed )

    # intialisation of the parameters
    precision_p_iteration = -1 .* ones( number_of_continuous_decision_variables, number_of_scenarios)
    iteration = 0 # iteration counter
    UB_optimal = Inf # initial value for the upper bound to enter the loop
    LB_optimal = - Inf # initial value for the lower bound to enter the loop
    f_rank = zeros(number_of_continuous_decision_variables, number_of_scenarios)

    results = zeros(4, max_number_of_iterations)


    original_problem, objective_Qs, objective_fs, objective_c, constraint_Qs = original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "yes", time_limit, seed)
    init_time  = time()

    while ( iteration == 0 ? true : ( (UB_optimal - LB_optimal)/LB_optimal > tolerance ) ) & ( time() - init_time < time_limit ) & (iteration < max_number_of_iterations)

        iteration  = iteration + 1

        dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p_iteration, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, time_limit - (time() - init_time), seed )
        optimize!(dynamic_RNMDT_problem)
        UB_new = objective_value(dynamic_RNMDT_problem)

        #storing the optimal values of the decision variables and auxiliary
        # variable w represeting the qudratic terms
        xr_new = value.(dynamic_RNMDT_problem[:x])
        yr_new = value.(dynamic_RNMDT_problem[:y])
        w_RNMDT_r_new = value.(dynamic_RNMDT_problem[:w_RNMDT])

        if (UB_new < UB_optimal)
            UB_optimal = UB_new
            xr_UB_optimal = xr_new
            yr_UB_optimal = yr_new
        end

        # setting strating values for the continuous variables
        #set_start_value.(original_problem[:y], yr_new)

        #fixing the values for the integer variables
        fix.(original_problem[:x], xr_new, force = true)

        # optimizing the original problems using the values for the decision
        # variables obtained from the obtimized relaxation problem
        optimize!(original_problem)
        LB_new  = objective_value(original_problem)

        if (LB_new > LB_optimal)
            LB_optimal = LB_new
            xr_LB_optimal = value.(original_problem[:x])
            yr_LB_optimal = value.(original_problem[:y])
        end

        # if iteration is not a multiple of N2
        if mod(iteration+1, N2) != 0
            [ f_rank[j, s] = sum( abs(constraint_Qs[s][r][i, j] * (w_RNMDT_r_new[i, j, s] - yr_new[i, s] * yr_new[j, s]))  for i = 1 : number_of_continuous_decision_variables, r = 1 : number_of_constraints) for j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios]
            for i  = 1 : min(N1, number_of_continuous_decision_variables * number_of_scenarios)
                precision_p_iteration[findall(x->x==maximum(f_rank),f_rank)[1]] = precision_p_iteration[findall(x->x==maximum(f_rank),f_rank)[1]] - 1
                f_rank[findall(x->x==maximum(f_rank),f_rank)[1]] = - Inf
            end

        else
            #find indexes with the biggest p
            indixes_max_p = findall(x -> x == maximum(precision_p_iteration), precision_p_iteration)
            precision_p_iteration[indixes_max_p] = precision_p_iteration[indixes_max_p] .- 1
        end



    results[:, iteration] = [UB_optimal, LB_optimal, (UB_optimal - LB_optimal)/LB_optimal, time() - init_time]
    end

    return results[:, 1 : iteration]

end



function dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, parallelised, seed)

    # intialisation of the parameters
    precision_p_iteration = -1 .* ones( number_of_continuous_decision_variables, number_of_scenarios)
    iteration = 0 # iteration counter
    UB_optimal = Inf # initial value for the upper bound to enter the loop
    LB_optimal = -Inf # initial value for the lower bound to enter the loop
    Max_number_of_iterations_for_LD_bundle = 100

    f_rank = zeros(number_of_continuous_decision_variables, number_of_scenarios)

    # generating the inital values for the center of gravity
    center_of_gravity_min = 0
    center_of_gravity_max = 0

    Random.seed!(seed)
    center_of_gravity_inital_value = Array{Any}(undef, number_of_scenarios - 1)
    [ center_of_gravity_inital_value[i] = center_of_gravity_min .+ (center_of_gravity_max - center_of_gravity_min)  .* rand(1,
        number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]


    results = zeros(6, max_number_of_iterations)
    x_results  = zeros(number_of_integer_decision_variables, number_of_scenarios, max_number_of_iterations)

    original_problem, objective_Qs, objective_fs, objective_c, constraint_Qs = original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "yes", time_limit, seed)
    init_time  = time()
    while ( iteration == 0 ? true : ( (UB_optimal - LB_optimal)/LB_optimal > tolerance ) ) & ( time() - init_time < time_limit ) & (iteration < max_number_of_iterations)

        iteration  = iteration + 1

        UB_new, xr_new, yr_new, w_RNMDT_r_new, lambda, number_of_the_serious_steps, center_of_gravity_inital_value = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p_iteration, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, Max_number_of_iterations_for_LD_bundle, center_of_gravity_inital_value, iteration, parallelised, seed)

        if (UB_new < UB_optimal)
            UB_optimal = UB_new
            xr_UB_optimal = xr_new
            yr_UB_optimal = yr_new
        end

        #xr = repeat( sum(x_values_LD_RNMDT, dims = 2)./number_of_scenarios, 1, number_of_scenarios )
        xr_new = repeat( Int.(round.( sum(xr_new, dims = 2)./number_of_scenarios, digits = 0)), 1, number_of_scenarios)

        x_results[:, :, iteration]  = xr_new

        # setting strating values for the continuous variables
        #set_start_value.(original_problem[:y], yr_new)

        #fixing the values for the integer variables
        fix.(original_problem[:x], xr_new)

        optimize!(original_problem)
        LB_new  = objective_value(original_problem)


        if (LB_new > LB_optimal)
            LB_optimal = LB_new
            xr_LB_optimal = value.(original_problem[:x])
            yr_LB_optimal = value.(original_problem[:y])
        end


        # if iteration is not a multiple of N2
        if mod(iteration+1, N2) != 0
            [ f_rank[j, s] = sum( abs(constraint_Qs[s][r][i, j] * (w_RNMDT_r_new[i, j, s] - yr_new[i, s] * yr_new[j, s]))  for i = 1 : number_of_continuous_decision_variables, r = 1 : number_of_constraints) for j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios]
            #f_rank = reshape(f_rank, 1, number_of_continuous_decision_variables + number_of_scenarios)
            #sorted_N1_indeces = sortperm(f_rank[1 : end], rev = true)[1 : min(N1, number_of_continuous_decision_variables)]
            for i  = 1 : min(N1, number_of_continuous_decision_variables * number_of_scenarios)
                precision_p_iteration[findall(x->x==maximum(f_rank),f_rank)[1]] = precision_p_iteration[findall(x->x==maximum(f_rank),f_rank)[1]] - 1
                f_rank[findall(x->x==maximum(f_rank),f_rank)[1]] = - Inf
            end

        else
            #find indexes with the biggest p
            indixes_max_p = findall(x -> x == maximum(precision_p_iteration), precision_p_iteration)
            precision_p_iteration[indixes_max_p] = precision_p_iteration[indixes_max_p] .- 1

        end


    results[:, iteration] = [UB_optimal, LB_optimal, (UB_optimal - LB_optimal)/LB_optimal, size(lambda, 1), number_of_the_serious_steps, time() - init_time]

    end

    return results[:, 1 : iteration]

end
