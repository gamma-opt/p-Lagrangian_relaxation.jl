include("Models_generation.jl")

include("p_lagrangian_decomposition_Budnle.jl")

function dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity )

    # intialisation of the parameters
    precision_p_iteration = -1 .* ones(1, number_of_continuous_decision_variables)
    iteration = 0 # iteration counter
    UB = Inf # initial value for the upper bound to enter the loop
    LB = - Inf # initial value for the lower bound to enter the loop
    f_rank = zeros(1, number_of_continuous_decision_variables)

    results = zeros(3, max_number_of_iterations)

    original_problem, objective_Qs = original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
    init_time  = time()

    while ( iteration == 0 ? true : ( (UB - LB)/LB > tolerance ) ) & ( time() - init_time < time_limit ) & (iteration < max_number_of_iterations)

        iteration  = iteration + 1

        dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p_iteration, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
        optimize!(dynamic_RNMDT_problem)
        UB = objective_value(dynamic_RNMDT_problem)

        #stroing the optimal values of the decision variables and auxiliary
        # variable w represeting the qudratic terms
        xr = value.(dynamic_RNMDT_problem[:x])
        yr = value.(dynamic_RNMDT_problem[:y])
        w_RNMDT_r = value.(dynamic_RNMDT_problem[:w_RNMDT])

        # setting strating values for the continuous variables
        set_start_value.(original_problem[:y], yr)

        #fixing the values for the integer variables
        fix.(original_problem[:x], xr)
        #set_start_value.(original_problem[:x], xr)

        # optimizing the original problems using the values for the decision
        # variables obtained from the obtimized relaxation problem
        optimize!(original_problem)
        LB  = objective_value(original_problem)

        #print(precision_p_iteration)
        # if iteration is not a multiple of N2
        print("before \n")
        print(precision_p_iteration)
        print("\n")
        if mod(iteration+1, N2) != 0
            [ f_rank[j] = sum( round(1/number_of_scenarios, digits = 2) * objective_Qs[s][i, j] * (w_RNMDT_r[i, j, s] - yr[i, s] * yr[j, s])  for i = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios) for j = 1 : number_of_continuous_decision_variables]
            sorted_N1_indeces = sortperm(f_rank[1 : end], rev = true)[1 : min(N1, number_of_continuous_decision_variables)]
            precision_p_iteration[sorted_N1_indeces] = precision_p_iteration[sorted_N1_indeces] .- 1
        else
            #find indexes with the biggest p
            indixes_max_p = findall(x -> x == maximum(precision_p_iteration), precision_p_iteration)
            precision_p_iteration[indixes_max_p] = precision_p_iteration[indixes_max_p] .- 1
        end
        print("after \n")
        print(precision_p_iteration)
        print("\n")


    results[:, iteration] = [UB, LB, (UB - LB)/LB]

    print("iteration: $(iteration) \n \n " )
    print("the difference: $((UB - LB)/LB) percent \n \n " )
    print("UB: $(UB) \n \n " )
    print("LB: $(LB) \n \n " )
    print("xr = \n")
    print(xr)
    print("\n")
    print("yr = \n")
    print(yr)
    print("\n")

    end
    final_time  = time() - init_time
    return results[:, 1 : iteration], final_time

end



function dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity )

    # intialisation of the parameters
    precision_p_iteration = -1 .* ones(1, number_of_continuous_decision_variables)
    iteration = 0 # iteration counter
    UB = Inf # initial value for the upper bound to enter the loop
    LB = -Inf # initial value for the lower bound to enter the loop
    Max_number_of_iterations_for_LD_bundle = 100

    f_rank = zeros(1, number_of_continuous_decision_variables)

    results = zeros(3, max_number_of_iterations)

    original_problem, objective_Qs, objective_fs, objective_c = original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
    init_time  = time()
    while ( iteration == 0 ? true : ( (UB - LB)/LB > tolerance ) ) & ( time() - init_time < time_limit ) & (iteration < max_number_of_iterations)

        iteration  = iteration + 1

        UB, x_values_LD_RNMDT, y_values_LD_RNMDT, w_RNMDT_values_LD_RNMDT = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p_iteration, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity, Max_number_of_iterations_for_LD_bundle)

        yr = y_values_LD_RNMDT

        w_RNMDT_r =  w_RNMDT_values_LD_RNMDT

        xr = repeat( sum(x_values_LD_RNMDT, dims = 2)./number_of_scenarios, 1, number_of_scenarios )
        #xr = repeat( Int.(round.( sum(x_values_LD_RNMDT, dims = 2)./number_of_scenarios, digits = 0)), 1, number_of_scenarios)

        # setting strating values for the continuous variables
        set_start_value.(original_problem[:y], yr)

        #fixing the values for the integer variables
        #fix.(original_problem[:x], xr)
        set_start_value.(original_problem[:x], xr)

        optimize!(original_problem)
        LB  = objective_value(original_problem)

        print("before \n")
        print(precision_p_iteration)
        print("\n")

        #print(precision_p_iteration)
        # if iteration is not a multiple of N2
        if mod(iteration+1, N2) != 0
            [ f_rank[j] = sum( round(1/number_of_scenarios, digits = 2) * objective_Qs[s][i, j] * (w_RNMDT_r[i, j, s] - yr[i, s] * yr[j, s])  for i = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios) for j = 1 : number_of_continuous_decision_variables]
            sorted_N1_indeces = sortperm(f_rank[1 : end], rev = true)[1 : min(N1, number_of_continuous_decision_variables)]
            precision_p_iteration[sorted_N1_indeces] = precision_p_iteration[sorted_N1_indeces] .- 1
        else
            #find indexes with the biggest p
            indixes_max_p = findall(x -> x == maximum(precision_p_iteration), precision_p_iteration)
            precision_p_iteration[indixes_max_p] = precision_p_iteration[indixes_max_p] .- 1
        end
        print("after \n")
        print(precision_p_iteration)
        print("\n")


    results[:, iteration] = [UB, LB, (UB - LB)/LB]

    print("iteration: $(iteration) \n \n " )
    print("the difference: $((UB - LB)/LB) percent \n \n " )
    print("UB: $(UB) \n \n " )
    print("LB: $(LB) \n \n " )
    print("xr = \n")
    print(xr)
    print("\n")
    print("yr = \n")
    print(yr)
    print("\n")
    end

    final_time  = time() - init_time

    return results[:, 1 : iteration], final_time

end



#for i = 1:5

    #RNMDT[i]  = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[i], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity )

#end


#using Plots
#plot(collect(1:length(LD_RNDT[1][1][1, :])), LD_RNDT[1][1][1, :], )
#plot!( collect(1:length(LD_RNDT[1][1][2, :])), LD_RNDT[1][1][2, :],)

#using XLSX

#for scenario = 1 : 5
    #df = DataFrame(UB = Float64[], LB_primal = Float64[], percent = Float64[], time = Float64[])


    #for i = 1 : length(RNMDT[scenario][1][1, :])
        #push!(df, (RNMDT[scenario][1][1, i], RNMDT[scenario][1][2, i], RNMDT[scenario][1][3, i], RNMDT[scenario][2]) )
    #end

    #XLSX.writetable("scenario_$scenario.xlsx", df, names(df))
#end


#test  = RNMDT_problem_generation()
