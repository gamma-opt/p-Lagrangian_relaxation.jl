include("Models_generation.jl")

function dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity )

    # intialisation of the parameters
    precision_p_iteration = -1 .* ones(1, number_of_continuos_decision_variables)
    iteration = 0 # iteration counter
    UB = Inf # initial value for the upper bound to enter the loop
    LB = -Inf # initial value for the lower bound to enter the loop

    original_problem = original_problem_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    while UB - LB > tolerance
        iteration  = iteration + 1
        dynamic_RNMDT_problem, objective_Qs, objective_fs = dynamic_precision_RNMDT_problem_generation(precision_p_iteration, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
        optimize!(dynamic_RNMDT_problem)
        UB = objective_value(dynamic_RNMDT_problem)
        xr = value.(dynamic_RNMDT_problem[:x])
        yr = value.(dynamic_RNMDT_problem[:y])
        w_RNMDT_r = value.(dynamic_RNMDT_problem[:w_RNMDT])

        set_start_value.(original_problem[:x], xr)
        fix.(original_problem[:y], yr)
        optimize!(original_problem)
        LB  = objective_value(original_problem)

        precision_p_iteration = precision_p_iteration .- 1

        #print(precision_p_iteration)

        if mod(iteration, N2) == 0
            [ f_rank[j] = sum( (1/number_of_scenarios) * objective_Qs[s][i, j] * (w_RNMDT_r[i, j, s] - xr[i, s] * xr[j, s])  for i = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios) for j = 1 : number_of_continuos_decision_variables]
            sorted_N1_indeces = sortperm(f_rank[1 : end], rev = true)[1 : N1]
            precision_p_iteration[sorted_N1_indeces] = precision_p_iteration[sorted_N1_indeces] .- 1
        end

    #print("the difference: $(UB - LB) \n \n " )
    #print("UB: $(UB) \n \n " )

    end

    return UB, iteration

end


number_of_scenarios = vcat( Array(2:10), Array(20:10:100) ) # number of scenarios for each of the instances
number_of_continuos_decision_variables = 10 # number of continuos variables per each secanrio
number_of_integer_decision_variables = 0 # number of integer variables per each scenario
number_of_constrains = 20 # number of the constraints per each scenario|
Qdensity = 1
precision_p = -18 .* ones(1, number_of_continuos_decision_variables)


dynamic_precision_RNMDT_algorithm(5, 3, 0.01, number_of_scenarios[10], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

original_problem = original_problem_generation(number_of_scenarios[10], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(original_problem)
objective_value(original_problem)


dynamic_RNMDT_problem, objective_Qs, objective_fs = dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios[2], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(dynamic_RNMDT_problem)
objective_value(dynamic_RNMDT_problem)


RNMDT_problem = RNMDT_problem_generation(-18, number_of_scenarios[2], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(RNMDT_problem)
objective_value(RNMDT_problem)
