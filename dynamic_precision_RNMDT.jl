include("Models_generation.jl")

#dynamic_precision_RNMDT_algorithm(N1, N2, precision_p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

#end

N1 = 5
N2 = 3

# initialisation of the parameters
precision_p_inital = -10 .* ones(1, number_of_continuos_decision_variables)
iteration = 0
UB = Inf

original_problem = original_problem_generation(number_of_scenarios[1], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(original_problem)
original_problem_objective_value = objective_value(original_problem)

RNMDT_problem = RNMDT_problem_generation(-20, number_of_scenarios[1], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(RNMDT_problem)
objective_value(RNMDT_problem)
xr = value.(RNMDT_problem[:x])
UB = sum( (1/number_of_scenarios[1]) *
    (
        sum( xr[i, s] * objective_Qs[s][i, j] * xr[j, s] for i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables)
        + sum( xr[i, s] * objective_fs[s][1, i] for i = 1 : number_of_continuos_decision_variables)
        + ( number_of_integer_decision_variables == 0 ? 0 : sum( yr[i, s] * objective_fs[s][2, i] for i = 1 : number_of_integer_decision_variables) )
    )
    for s in 1:number_of_scenarios[1])



dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p_inital, number_of_scenarios[1], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(dynamic_RNMDT_problem)
LB = objective_value(dynamic_RNMDT_problem)
xr = value.(dynamic_RNMDT_problem[:x])
yr = value.(dynamic_RNMDT_problem[:y])
w_RNMDT_r = value.(dynamic_RNMDT_problem[:w_RNMDT])

UB = sum( (1/number_of_scenarios[1]) *
    (
        sum( xr[i, s] * objective_Qs[s][i, j] * xr[j, s] for i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables)
        + sum( xr[i, s] * objective_fs[s][1, i] for i = 1 : number_of_continuos_decision_variables)
        + ( number_of_integer_decision_variables == 0 ? 0 : sum( yr[i, s] * objective_fs[s][2, i] for i = 1 : number_of_integer_decision_variables) )
    )
    for s in 1:number_of_scenarios[1])
