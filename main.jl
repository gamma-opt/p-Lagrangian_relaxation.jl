using DataFrames, CSV, Plots

include("dynamic_precision_based_algorithms.jl")

number_of_scenarios = vcat( Array(2:10), Array(20:10:100) ) # number of scenarios for each of the instances
number_of_continuous_decision_variables = 10 # number of continuous variables per each secanrio
number_of_integer_decision_variables = 0 # number of integer variables per each scenario
number_of_constrains = 20 # number of the constraints per each scenario|
Qdensity = 1

precision_p = -1 .* ones(1, number_of_continuous_decision_variables)

max_number_of_iterations = 20
time_limit = 7200

# original problem
original_problem, objective_Qs, objective_fs = original_problem_generation(number_of_scenarios[1], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(original_problem)
objective_value(original_problem)

# RNMDT
#dynamic_precision_RNMDT_algorithm(5, 3, 0.1, time_limit, number_of_scenarios[1], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

# LD + RNDMT
#dynamic_precision_RNMDT_algorithm(5, 3, 0.1, time_limit, number_of_scenarios[2], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

# single iteration of LD + RNDMT
UB, xr, yr, wr = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios[1], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, 100)
