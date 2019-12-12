using DataFrames, CSV, Plots

include("p_lagrangian_decomposition_Budnle.jl")

#intialization of the parameters for the bundle method
number_of_iterations = 20 # max number of iterations for the LAgrangian decomposition

p = -7 # precision factor for the mixed integer relaxation

number_of_scenarios = vcat( Array(2:10), Array(20:10:100) ) # number of scenarios for each of the instances
number_of_continuos_decision_variables = 10 # number of continuos variables per each secanrio
number_of_integer_decision_variables = 0 # number of integer variables per each scenario
number_of_constrains = 20 # number of the constraints per each scenario|
Qdensity = 1

continuous_instances = DataFrame(instance = Int[], original_model_solution = Float64[], original_model_time = Float64[], RNMDT_model_solution = Float64[], RNMDT_model_time = Float64[], p_RNMDT_model_solution = Float64[], p_RNMDT_model_time = Float64[] )

mixed_integer_instances = DataFrame(instance = Int[], RNMDT_model_solution = Float64[], RNMDT_model_time = Float64[], p_RNMDT_model_solution = Float64[], p_RNMDT_model_time = Float64[] )

for i = 1 : length(number_of_scenarios) # length(number_of_scenarios)

    original_problem = original_problem_generation(number_of_scenarios[i], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
    original_initial_time = time()
    optimize!(original_problem)
    original_final_time = time() - original_initial_time
    original_objective_value = objective_value(original_problem)

    RNMDT_problem = RNMDT_problem_generation(p, number_of_scenarios[i], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
    RNMDT_initial_time = time()
    optimize!(RNMDT_problem)
    RNMDT_final_time = time() - RNMDT_initial_time
    RNMDT_objective_value = objective_value(RNMDT_problem)

    result = Lagrangian_decomposition_bundle(p, number_of_scenarios[i], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

    push!(continuous_instances, (i, original_objective_value, original_final_time, RNMDT_objective_value, RNMDT_final_time, result[1][end], result[2]))

end

number_of_integer_decision_variables = 0 # number of integer variables per each scenario

for i = 1 : length(number_of_scenarios) # length(number_of_scenarios)

    RNMDT_problem = RNMDT_problem_generation(p, number_of_scenarios[i], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
    RNMDT_initial_time = time()
    optimize!(RNMDT_problem)
    RNMDT_final_time = time() - RNMDT_initial_time
    RNMDT_objective_value = objective_value(RNMDT_problem)

    result = Lagrangian_decomposition_bundle(p, number_of_scenarios[i], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

    push!( mixed_integer_instances, (i, RNMDT_objective_value, RNMDT_final_time, result[1][end], result[2]))

end

showall(continuous_instances)

showall(mixed_integer_instances)

CSV.write("continuos_instances.csv", continuous_instances)

CSV.write("mixed_integer_instances.csv", mixed_integer_instances)


original_problem = original_problem_generation(number_of_scenarios[18], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
original_initial_time = time()
optimize!(original_problem)
original_final_time = time() - original_initial_time
original_objective_value = objective_value(original_problem)

p = -10
RNMDT_problem = RNMDT_problem_generation(p, number_of_scenarios[1], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
RNMDT_initial_time = time()
optimize!(RNMDT_problem)
RNMDT_final_time = time() - RNMDT_initial_time
RNMDT_objective_value = objective_value(RNMDT_problem)
value.(RNMDT_problem[:x])





precision_p = p .* ones(1, number_of_continuos_decision_variables)
dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios[1], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
dynamic_RNMDT_initial_time = time()
optimize!(dynamic_RNMDT_problem)
dynamic_RNMDT_final_time = time() - dynamic_RNMDT_initial_time
dynamic_RNMDT_objective_value = objective_value(dynamic_RNMDT_problem)
value.(dynamic_RNMDT_problem[:x])

result = Lagrangian_decomposition_bundle(p, number_of_scenarios[18], number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

plot(Array(1:length(result[1])), result[1], xlabel = "iteration", ylabel = "upper bound", title = "Bundle method insipred lagrangian decomposition")
plot!(Array(1:length(result[1])), result[1], seriestype=:scatter)
savefig("Bundle method")
