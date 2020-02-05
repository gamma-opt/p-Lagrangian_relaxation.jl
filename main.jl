
using DataFrames, CSV, Plots, Gurobi, Random, XLSX

include("dynamic_precision_based_algorithms.jl")

number_of_scenarios = vcat( Array(2:10), Array(20:10:100) ) # number of scenarios for each of the instances
number_of_continuous_decision_variables = 10 # number of continuous variables per each secanrio
number_of_integer_decision_variables = 4 # number of integer variables per each scenario
number_of_constrains = 20 # number of the constraints per each scenario|
Qdensity = 1.0

instance = 2
#---------------------------original problem------------------------------------
original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
JuMP.optimize!(original_problem)
objective_value(original_problem)
value.(original_problem[:y])
value.(original_problem[:x])
has_duals(original_problem)

#might not work for some of the solvers
#dual_objective_value(original_problem)
#-----------------------precision-based RNMDT-----------------------------------

# single model generation and optimisation with fixed p
p = -1
precision_p = p .* ones(1, number_of_continuous_decision_variables)
dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
JuMP.optimize!(dynamic_RNMDT_problem)
objective_value(dynamic_RNMDT_problem)

value.(dynamic_RNMDT_problem[:x])
value.(dynamic_RNMDT_problem[:y])

# dynamic precision algorithm test
N1 = 8
N2 = 3
tolerance = 0.0001
time_limit = 7200
max_number_of_iterations = 100

dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)

#original RNMDT
p = - 3
orginial_RNDMT_problem = RNMDT_problem_generation(p, number_of_scenarios[2], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)
optimize!(orginial_RNDMT_problem)
objective_value(orginial_RNDMT_problem)
value.(orginial_RNDMT_problem[:x])

#-----------------------dynamic precision-based LD + RNDMT----------------------
subproblem = dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios[1], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
optimize!(subproblem[1])
objective_value(subproblem[1])
objective_value(subproblem[2])

# single iteration of LD + RNDMT
p = -1
precision_p = p .* ones(1, number_of_continuous_decision_variables)
max_number_of_iterations = 100
UB, xr, yr, wr, lambda = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity, max_number_of_iterations)

# dynamic precision algorithm test
N1 = 8
N2 = 3
tolerance = 0.00001
time_limit = 7200
max_number_of_iterations = 100

dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity )

#attempt for the subgradient
#dual_objective_value_at_lagrangian, decision_variables_values_for_each_scenario, decision_variables_values_for_each_scenario, RNMDT_quadraticity_variables_w_for_each_scenario, vector_of_lambda_lagrangian = dynamic_precision_based_Lagrangian_decomposition_subgradient(precision_p, number_of_scenarios[2], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, max_number_of_iterations)
#plt = plot(collect(1:max_number_of_iterations), dual_objective_value_at_lagrangian')
#savefig(plt, "subgradient_attemp.pdf")

#-------------------------attempt for the plots---------------------------------

#using Plots

plot1 = zeros(1, size(lambda,1))
[plot1[1, i] =  lambda[i,1][1] for i = 1 : size(lambda,1)]
pl = plot(collect(1:size(lambda,1)), plot1', label = "lambda 1")

plot2 = zeros(1, size(lambda,1))
[plot2[1, i] =  lambda[i,2][1] for i = 1 : size(lambda,1)]
plot!(collect(1:size(lambda,1)), plot2', label = "lambda 2")

plot3 = zeros(1, size(lambda,1))
[plot3[1, i] =  lambda[i,3][1] for i = 1 : size(lambda,1)]
plot!(collect(1:size(lambda,1)), plot3', label = "Lambda")

title!("m = 0.1, d = 0.00005")
xaxis!("iteration")
yaxis!("value of lagrangian multiplier")

savefig(pl, "LD_RNMDT_1_variable_4_scenarios")

#-------------------------------------------------------------------------------

p = -3
precision_p = p .* ones(1, number_of_continuous_decision_variables)
df_single_iteration = DataFrame(instance = Int64[], Orig = Float64[], RNMDT = Float64[], LD_RNDMT = Float64[], num_ser_steps = Int64[])

for instance = 13:18

    original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
    JuMP.optimize!(original_problem)

    dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)
    JuMP.optimize!(dynamic_RNMDT_problem)

    UB, xr, yr, wr, lambda, number_of_serious_steps = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity, max_number_of_iterations)

    push!(df_single_iteration, (instance, objective_value(original_problem), objective_value(dynamic_RNMDT_problem), UB, number_of_serious_steps))

end

XLSX.writetable("single_iteration.xlsx", collect(DataFrames.eachcol(df_single_iteration)), DataFrames.names(df_single_iteration))

#-------------------------------------------------------------------------------

N1 = 8
N2 = 3
tolerance = 0.00001
time_limit = 7200
max_number_of_iterations = 100

df_LD_RNMDT = DataFrame(instance  = Int64[], iteration = Int64[], UB = Float64[], LB = Float64[], Gap  = Float64[], time = Float64[])
df_RNMDT = DataFrame(instance  = Int64[], iteration = Int64[], UB = Float64[], LB = Float64[], Gap  = Float64[], time = Float64[])

for instance = 1:18

    dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity)

    dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constrains, Qdensity )

    for i = 1:size(dynamic_LD_RNMDT_results[1],2)
        push!(df_LD_RNMDT, (instance, i, dynamic_LD_RNMDT_results[1][1, i], dynamic_LD_RNMDT_results[1][2, i], dynamic_LD_RNMDT_results[1][3, i], dynamic_LD_RNMDT_results[2]))
    end

    for i = 1:size(dynamic_RNMDT_results[1],2)
        push!(df_RNMDT, (instance, i, dynamic_RNMDT_results[1][1, i], dynamic_RNMDT_results[1][2, i], dynamic_RNMDT_results[1][3, i], dynamic_RNMDT_results[2]))
    end

end

XLSX.writetable("LD_RNMDT.xlsx", collect(DataFrames.eachcol(df_LD_RNMDT)), DataFrames.names(df_LD_RNMDT))
XLSX.writetable("RNMDT.xlsx", collect(DataFrames.eachcol(df_RNMDT)), DataFrames.names(df_RNMDT))
