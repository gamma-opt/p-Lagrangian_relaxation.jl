
using DataFrames, CSV, Plots, Gurobi, Random, XLSX, Ipopt, Suppressor, LinearAlgebra, SparseArrays

include("dynamic_precision_based_algorithms.jl")
number_of_scenarios = vcat( Array(2:10), Array(20:10:100) ) # number of scenarios for each of the instances
number_of_continuous_decision_variables = 10 # number of continuous variables per each secanrio
number_of_integer_decision_variables = 4 # number of integer variables per each scenario
number_of_constraints = 20 # number of the constraints per each scenario|
Qdensity = 0.8

instance = 18

macro warmup(number_of_workers)
        @sync @distributed for i = 1 : number_of_workers
            warmup_model = Model(with_optimizer(Gurobi.Optimizer))
            @variable(warmup_model, x )
            @constraint(warmup_model, 0 <= x <= 10)
            @objective(warmup_model, Max, x^2)
            optimize!(warmup_model)
        end

end
#---------------------------original problem------------------------------------
original_problem1, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[7], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no")

JuMP.optimize!(original_problem1)
objective_value(original_problem1)

original_problem2, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[8], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no")

JuMP.optimize!(original_problem2)
objective_value(original_problem2)


original_problem3, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[9], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no")

JuMP.optimize!(original_problem3)
objective_value(original_problem3)

value.(original_problem[:y])
value.(original_problem[:x])
has_duals(original_problem)

#might not work for some of the solvers
#dual_objective_value(original_problem)
#-----------------------precision-based RNMDT-----------------------------------

# single model generation and optimisation with fixed p
for p  = -1 : -1 : -1
#p = -4
    precision_p = p .* ones(number_of_continuous_decision_variables, number_of_scenarios[instance])
    dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p,  number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
    #stime = time()
    JuMP.optimize!(dynamic_RNMDT_problem)
    #etime = time() - stime
    print("$(objective_value(dynamic_RNMDT_problem))\n\n")
end

value.(dynamic_RNMDT_problem[:x])
value.(dynamic_RNMDT_problem[:y])

# dynamic precision algorithm test
N1 = 10
N2 = 3
tolerance = 0.001
time_limit = 7200
max_number_of_iterations = 15

stime = time()
dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, 100, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
etime = time() - stime

df = DataFrame(UB = dynamic_RNMDT_results[1][1, :], LB = dynamic_RNMDT_results[1][2, :], gap =  dynamic_RNMDT_results[1][3, :] )
XLSX.writetable("instance_$instance.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))
#original RNMDT
p = - 3
orginial_RNDMT_problem = RNMDT_problem_generation(p, number_of_scenarios[2], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constraints, Qdensity)
optimize!(orginial_RNDMT_problem)
objective_value(orginial_RNDMT_problem)
value.(orginial_RNDMT_problem[:x])

#-----------------------dynamic precision-based LD + RNDMT----------------------
subproblem = dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
optimize!(subproblem[3])
objective_value(subproblem[1])
objective_value(subproblem[3])

# single iteration of LD + RNDMT
rmprocs(workers()) #removing all existing workers
addprocs(6) #add the number of workers equal to the number of scenarios (otherwise there might be done redundant work by means of copying the array of models to every process and updating objective everywhere)
@everywhere using JuMP, Gurobi, Random, SparseArrays, LinearAlgebra
@warmup(6)
@everywhere include("models_generation.jl")

p = -
precision_p = p .* ones(number_of_continuous_decision_variables, number_of_scenarios[instance])
max_number_of_iterations = 100

# generating the inital values for the center of gravity
center_of_gravity_min = 0
center_of_gravity_max = 10
Random.seed!(0)
center_of_gravity_inital_value = Array{Any}(undef, number_of_scenarios[instance] - 1)
[ center_of_gravity_inital_value[i] = center_of_gravity_min .+ (center_of_gravity_max - center_of_gravity_min)  .* rand(1,
    number_of_integer_decision_variables)
        for i = 1 : number_of_scenarios[instance] - 1 ]

init = time()
UB, xr, yr, wr, lambda, number_of_the_serious_steps, center = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations, center_of_gravity_inital_value, 1)
final = time() - init
df = DataFrame(iter = collect(1:length(UB)), UB = UB )
XLSX.writetable("instance_$instance SS_1000 ES 10.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

rmprocs(workers()) #removing all existing workers
addprocs(3) #add the number of workers equal to the number of scenarios (otherwise there might be done redundant work by means of copying the array of models to every process and updating objective everywhere)
@everywhere using JuMP, Gurobi
@warmup(3)
@everywhere include("models_generation.jl")

p = -4
precision_p = p .* ones( number_of_continuous_decision_variables, 100)
max_number_of_iterations = 100

# generating the inital values for the center of gravity
center_of_gravity_min = 0
center_of_gravity_max = 10
Random.seed!(0)
center_of_gravity_inital_value = Array{Any}(undef, number_of_scenarios[instance] - 1)
[ center_of_gravity_inital_value[i] = center_of_gravity_min .+ (center_of_gravity_max - center_of_gravity_min)  .* rand(1,
    number_of_integer_decision_variables)
        for i = 1 : number_of_scenarios[instance] - 1 ]

center_of_gravity_inital_value = center1
init1 = time()
UB1, xr1, yr1, wr1, lambda1, number_of_the_serious_steps1, center1= dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, 100, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations, center_of_gravity_inital_value, 1)
final1 = time() - init1

original_problem, objective_Qs, objective_fs, objective_c, constraint_Qs = original_problem_generation(100, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "yes")

xr = repeat( Int.(round.( sum(xr1, dims = 2)./100, digits = 0)), 1, 100)
xr = repeat([50, 0 , 39, 16], 1, 100)
# setting strating values for the continuous variables
#set_start_value.(original_problem[:y], yr1)

#fixing the values for the integer variables
fix.(original_problem[:x], xr)
#set_start_value.(original_problem[:x], xr)

optimize!(original_problem)
LB  = objective_value(original_problem)

# dynamic precision algorithm test
N1 = 100
N2 = 3
tolerance = 0.001
time_limit = 7200
max_number_of_iterations = 15

rmprocs(workers()) #removing all existing workers
addprocs(number_of_scenarios[instance]) #add the number of workers equal to the number of scenarios (otherwise there might be done redundant work by means of copying the array of models to every process and updating objective everywhere)
@everywhere using JuMP, Gurobi
@everywhere include("dynamic_precision_based_algorithms.jl")

stime = time()
dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm_random_center(N1, N2, tolerance, time_limit, max_number_of_iterations, 100, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity )
etime = time() - stime
df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results[1],2)), UB = dynamic_LD_RNMDT_results[1][1, :], LB = dynamic_LD_RNMDT_results[1][2, :], gap  = dynamic_LD_RNMDT_results[1][3, :], num_iter = dynamic_LD_RNMDT_results[1][4, :], act_step = dynamic_LD_RNMDT_results[1][5, :])
XLSX.writetable("instance_$instance _rc.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

stime = time()
dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, 100, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity )
etime = time() - stime
df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results[1],2)), UB = dynamic_LD_RNMDT_results[1][1, :], LB = dynamic_LD_RNMDT_results[1][2, :], gap  = dynamic_LD_RNMDT_results[1][3, :], num_iter = dynamic_LD_RNMDT_results[1][4, :], act_step = dynamic_LD_RNMDT_results[1][5, :])
XLSX.writetable("instance_$instance _sc.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

#attempt for the subgradient
#dual_objective_value_at_lagrangian, decision_variables_values_for_each_scenario, decision_variables_values_for_each_scenario, RNMDT_quadraticity_variables_w_for_each_scenario, vector_of_lambda_lagrangian = dynamic_precision_based_Lagrangian_decomposition_subgradient(precision_p, number_of_scenarios[2], number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations)
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

    original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
    JuMP.optimize!(original_problem)

    dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
    JuMP.optimize!(dynamic_RNMDT_problem)

    UB, xr, yr, wr, lambda, number_of_serious_steps = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations)

    push!(df_single_iteration, (instance, objective_value(original_problem), objective_value(dynamic_RNMDT_problem), UB, number_of_serious_steps))

end

XLSX.writetable("single_iteration.xlsx", collect(DataFrames.eachcol(df_single_iteration)), DataFrames.names(df_single_iteration))

#-------------------------------------------------------------------------------

N1 = 5
N2 = 10
tolerance = 0.0001
time_limit = 7200
max_number_of_iterations = 100

df_LD_RNMDT = DataFrame(instance  = Int64[], iteration = Int64[], UB = Float64[], LB = Float64[], Gap  = Float64[], time = Float64[])
df_RNMDT = DataFrame(instance  = Int64[], iteration = Int64[], UB = Float64[], LB = Float64[], Gap  = Float64[], time = Float64[])

for instance = 14:18

    dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)

    dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity )

    for i = 1:size(dynamic_LD_RNMDT_results[1],2)
        push!(df_LD_RNMDT, (instance, i, dynamic_LD_RNMDT_results[1][1, i], dynamic_LD_RNMDT_results[1][2, i], dynamic_LD_RNMDT_results[1][3, i], dynamic_LD_RNMDT_results[2]))
    end

    for i = 1:size(dynamic_RNMDT_results[1],2)
        push!(df_RNMDT, (instance, i, dynamic_RNMDT_results[1][1, i], dynamic_RNMDT_results[1][2, i], dynamic_RNMDT_results[1][3, i], dynamic_RNMDT_results[2]))
    end

end

XLSX.writetable("LD_RNMDT.xlsx", collect(DataFrames.eachcol(df_LD_RNMDT)), DataFrames.names(df_LD_RNMDT))
XLSX.writetable("RNMDT.xlsx", collect(DataFrames.eachcol(df_RNMDT)), DataFrames.names(df_RNMDT))

#-------------------------------------------------------------------------------

df_RNMDT_time = DataFrame(number_of_scenarios  = Int64[], objective_value = Float64[], time = Float64[])
df_LD_RNMDT_time = DataFrame(number_of_scenarios  = Int64[], objective_value = Float64[], num_iter = Int64[], time = Float64[])

center_of_gravity_min = 0
center_of_gravity_max = 10

max_number_of_iterations = 100

p = -1
precision_p = p .* ones(1, number_of_continuous_decision_variables)

for i = 600:100:1000

    #dynamic_RNMDT_problem = dynamic_precision_RNMDT_problem_generation(precision_p, i, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity)
    #stime = time()
    #JuMP.optimize!(dynamic_RNMDT_problem)
    #etime = time() - stime

    #push!(df_RNMDT_time, (i, objective_value(dynamic_RNMDT_problem), etime))

    Random.seed!(0)
    center_of_gravity_inital_value = Array{Any}(undef, i - 1)
    [ center_of_gravity_inital_value[i] = center_of_gravity_min .+ (center_of_gravity_max - center_of_gravity_min)  .* rand(1,
        number_of_integer_decision_variables)
            for i = 1 : i - 1 ]

    stime = time()
    UB, xr, yr, wr, lambda, number_of_serious_steps = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, i, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations, center_of_gravity_inital_value)
    etime = time() - stime
    push!(df_LD_RNMDT_time, (i, UB, size(lambda,1), etime))

end

XLSX.writetable("time_LD_RNMDT2.xlsx", collect(DataFrames.eachcol(df_LD_RNMDT_time)), DataFrames.names(df_LD_RNMDT_time))
#XLSX.writetable("time_RNMDT.xlsx", collect(DataFrames.eachcol(df_RNMDT_time)), DataFrames.names(df_RNMDT_time))

Random.seed!(0)
center_of_gravity_inital_value = Array{Any}(undef, 300- 1)
[ center_of_gravity_inital_value[i] = center_of_gravity_min .+ (center_of_gravity_max - center_of_gravity_min)  .* rand(1,
    number_of_integer_decision_variables)
        for i = 1 : 300 - 1 ]


rmprocs(workers()) #removing all existing workers
addprocs(6) #add the number of workers equal to the number of scenarios (otherwise there might be done redundant work by means of copying the array of models to every process and updating objective everywhere)
@everywhere using JuMP, Gurobi, Random, SparseArrays, LinearAlgebra
@warmup(6)
@everywhere include("models_generation.jl")

stime = time()
UB, xr, yr, wr, lambda, number_of_serious_steps = dynamic_precision_based_Lagrangian_decomposition_bundle(precision_p, 300, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, max_number_of_iterations, center_of_gravity_inital_value)
etime = time() - stime
