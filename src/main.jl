
using DataFrames, JuMP, CSV, Plots, Gurobi, Random, XLSX, Ipopt, Suppressor, LinearAlgebra, SparseArrays, SharedArrays, Distributed, Dates

include("parameters_generation.jl")
include("models_generation.jl")
include("p_lagrangian_decomposition_bundle.jl")
include("dynamic_precision_based_algorithms.jl")

number_of_scenarios = Array(5 : 5 : 25)  # number of scenarios for each of the instances
vector_number_of_continuous_decision_variables = Array(20 : 10 : 40) # number of continuous variables per each secanrio
vector_number_of_integer_decision_variables = Array(5 : 5 : 15) # number of integer variables per each scenario
vector_number_of_constraints = Array(45 : 15 : 75 ) # number of the constraints per each scenario|
global Qdensity = 0.8



#-------------------------------------------------------------------------------

#global N1 = 5
global N2 = 5
global tolerance = 0.001
global time_limit = 3600
global max_number_of_iterations = 100

global link = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_Decomposition-/Results/experiments_5.03.20_1try/"
global loglink_par = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_Decomposition-/Results/experiments_5.03.20_1try/logfiles_non_par/"
global loglink_par_bundle = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_Decomposition-/Results/experiments_5.03.20_1try/logfiles_non_par_bundle/"



for seed  = 0 : 1 : 0
    for instance = 3:3
        for instance_size = 2:2

            number_of_continuous_decision_variables = vector_number_of_continuous_decision_variables[instance_size]
            number_of_integer_decision_variables = vector_number_of_integer_decision_variables[instance_size]
            number_of_constraints = vector_number_of_constraints[instance_size]

            #original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no", time_limit, seed)
            #optimize!(original_problem)
            #df = DataFrame( optimal_solution = objective_value(original_problem), gap = MOI.get(original_problem, MOI.RelativeGap()) )
            #XLSX.writetable(link * "original_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            global N1 = Int(round( 0.1*(number_of_continuous_decision_variables)*number_of_scenarios[instance], digits = 0 ) )

            #dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)
            #df = DataFrame(iter = collect(1:size(dynamic_RNMDT_results, 2)), UB = dynamic_RNMDT_results[1, :], LB = dynamic_RNMDT_results[2, :], gap =  dynamic_RNMDT_results[3, :],  time = dynamic_RNMDT_results[4, :])
            #XLSX.writetable(link * "dynamic_RNDMT_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            dynamic_LD_RNMDT_results_non_par = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "non_parallelised", seed)
            df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results_non_par, 2)), UB = dynamic_LD_RNMDT_results_non_par[1, :], LB = dynamic_LD_RNMDT_results_non_par[2, :], gap  = dynamic_LD_RNMDT_results_non_par[3, :], num_iter = dynamic_LD_RNMDT_results_non_par[4, :], act_step = dynamic_LD_RNMDT_results_non_par[5, :], time = dynamic_LD_RNMDT_results_non_par[6, :])
            XLSX.writetable(link * "dynamic_LD_RNMDT_not_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_eigth.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            #dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm( N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "parallelised", seed )
            #df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results, 2)), UB = dynamic_LD_RNMDT_results[1, :], LB = dynamic_LD_RNMDT_results[2, :], gap  = dynamic_LD_RNMDT_results[3, :], num_iter = dynamic_LD_RNMDT_results[4, :], act_step = dynamic_LD_RNMDT_results[5, :], time = dynamic_LD_RNMDT_results[6, :])
            #XLSX.writetable(link * "dynamic_LD_RNMDT_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_fiftingth_fixed_bundle_subproblem_to_1_thread.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

        end
    end
end

global loglink_par = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_Decomposition-/Results/experiments_5.03.20_1try/logfiles_par/"
global loglink_par_bundle = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_Decomposition-/Results/experiments_5.03.20_1try/logfiles_par_bundle/"

for seed  = 0 : 1 : 0
    for instance = 3:3
        for instance_size = 2:2

            number_of_continuous_decision_variables = vector_number_of_continuous_decision_variables[instance_size]
            number_of_integer_decision_variables = vector_number_of_integer_decision_variables[instance_size]
            number_of_constraints = vector_number_of_constraints[instance_size]

            #original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no", time_limit, seed)
            #optimize!(original_problem)
            #df = DataFrame( optimal_solution = objective_value(original_problem), gap = MOI.get(original_problem, MOI.RelativeGap()) )
            #XLSX.writetable(link * "original_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            global N1 = Int(round( 0.1*(number_of_continuous_decision_variables)*number_of_scenarios[instance], digits = 0 ) )

            #dynamic_RNMDT_results = dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)
            #df = DataFrame(iter = collect(1:size(dynamic_RNMDT_results, 2)), UB = dynamic_RNMDT_results[1, :], LB = dynamic_RNMDT_results[2, :], gap =  dynamic_RNMDT_results[3, :],  time = dynamic_RNMDT_results[4, :])
            #XLSX.writetable(link * "dynamic_RNDMT_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            #dynamic_LD_RNMDT_results_non_par = dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "non_parallelised", seed)
            #df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results_non_par, 2)), UB = dynamic_LD_RNMDT_results_non_par[1, :], LB = dynamic_LD_RNMDT_results_non_par[2, :], gap  = dynamic_LD_RNMDT_results_non_par[3, :], num_iter = dynamic_LD_RNMDT_results_non_par[4, :], act_step = dynamic_LD_RNMDT_results_non_par[5, :], time = dynamic_LD_RNMDT_results_non_par[6, :])
            #XLSX.writetable(link * "dynamic_LD_RNMDT_not_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_eigth.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            dynamic_LD_RNMDT_results = dynamic_precision_LD_RNMDT_algorithm( N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "parallelised", seed )
            df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results, 2)), UB = dynamic_LD_RNMDT_results[1, :], LB = dynamic_LD_RNMDT_results[2, :], gap  = dynamic_LD_RNMDT_results[3, :], num_iter = dynamic_LD_RNMDT_results[4, :], act_step = dynamic_LD_RNMDT_results[5, :], time = dynamic_LD_RNMDT_results[6, :])
            XLSX.writetable(link * "dynamic_LD_RNMDT_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_fiftingth_fixed_bundle_subproblem_to_1_thread.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

        end
    end
end
