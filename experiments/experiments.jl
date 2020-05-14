
using JuMP, DataFrames, XLSX, Dates, p_Lagrangian_relaxation


number_of_scenarios = Array(5 : 5 : 25)  # number of scenarios for each of the instances
vector_number_of_continuous_decision_variables = Array(20 : 10 : 40) # number of continuous variables per each secanrio
vector_number_of_integer_decision_variables = Array(5 : 5 : 15) # number of integer variables per each scenario
vector_number_of_constraints = Array(45 : 15 : 75 ) # number of the constraints per each scenario|
global Qdensity = 0.8



#-------------------------------------------------------------------------------

# defining the parameters for the dynamic precision-based algorithm
global N2 = 5
global tolerance = 0.001
global time_limit = 200
global max_number_of_iterations = 100

# the parameters for the Bundle method are defined inside the dp_Lagrangian_relaxation.dynamic_precision_LD_RNMDT_algorithm function

# write here the link where the experiments results will be stored
global link = "/Users/nikitabelyak/Dropbox (Aalto)/p_Lagrangian_relaxation/experiments/"

for seed  = 0 : 1 : 4 # random seed
    for instance = 1:5 # number of scnearios: 5,10,15,20,25
        for instance_size = 1:3 #instance size: small, medium and large

            number_of_continuous_decision_variables = vector_number_of_continuous_decision_variables[instance_size]
            number_of_integer_decision_variables = vector_number_of_integer_decision_variables[instance_size]
            number_of_constraints = vector_number_of_constraints[instance_size]

            # Full space MIQCQP
            original_problem, objective_Qs, objective_fs, constraint_Qs, constraint_fs, objective_c = p_Lagrangian_relaxation.original_problem_generation(number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "no", time_limit, seed)
            optimize!(original_problem)
            df = DataFrame( optimal_solution = objective_value(original_problem), gap = MOI.get(original_problem, MOI.RelativeGap()) )
            XLSX.writetable(link * "original_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            # calculating N1 as 10% of the number of continuous decision variables
            global N1 = Int(round( 0.1*(number_of_continuous_decision_variables)*number_of_scenarios[instance], digits = 0 ) )

            # RNMDT relaxation solved using the dynamic precision-based
            dynamic_RNMDT_results = p_Lagrangian_relaxation.dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)
            df = DataFrame(iter = collect(1:size(dynamic_RNMDT_results, 2)), UB = dynamic_RNMDT_results[1, :], LB = dynamic_RNMDT_results[2, :], gap =  dynamic_RNMDT_results[3, :],  time = dynamic_RNMDT_results[4, :])
            XLSX.writetable(link * "dynamic_RNDMT_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_seventh.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            # p-Lagrangian relaxation solved using the dynamic precision-based (sequential version)
            dynamic_LD_RNMDT_results_non_par = p_Lagrangian_relaxation.dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "non_parallelised", seed)
            df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results_non_par, 2)), UB = dynamic_LD_RNMDT_results_non_par[1, :], LB = dynamic_LD_RNMDT_results_non_par[2, :], gap  = dynamic_LD_RNMDT_results_non_par[3, :], num_iter = dynamic_LD_RNMDT_results_non_par[4, :], act_step = dynamic_LD_RNMDT_results_non_par[5, :], time = dynamic_LD_RNMDT_results_non_par[6, :])
            XLSX.writetable(link * "dynamic_LD_RNMDT_not_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_eigth.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

            # p-Lagrangian relaxation solved using the dynamic precision-based (parallelised version)
            dynamic_LD_RNMDT_results = p_Lagrangian_relaxation.dynamic_precision_LD_RNMDT_algorithm( N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios[instance], number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, "parallelised", seed )
            df = DataFrame(iter = collect(1:size(dynamic_LD_RNMDT_results, 2)), UB = dynamic_LD_RNMDT_results[1, :], LB = dynamic_LD_RNMDT_results[2, :], gap  = dynamic_LD_RNMDT_results[3, :], num_iter = dynamic_LD_RNMDT_results[4, :], act_step = dynamic_LD_RNMDT_results[5, :], time = dynamic_LD_RNMDT_results[6, :])
            XLSX.writetable(link * "dynamic_LD_RNMDT_par_instance_instance_$(number_of_scenarios[instance])_scen_$(number_of_continuous_decision_variables)_cont_var_$(number_of_integer_decision_variables)_int_var_$(number_of_constraints)_constraints_($seed)_seed_$(Dates.today())_fiftingth_fixed_bundle_subproblem_to_1_thread.xlsx", collect(DataFrames.eachcol(df)), DataFrames.names(df))

        end
    end
end
