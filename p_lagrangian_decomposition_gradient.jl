using SharedArrays, Distributed
include("models_generation.jl")

# maximum number_of number_of_iterations
number_of_iterations = 300


# auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

# lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
[ vector_of_lambda_lagrangian[1, i] = -20.5 .+ 100.0 .* rand(1,
    number_of_continuous_decision_variables + number_of_integer_decision_variables)
        for i = 1 : number_of_scenarios - 1 ]

# indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
x_indices = 1 : number_of_continuous_decision_variables

# indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
y_indices = number_of_continuous_decision_variables + 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables

# vector that contains decision variables written in a column (x and y in this case)
# (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
decision_variables_values_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

# subgradient vector at each iteration
subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

result_decomposition = Array{Float64}(undef, 2, number_of_iterations) # first row - Upper bound, second row - Lower bound, the columns are number_of_iterations

#array of the models for every scenario with lambda inital
subproblems = subproblem


UB_decomposition = 0
LB_decomposition = 0

theta_decomposition = 0.0001

initial_time = time()

for iteration = 1:number_of_iterations - 1

     UB_decomposition = 0
     LB_decomposition = 0

        @sync @distributed for s = 1 : number_of_scenarios

        #objective_update
        @objective( subproblems[s], Max,
                (1 / number_of_scenarios) * ( sum( objective_Qs[s][i, i] * subproblems[s][:w_RNMDT][i, i]
                    for i = 1 : number_of_continuous_decision_variables )
                + 2 * sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                    for i = 1 : number_of_continuous_decision_variables,
                    j = i+1 : number_of_continuous_decision_variables)
                + sum( ( subproblems[s][:x] .* objective_fs[s][1, :] )
                    .+ ( subproblems[s][:y] .* objective_fs[s][2, :] ) )
                    + objective_fs[s][3, 1])
                +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[x_indices] .* subproblems[s][:x] )

                +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[y_indices] .* subproblems[s][:y] )
                )

            status = optimize!(subproblems[s])

            UB_decomposition = UB_decomposition + objective_value(subproblems[s])


            decision_variables_values_for_each_scenario[ x_indices, s ] = value.(subproblems[s][:x])
            decision_variables_values_for_each_scenario[ y_indices, s ] = value.(subproblems[s][:y])


            end
        LB_decomposition = 0.5 * UB_decomposition

        print(" iteration = $iteration, difference = $(UB_decomposition[1] - LB_decomposition[1]) \n ")
        result_decomposition[1, iteration] = UB_decomposition
        result_decomposition[2, iteration] = LB_decomposition


        [ subgradient_vector[iteration,  s - 1 ] =  decision_variables_values_for_each_scenario[:, 1] - decision_variables_values_for_each_scenario[ :, s] for  s in 2 : number_of_scenarios ]

        #[lambda_lagrangian[n,j,y,m] = lambda_lagrangian[n,j,y,m] - (sum(s_decomposition[n,j,y,:].^2)>0 ? theta_decomposition * (UB_decomposition[1] - LB_decomposition[1])/sum(s_decomposition[n,j,y,:].^2)*s_decomposition[n,j,y,m-1] : 0 ) for n in Nn, j in Jn, y in Yn, m in Mn[2:end] ]
        sum_of_the_squared_subgradeint_coordinates = sum( subgradient_vector[iteration, j] .^ 2 for j = 1:number_of_scenarios - 1 )
        [ vector_of_lambda_lagrangian[ iteration + 1, s] =
            vector_of_lambda_lagrangian[ iteration, s] .-
                [ sum_of_the_squared_subgradeint_coordinates[j] > 0 ? theta_decomposition * (UB_decomposition - LB_decomposition) * subgradient_vector[iteration, s][j] / sum_of_the_squared_subgradeint_coordinates[j] : 0
                    for j = 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables]'
                        for s = 1 : number_of_scenarios - 1]


        #[ lambda_lagrangian[n,j,y,m] = lambda_lagrangian[n,j,y,m] - ( sum( s_decomposition[n,j,y,:].^2 ) > 0 ? theta_decomposition * (UB_decomposition[1] - LB_decomposition[1]) / sum( s_decomposition[n,j,y,:] .^ 2 ) * s_decomposition[n,j,y,m] : 0 ) for n in Nn, j in Jn, y in Yn, m in 1:(size(Mn,1)-1) ]

        if mod(iteration,10) == 0
            global theta_decomposition = theta_decomposition*0.7
        end

        #theta_decomposition = theta_decomposition/(beta_decompostion*iteration)
    end
    final_time = time()-initial_time

function dynamic_precision_based_Lagrangian_decomposition_subgradient(precision_p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)

        # generating the parameters
        constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

        # limits for the initial values of the Lagrangian multipliers
        lagrangian_multipliers_max = 10
        lagrangian_multipliers_min = 0

        # auxiliary function for Lagrangian lagrangian_multipliers_representing_variableltipliers update
        f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1 : end]) : - lambda_lagrangian[dec_index-1])

        # lagrangian relaxation variables for the x and y non anticipativity conditions written in the column, for each iteration
        vector_of_lambda_lagrangian = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)
        [ vector_of_lambda_lagrangian[1, i] = lagrangian_multipliers_min .+ (lagrangian_multipliers_max - lagrangian_multipliers_min)  .* rand(1,
            number_of_continuous_decision_variables + number_of_integer_decision_variables)
                for i = 1 : number_of_scenarios - 1 ]

        # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
        x_indices = 1 : number_of_continuous_decision_variables

        # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
        y_indices = number_of_continuous_decision_variables + 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables

        # dual function at the lagragian multiplers' vector at correspondent iteration
        dual_objective_value_at_lagrangian = Array{Float64}(undef, 1, number_of_iterations)

        # vector that contains decision variables written in a column (x and y in this case)
        # (each row represnets the components of the correspondent lagrangian lagrangian_multipliers_representing_variableltiplier)
        decision_variables_values_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables + number_of_integer_decision_variables, number_of_scenarios)

        # the vector that contains quadraticity representing variables wij for each of the scenarios
        RNMDT_quadraticity_variables_w_for_each_scenario = SharedArray{Float64}(number_of_continuous_decision_variables, number_of_continuous_decision_variables, number_of_scenarios)

        # subgradient vector at each iteration
        subgradient_vector = Array{Any}(undef, number_of_iterations, number_of_scenarios - 1)

        #array of the subproblems for every scenario with lambda inital
        subproblems = dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios, number_of_continuous_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

        # upper bound for the original problem obtained by solving relaxed subproblems
        # and summing up the values of the objective functions
        ub_of_original_problem = SharedArray{Float64}(1,1)

        step_size = 0.3

        iteration = 0 #strating counter

        eps_stop = 1e-4 # stopping criteria

        initial_time = time()
        while (iteration < number_of_iterations) & ((iteration>6) ? (norm(dual_objective_value_at_lagrangian[iteration - 5:iteration-1].- dual_objective_value_at_lagrangian[iteration-6:iteration-2])^2 >= eps_stop) : true)

            iteration += 1

            ub_of_original_problem[1] = 0

            @sync @distributed for s in 1 : number_of_scenarios

                #objective_update
                @objective( subproblems[s], Max,
                    (1/number_of_scenarios) *
                        (
                        sum( objective_Qs[s][i, i] * subproblems[s][:w_RNMDT][i, i]
                            for i = 1 : number_of_continuous_decision_variables )
                        + 2 * sum(objective_Qs[s][i, j] * subproblems[s][:w_RNMDT][i, j]
                            for i = 1 : number_of_continuous_decision_variables,
                                j = i+1 : number_of_continuous_decision_variables)
                        + sum( subproblems[s][:x][i] * objective_fs[s][1, i]  for i = 1:number_of_continuous_decision_variables)
                        + sum( subproblems[s][:y][j] * objective_fs[s][2, j]  for j = 1:number_of_integer_decision_variables)
                        )
                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[x_indices] .* subproblems[s][:x] )

                        +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[iteration, :], s )[y_indices] .* subproblems[s][:y] )
                )

                status = optimize!(subproblems[s])
                obj_value = objective_value(subproblems[s])
                ub_of_original_problem[1] = ub_of_original_problem[1] + obj_value

                decision_variables_values_for_each_scenario[ x_indices, s ] = value.(subproblems[s][:x])
                decision_variables_values_for_each_scenario[ y_indices, s ] = value.(subproblems[s][:y])
                RNMDT_quadraticity_variables_w_for_each_scenario[ :, :, s] = value.(subproblems[s][:w_RNMDT])
            end

            dual_objective_value_at_lagrangian[iteration] = ub_of_original_problem[1]

            [ subgradient_vector[iteration,  s - 1] =  decision_variables_values_for_each_scenario[ :, 1] - decision_variables_values_for_each_scenario[ :, s] for  s in 2: number_of_scenarios ]

            if iteration < number_of_iterations
                sum_of_the_squared_subgradeint_coordinates = sum( subgradient_vector[iteration, j] .^ 2 for j = 1:number_of_scenarios - 1 )
                [ vector_of_lambda_lagrangian[ iteration + 1, s] =
                    vector_of_lambda_lagrangian[ iteration, s] .-
                        [ sum_of_the_squared_subgradeint_coordinates[j] > 0 ? step_size * (0.7*dual_objective_value_at_lagrangian[iteration]) * subgradient_vector[iteration, s][j] / sum_of_the_squared_subgradeint_coordinates[j] : 0
                            for j = 1 : number_of_continuous_decision_variables + number_of_integer_decision_variables]'
                                for s = 1 : number_of_scenarios - 1]
            end

            print("iteration  =  $iteration \n \n")

            print("lag. mul. values  = \n")
            print(vector_of_lambda_lagrangian[iteration, :] )
            print("\n \n")

            print("Optimal value of dual objective with fixed lag. mul. values: $(dual_objective_value_at_lagrangian[iteration]) \n \n")

            print("decision_variables: \n")
            print(decision_variables_values_for_each_scenario[x_indices,:])
            print("\n \n")
            if mod(iteration,10) == 0
                step_size = step_size*0.7
            end

    end

    return dual_objective_value_at_lagrangian, decision_variables_values_for_each_scenario[x_indices, :], decision_variables_values_for_each_scenario[y_indices, :], RNMDT_quadraticity_variables_w_for_each_scenario, vector_of_lambda_lagrangian[ 1 : iteration, :]

end
