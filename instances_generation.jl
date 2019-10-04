using LinearAlgebra, SparseArrays, JuMP, Gurobi,  Random

number_of_scenarios = 3
number_of_continuos_decision_variables = 4
number_of_integer_decision_variables = number_of_continuos_decision_variables
number_of_constrains = 2 # number of the constraints for each of the scenario
Qdensity = 0.6 # dencity of the matrices
Max_value_for_matrix_elements = 5
#---------------generating constraints for the scenarios------------------------

# auxiliary function for generating quadratic matrices for cosnstraints and obejctive with predefined densiity
function matrix_generation(density, dimention, max_range, PSD)

    # creating the matrix of predefined density,
    # taking into account that this density is defined for the case when a matrix is symmetrical
    # we generate triangular matrix in such a way that if it is represented as a symmetric one
    # (dividing upper triangular part by 2 and reflecting it towards the diagonal)
    # it will have that predefined density

    # separately generating diagonal elements not equal to zero (to be able to create PSD after)
    diagonal_elements  = (max_range) .* round.(rand(1, dimention), digits = 1)

    # calculating the number of non-zero elements in the upper diagonal part
    # depending on the predefined density
    # and taking into account that diagonal elements are non-zero
    number_of_other_non_zero_elements = Int(round((density * dimention^2 - dimention) / 2))

    # creating one-dimensional array containing an abovementioned number of non zero elements,
    # filling the rest with zeros and shuffling the array
    upper_triangular_elements  = shuffle!([(max_range .* round.(rand(1, number_of_other_non_zero_elements), digits = 1 )) zeros(1, Int(dimention*(dimention-1)/2) - number_of_other_non_zero_elements) ] )

    # creating final matrix with zero elements
    Q = zeros(dimention, dimention)

    # using above mentioned diagonal elements and upper triangular part
    # rewriting zeros in the upper triangular part accordingly
    for i = 1:dimention
        for j = i:dimention
            if i == j # if the element is on the diagonal using an array of diagonal elements
                Q[i,i] = diagonal_elements[i]
            else
                Q[i,j] = upper_triangular_elements[1]/2 # otherwise using one elemnts from upper triangular elements
                Q[j,i] = upper_triangular_elements[1]/2 # otherwise using one elemnts from upper triangular elements
                upper_triangular_elements = upper_triangular_elements[2:end] # deleting that element from the set
            end
        end
    end


    if PSD == "yes"
        for i = 1:dimention
            #Q[i,i] = Q[i,i] > sum(Q[i,i+1:end]) ? Q[i,i] : sum(Q[i,i+1:end]) + 1
            Q[i,i] = Q[i,i] > (sum(Q[i, 1:i-1]) + sum(Q[i, i+1:end])) ? Q[i,i] : sum(Q[i,i+1:end])*2 + 100
        end
    end

    return Q
end


# generating matrices Qsi for the left hand side of the contraint for each of the scenario
constraint_Qs = Array{Any}(undef, 1, number_of_scenarios)

#[ constraint_Qs[i] = [ Matrix(Symmetric(sprand(number_of_continuos_decision_variables, number_of_continuos_decision_variables, Qdensity)))
    #for j = 1 : number_of_constrains] for i = 1:number_of_scenarios ]
[ constraint_Qs[i] = [ matrix_generation(Qdensity, number_of_integer_decision_variables, Max_value_for_matrix_elements, "yes")
        for j = 1 : number_of_constrains] for i = 1:number_of_scenarios ]

# generating affine functions' coefficients for the left hand side of the constraint for each of the scenario
constraint_Fs = Array{Any}(undef, 1, number_of_scenarios)

[ constraint_Fs[i] = [ rand(2, number_of_continuos_decision_variables)
    for j = 1:number_of_constrains ] for i = 1:number_of_scenarios ]
# first row - x_coeficients (continuous variables)
# second row - y_coeficients (iteger variables)

#---------------generating obejctive fucntions for the scenarios----------------

# generating matrices Qsi for the objecyive for each of the scenario
objective_Qs = Array{Any}(undef, 1, number_of_scenarios)

#[ objective_Qs[i] =  Matrix(Symmetric(sprand(number_of_continuos_decision_variables, number_of_continuos_decision_variables, Qdensity)))
     #for i = 1:number_of_scenarios ]
[ objective_Qs[i] = matrix_generation(Qdensity, number_of_integer_decision_variables, Max_value_for_matrix_elements, "yes")
     for i = 1:number_of_scenarios ]

# generating linear functions' coefficients for the objective for each of the scenario
objective_Fs = Array{Any}(undef, 1, number_of_scenarios)

[ objective_Fs[i] = [rand(2, number_of_continuos_decision_variables) ; [rand(1,1) zeros(1, number_of_continuos_decision_variables-1)] ] for i = 1:number_of_scenarios  ]
# first row - x_coeficients (continuous variables)
# second row - y_coeficients (iteger variables)
# third row  - constant

#--------------generating JuMP subproblems--------------------------------------

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

vector_of_lambda_lagrangian = Array{Any}(undef, 1, number_of_scenarios-1)
[ vector_of_lambda_lagrangian[i] = 0.5 .+ 1.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:number_of_scenarios-1 ]

# creatinf the array of subproblems
subproblem  = Array{Any}(undef, 1, number_of_scenarios)

# formulating the subproblems
for s = 1:number_of_scenarios

    global subproblem[s] = Model(with_optimizer(Gurobi.Optimizer, BarQCPConvTol = 0.4, PSDTol = 100))
    @variable(subproblem[s], x[1 : number_of_continuos_decision_variables])
    @variable(subproblem[s], y[1 : number_of_integer_decision_variables], Int )
    @objective(subproblem[s], Max,  x' * objective_Qs[s] * x + sum( ( x .* objective_Fs[s][1, :] ) .+ ( y .* objective_Fs[s][2, :] ) )   + objective_Fs[s][3, 1] )
         # +  sum( f_lambda_lagrangian(vector_of_lambda_lagrangian, s ) .* x ) )

    for i = 1:number_of_constrains
        @constraint(subproblem[s], x' * constraint_Qs[s][i] * x + sum( ( x .* constraint_Fs[s][i][1, :] ) .+ ( y .* constraint_Fs[s][i][2, :] ) ) <= 0 )
    end

end

optimize!(subproblem[1])
