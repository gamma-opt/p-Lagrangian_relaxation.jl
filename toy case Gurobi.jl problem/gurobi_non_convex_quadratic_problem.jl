using JuMP, Gurobi, Random, LinearAlgebra
Random.seed!(0)
# The formulation of the problem is attached in a separate .pdf file

#----------------Setting the parameters-----------------------------------------

number_of_integer_decision_variables = 2
number_of_continuous_decision_variables = 5
number_of_scenarios = 3
number_of_constrains = 1

# the min and max limits for the matrices elements generation
Max_value_for_matrix_elements = 100
Min_value_for_matrix_elements = 0

# the max limit for the budget for each of the constaints generation
Max_value_for_affine_constraint = 10000000

# the max and min limits for the decision variables bounds generation
x_limits = [0 100] # max and min values for the continuous variables' boundaries
y_limits = [0 100] # max and min values for the integer variables' boundaries

#---------------generating the bounds for the decision variables----------------

x_boundaries = [ x_limits[1]*ones( number_of_integer_decision_variables, 1 ) rand(
    Int( round( (x_limits[2] - x_limits[1]) / 2 ) ) : x_limits[2],
    number_of_integer_decision_variables, 1 ) ]

y_boundaries = [ y_limits[1]*ones( number_of_continuous_decision_variables, 1 ) rand(
    Int( round( (y_limits[2] - y_limits[1]) / 2 ) ) : y_limits[2],
    number_of_continuous_decision_variables, 1 ) ]

#---------------generating the matrices for the constraints---------------------

# generating matrices Qsi for the left-hand side of the constraint for each of the scenarios
constraint_Qs = Array{Any}(undef, 1, number_of_scenarios)

[ constraint_Qs[i] = [ round.(Min_value_for_matrix_elements .+ (Max_value_for_matrix_elements - Min_value_for_matrix_elements) .* Matrix(Symmetric(rand(number_of_continuous_decision_variables, number_of_continuous_decision_variables))), digits = 1)
    for j = 1 : number_of_constrains] for i = 1:number_of_scenarios ]

# generating affine functions' coefficients for the left-hand side of the constraint for each of the scenarios
constraint_fs = Array{Any}(undef, 1, number_of_scenarios)

[ constraint_fs[i] = [ round.([(Min_value_for_matrix_elements + (Max_value_for_matrix_elements - Min_value_for_matrix_elements)) .* [rand(1, number_of_integer_decision_variables) zeros(1, ( number_of_continuous_decision_variables > number_of_integer_decision_variables ) ? (number_of_continuous_decision_variables - number_of_integer_decision_variables) : 0 ) ];
                        (Min_value_for_matrix_elements + (Max_value_for_matrix_elements - Min_value_for_matrix_elements)) .* [ rand(1, number_of_continuous_decision_variables) zeros(1, ( number_of_integer_decision_variables > number_of_continuous_decision_variables ) ? (number_of_integer_decision_variables - number_of_continuous_decision_variables) : 0 ) ];
                        -Max_value_for_affine_constraint .* [rand(1,1) zeros(1, ( number_of_continuous_decision_variables > number_of_integer_decision_variables ) ? (number_of_continuous_decision_variables - 1) : (number_of_integer_decision_variables - 1))] ], digits = 1)
    for j = 1:number_of_constrains] for i = 1:number_of_scenarios ]
# first row - x_coeficients (integer variables)
# second row - y_coeficients (continuous variables)
# third row  - constant

#---------------generating the matrices for the objective-----------------------

# generating matrices Qsi for the objective function  for each of the scenarios
objective_Qs = Array{Any}(undef, 1, number_of_scenarios)

[ objective_Qs[i] = round.(Min_value_for_matrix_elements .+ (Max_value_for_matrix_elements - Min_value_for_matrix_elements) .* Matrix(Symmetric(rand(number_of_continuous_decision_variables, number_of_continuous_decision_variables))), digits = 1)
    for i = 1:number_of_scenarios ]

# generating linear functions' coefficients for the objective function for each of the scenarios
objective_fs = Array{Any}(undef, 1, number_of_scenarios)

[ objective_fs[i] = round.( Min_value_for_matrix_elements .+ (Max_value_for_matrix_elements - Min_value_for_matrix_elements) .*  rand(1, number_of_continuous_decision_variables),  digits = 1) for i = 1:number_of_scenarios  ]
# first row - x_coeficients (continuous variables)
# second row - y_coeficients (iteger variables)

objective_c = round.( Min_value_for_matrix_elements .+ (Max_value_for_matrix_elements - Min_value_for_matrix_elements) .*  rand(1, number_of_integer_decision_variables),  digits = 1)

#----------------formulating the model------------------------------------------

original_problem = Model(with_optimizer(Gurobi.Optimizer))

# integer decision variables
@variable(original_problem, x[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int)

# continuous decision variables
@variable(original_problem, y[ 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ])

# objective
@objective(original_problem, Max,
    sum( round( (1/number_of_scenarios), digits = 2) *
        (
            sum( y[i, s] * objective_Qs[s][i, j] * y[j, s] for i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables)
            + sum( x[i, s] * objective_c[i] for i = 1 : number_of_integer_decision_variables)
            + sum( y[i, s] * objective_fs[s][i] for i = 1 : number_of_continuous_decision_variables)
        )
    for s in 1:number_of_scenarios)
    )

# quadratic constrints
@constraint(original_problem, [s = 1 : number_of_scenarios, i = 1 : number_of_constrains ],
    sum( y[j, s] * constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : number_of_continuous_decision_variables, k = 1: number_of_continuous_decision_variables)
    + sum( x[j, s] * constraint_fs[s][i][1, j] for j = 1 : number_of_integer_decision_variables)
    + sum( y[j, s] * constraint_fs[s][i][2, j] for j = 1:  number_of_continuous_decision_variables)
    + constraint_fs[s][i][3, 1] <= 0 )

# box constraints for integer decision variables
@constraint(original_problem, [ s = 1 : number_of_scenarios ],
    x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

# box constraints for continuous decision variables
@constraint(original_problem, [s = 1 : number_of_scenarios],
    y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

# non-anticipativity conditions
@constraint( original_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_integer_decision_variables],
    x[i, s] - x[i, 1] == 0 )

#-----------------attempt to optimise-------------------------------------------

optimize!(original_problem)
