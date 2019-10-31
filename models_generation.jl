using Random
Random.seed!(0)

include("parameters_generation.jl")

#--------------generating original JuMP problem---------------------------------

original_problem = Model(with_optimizer(Gurobi.Optimizer))

# continuous decision variables
@variable(original_problem, x[ 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ])

# integer decision variables
@variable(original_problem, y[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int )

# objective
@objective(original_problem, Max,  sum( x[:, s]' * objective_Qs[s] * x[:, s] + sum( ( x[:, s] .* objective_fs[s][1, :] ) .+ ( y[:, s] .* objective_fs[s][2, :] ) ) + objective_fs[s][3, 1] for s in 1:number_of_scenarios) )

# constraints 4-2
@constraint(original_problem, [s = 1 : number_of_scenarios, i = 1 : number_of_constrains ],
    x[:, s]' * constraint_Qs[s][i] * x[:, s] + sum( ( x[:, s] .* constraint_fs[s][i][1, :] ) .+ ( y[:, s] .* constraint_fs[s][i][2, :] ) ) - auxilary_constant_for_affine_part <= 0 )

# constraint 4-3
@constraint(original_problem, [ s = 1 : number_of_scenarios ],
    x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

# constraint 4-4
@constraint(original_problem, [s = 1 : number_of_scenarios],
    y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

# constraint 4-5
# non-anticipativity condition from original problem 4-5
@constraint( original_problem, sum( constraint_A1[s] .* repeat(x[:, s], number_of_scenarios - 1) for s = 1 : number_of_scenarios ) .== constraint_b1 )
@constraint( original_problem, sum( constraint_B1[s] .* repeat(y[:, s], number_of_scenarios - 1) for s = 1 : number_of_scenarios ) .== constraint_b1 )

#--------------generating mixed-integer based relaxations-----------------------

# all the commments in this section refer to the paper
# Enhancing the normalized multiparametric disaggregation
# technique for mixed-integer quadratic programming

p = -2

RNMDT_problem = Model(with_optimizer(Gurobi.Optimizer))

# continuous decision variables
@variable(RNMDT_problem, x[ 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ])

# integer decision variables
@variable(RNMDT_problem, y[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int)

@variables RNMDT_problem begin
    #RNMDT variables
    x_heat_RNMDT[ 1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables, p : -1, 1 : number_of_scenarios ] # x_heat
    delta_x_RNMDT[ 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ] # delta_x_heat
    z_RNMDT[ 1 : number_of_continuos_decision_variables, p : -1, 1 : number_of_scenarios ], Bin
    w_RNMDT[ 1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ]
    delta_w_RNMDT[1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ]
end

@objective( RNMDT_problem, Max,
    sum(
        sum( objective_Qs[s][i, i] * w_RNMDT[i, i, s]
            for i = 1 : number_of_continuos_decision_variables )
        + 2 * sum(objective_Qs[s][i, j] * w_RNMDT[i, j, s]
            for i = 1 : number_of_continuos_decision_variables,
                j = i+1 : number_of_continuos_decision_variables)
        + sum( ( x[:, s] .* objective_fs[s][1, :] )
            .+ ( y[:, s] .* objective_fs[s][2, :] ) )
            + objective_fs[s][3, 1]
    for s = 1 : number_of_scenarios) ) # 26



 @constraint( RNMDT_problem, [ s = 1 : number_of_scenarios, r = 1 : number_of_constrains ],
    sum( constraint_Qs[s][r][i, i] * w_RNMDT[i, i, s]
        for i = 1 : number_of_continuos_decision_variables )
        + 2 * sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j, s]
     for i = 1 : number_of_continuos_decision_variables,
         j = i+1 : number_of_continuos_decision_variables)
 + sum( ( x[:, s] .* constraint_fs[s][r][1, :] )
     .+ ( y[:, s] .* constraint_fs[s][r][2, :] ) )
     - auxilary_constant_for_affine_part <= 0
    ) # 27

@constraint( RNMDT_problem, [ j  = 1 : number_of_continuos_decision_variables, s  = 1 : number_of_scenarios],
    x[j, s] == (x_boundaries[j, 2] - x_boundaries[j, 1])  *
        ( sum( 2.0^l * z_RNMDT[j, l, s] for l = p : -1) + delta_x_RNMDT[j, s] ) ) # 28

@constraint( RNMDT_problem, [ i  = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios ],
    w_RNMDT[i, j, s] == (x_boundaries[j, 2] - x_boundaries[j, 1]) *
        ( sum( 2.0^l * x_heat_RNMDT[i, j, l, s] for l = p : -1) + delta_w_RNMDT[i, j, s] ) ) # 29

@constraint( RNMDT_problem, [ j  = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios, ],
    0 <= delta_x_RNMDT[j, s] <= 2.0^p ) # 30

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios],
    2.0^p * ( x[i, s] - x_boundaries[i, 2] ) + x_boundaries[i, 2] * delta_x_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) # 31

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios],
    delta_w_RNMDT[i, j, s] <= 2.0^p * ( x[i, s] - x_boundaries[i, 1] ) + x_boundaries[i, 1] * delta_x_RNMDT[j, s]) # 31

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios ],
    x_boundaries[i, 1] * delta_x_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) #32

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, s = 1 : number_of_scenarios],
    delta_w_RNMDT[i, j, s] <= x_boundaries[i, 2] * delta_x_RNMDT[j, s]) # 32

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1, s = 1 : number_of_scenarios],
    x_boundaries[i, 1] * z_RNMDT[j, l, s]  <= x_heat_RNMDT[i, j, l, s]) # 33

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1, s = 1 : number_of_scenarios],
    x_heat_RNMDT[i, j, l, s] <= x_boundaries[i, 2] * z_RNMDT[j, l, s]) # 33

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1, s = 1:number_of_scenarios],
    x_boundaries[i, 1] * (1 - z_RNMDT[j, l, s]) <= x[i,s] - x_heat_RNMDT[i, j, l, s] ) # 34

@constraint( RNMDT_problem, [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1, s = 1:number_of_scenarios],
    x[i,s] - x_heat_RNMDT[i, j, l, s] <= x_boundaries[i, 2] * (1 - z_RNMDT[j, l, s]) ) # 34

# box constraints from original problem 4-3
@constraint( RNMDT_problem, [ s = 1 : number_of_scenarios ],
    x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

# box constraints from original problem 4-4
@constraint( RNMDT_problem, [s = 1 : number_of_scenarios],
    y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

# non-anticipativity condition from original problem 4-5
@constraint( RNMDT_problem, sum( constraint_A1[s] .* repeat(x[:, s], number_of_scenarios - 1) for s = 1 : number_of_scenarios ) .== constraint_b1 )
@constraint( RNMDT_problem, sum( constraint_B1[s] .* repeat(y[:, s], number_of_scenarios - 1) for s = 1 : number_of_scenarios ) .== constraint_b1 )

#--------------generating JuMP subproblems--------------------------------------

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

vector_of_lambda_lagrangian = Array{Any}(undef, 1, number_of_scenarios-1)
[ vector_of_lambda_lagrangian[i] = 0.5 .+ 1.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:number_of_scenarios-1 ]

vector_of_mu_lagrangian = Array{Any}(undef, 1, number_of_scenarios-1)
[ vector_of_mu_lagrangian[i] = 0.5 .+ 1.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:number_of_scenarios-1 ]


# creating the array of subproblems
subproblem  = Array{Any}(undef, 1, number_of_scenarios)

# formulating the subproblems
for s = 1:number_of_scenarios

    global subproblem[s] = Model(with_optimizer(Gurobi.Optimizer,  OutputFlag = 0))

    # continuous decision variables
    @variable( subproblem[s], x[1 : number_of_continuos_decision_variables] )

    # integer decision variables
    @variable( subproblem[s], y[1 : number_of_integer_decision_variables], Int )

    @variables subproblem[s] begin
        #RNMDT variables
        x_heat_RNMDT[ 1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables, p : -1 ] # x_heat
        delta_x_RNMDT[ 1 : number_of_continuos_decision_variables ] # delta_x_heat
        z_RNMDT[ 1 : number_of_continuos_decision_variables, p : -1 ], Bin
        w_RNMDT[ 1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables ]
        delta_w_RNMDT[1 : number_of_continuos_decision_variables, 1 : number_of_continuos_decision_variables ]
    end

    # objective with relaxed terms
    @objective( subproblem[s], Max,
            sum( objective_Qs[s][i, i] * w_RNMDT[i, i]
                for i = 1 : number_of_continuos_decision_variables )
            + 2 * sum(objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuos_decision_variables,
                    j = i+1 : number_of_continuos_decision_variables)
            + sum( ( x .* objective_fs[s][1, :] )
                .+ ( y .* objective_fs[s][2, :] ) )
                + objective_fs[s][3, 1]
            +  sum( f_lambda_lagrangian(vector_of_lambda_lagrangian, s ) .* x )

            +  sum( f_lambda_lagrangian(vector_of_mu_lagrangian, s ) .* y )
                )

    @constraint( subproblem[s], [ r = 1 : number_of_constrains ],
            sum( constraint_Qs[s][r][i, i] * w_RNMDT[i, i]
                for i = 1 : number_of_continuos_decision_variables )
            + 2 * sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuos_decision_variables,
                j = i+1 : number_of_continuos_decision_variables)
            + sum( ( x .* constraint_fs[s][r][1, :] )
                .+ ( y .* constraint_fs[s][r][2, :] ) )
                - auxilary_constant_for_affine_part <= 0
            ) # 27

    @constraint( subproblem[s], [ j  = 1 : number_of_continuos_decision_variables],
        x[j] == (x_boundaries[j, 2] - x_boundaries[j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = p : -1) + delta_x_RNMDT[j] ) ) # 28

    @constraint( subproblem[s], [ i  = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables ],
        w_RNMDT[i, j] == (x_boundaries[j, 2] - x_boundaries[j, 1]) *
            ( sum( 2.0^l * x_heat_RNMDT[i, j, l] for l = p : -1) + delta_w_RNMDT[i, j] ) ) # 29

    @constraint( subproblem[s], [ j  = 1 : number_of_continuos_decision_variables ],
        0 <= delta_x_RNMDT[j] <= 2.0^p ) # 30

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables],
        2.0^p * ( x[i] - x_boundaries[i, 2] ) + x_boundaries[i, 2] * delta_x_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables],
        delta_w_RNMDT[i, j] <= 2.0^p * ( x[i] - x_boundaries[i, 1] ) + x_boundaries[i, 1] * delta_x_RNMDT[j]) # 31

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables],
        x_boundaries[i, 1] * delta_x_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables],
        delta_w_RNMDT[i, j] <= x_boundaries[i, 2] * delta_x_RNMDT[j]) # 32

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1],
        x_boundaries[i, 1] * z_RNMDT[j, l]  <= x_heat_RNMDT[i, j, l]) # 33

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1],
        x_heat_RNMDT[i, j, l] <= x_boundaries[i, 2] * z_RNMDT[j, l]) # 33

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1],
        x_boundaries[i, 1] * (1 - z_RNMDT[j, l]) <= x[i] - x_heat_RNMDT[i, j, l] ) # 34

    @constraint( subproblem[s], [ i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables, l = p : -1],
        x[i] - x_heat_RNMDT[i, j, l] <= x_boundaries[i, 2] * (1 - z_RNMDT[j, l]) ) # 34

    # box constraints from original problem 4-3
    @constraint( subproblem[s],
        x_boundaries[:, 1] .<= x .<= x_boundaries[:, 2])

    # box constraints from original problem 4-4
    @constraint( subproblem[s],
        y_boundaries[:, 1] .<= y .<= y_boundaries[:, 2])

end



# auxiliary function for printing the results
function auxiliary_print(model)
    optimize!(model)
    x = value.(model[:x])
    y = value.(model[:y])

    objective_value = sum( x[:, s]' * objective_Qs[s] * x[:, s] + sum( ( x[:, s] .* objective_fs[s][1, :] ) .+ ( y[:, s] .* objective_fs[s][2, :] ) )  + objective_fs[s][3, 1] for s in 1:number_of_scenarios)
    print("\n")
    print("objective = $objective_value \n")
    print("\n")

    for s = 1: number_of_scenarios
        for i = 1:number_of_constrains
            result = ( x[:, s]' * constraint_Qs[s][i] * x[:, s]  + sum( ( x[:, s] .* constraint_fs[s][i][1, :] ) .+ ( y[:, s] .* constraint_fs[s][i][2, :] ) )  - auxilary_constant_for_affine_part )
            print("scenario $s, constraint $i, result = $result\n \n")
        end
    end
end

#print( sum( constraint_A1[s] * x[:, s] + constraint_B1[s] * y[:, s] for s = 1 : number_of_scenarios) .== constraint_b1 )

#for s = 1: number_of_scenarios
    #for i = 1:number_of_constrains
        #print("scenario $s, constraint $i, PSD: $(isposdef(constraint_Qs[s][i]))  ")
        #print("matrix = $(constraint_Qs[s][i]) \n \n")
    #end
#end

#for s = 1: number_of_scenarios
        #print("scenario $s,  PSD: $(isposdef(objective_Qs[s])) \n")
        #print("matrix = $(objective_Qs[s]) \n \n")

    #end

#auxiliary_print(RNMDT_problem)

#value.(RNMDT_problem[:x])
#value.(RNMDT_problem[:y])
#value.(RNMDT_problem[:w_RNMDT])
#auxiliary_print(original_problem)
