using JuMP, Gurobi, Ipopt, Cbc

include("parameters_generation.jl")

#--------------generating original JuMP problem---------------------------------

function original_problem_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    original_problem = Model(with_optimizer(Ipopt.Optimizer))

    # continuous decision variables
    @variable(original_problem, x[ 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ])

    # integer decision variables
    @variable(original_problem, y[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int )

    # objective
    @objective(original_problem, Max,
        sum( (1/number_of_scenarios) *
            (
                sum( x[i, s] * objective_Qs[s][i, j] * x[j, s] for i = 1 : number_of_continuos_decision_variables, j = 1 : number_of_continuos_decision_variables)
                + sum( x[i, s] * objective_fs[s][1, i] for i = 1 : number_of_continuos_decision_variables)
                + sum( y[i, s] * objective_fs[s][2, i] for i = 1 : number_of_integer_decision_variables)
            )
        for s in 1:number_of_scenarios)
        )

    # constraints 4-2
    @constraint(original_problem, [s = 1 : number_of_scenarios, i = 1 : number_of_constrains ],
        sum( x[j, s]' * constraint_Qs[s][i][j,k] * x[k, s] for j = 1 : number_of_continuos_decision_variables, k = 1: number_of_continuos_decision_variables)
        + sum( x[j, s] * constraint_fs[s][i][1, j] for j = 1 : number_of_continuos_decision_variables)
        + sum( y[j, s] * constraint_fs[s][i][2, j] for j = 1:  number_of_integer_decision_variables)
        + constraint_fs[s][i][3, 1] <= 0 )

    # constraint 4-3
    @constraint(original_problem, [ s = 1 : number_of_scenarios ],
        x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

    # constraint 4-4
    @constraint(original_problem, [s = 1 : number_of_scenarios],
        y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

    # constraint 4-5
    # non-anticipativity condition from original problem 4-5
    @constraint( original_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_continuos_decision_variables],
        x[i, s] - x[i, 1] == 0 )
    @constraint( original_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_integer_decision_variables],
        y[i, s] - y[i, 1] == 0 )

    return original_problem

end

#--------------generating mixed-integer based relaxations-----------------------

function RNMDT_problem_generation(p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # all the commments in this section refer to the paper
    # Enhancing the normalized multiparametric disaggregation
    # technique for mixed-integer quadratic programming

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

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
        sum( (1/number_of_scenarios) *
            (
                sum( objective_Qs[s][i, i] * w_RNMDT[i, i, s]
                    for i = 1 : number_of_continuos_decision_variables )
                + 2 * sum(objective_Qs[s][i, j] * w_RNMDT[i, j, s]
                    for i = 1 : number_of_continuos_decision_variables,
                        j = i+1 : number_of_continuos_decision_variables)
                + sum( x[i, s] * objective_fs[s][1, i]  for i = 1:number_of_continuos_decision_variables)
                + sum( y[j, s] * objective_fs[s][2, j]  for j = 1:number_of_integer_decision_variables)
            )
        for s = 1 : number_of_scenarios)
    ) # 26



    @constraint( RNMDT_problem, [ s = 1 : number_of_scenarios, r = 1 : number_of_constrains ],
        sum( constraint_Qs[s][r][i, i] * w_RNMDT[i, i, s]
            for i = 1 : number_of_continuos_decision_variables )
        + 2 * sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j, s]
            for i = 1 : number_of_continuos_decision_variables,
            j = i+1 : number_of_continuos_decision_variables)
        + sum( x[i, s] * constraint_fs[s][r][1, i] for i = 1:number_of_continuos_decision_variables)
        + sum( y[j, s] * constraint_fs[s][r][2, j] for j = 1:number_of_integer_decision_variables )
        + constraint_fs[s][r][3, 1] <= 0
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
    @constraint( RNMDT_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_continuos_decision_variables],
        x[i, s] - x[i, 1] == 0 )
    @constraint( RNMDT_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_integer_decision_variables],
        y[i, s] - y[i, 1] == 0 )

    return RNMDT_problem

end

#--------------generating JuMP subproblems--------------------------------------

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

function LD_RNDMT_problem_generation(p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)

    # randomized lagrangian relaxation variables for the x and y non anticipativity conditions written in the column
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_scenarios - 1)
    [ vector_of_lambda_lagrangian[i] = 0 .+ 20.0 .* rand(1,
        number_of_continuos_decision_variables + number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the x variables
    x_indices = 1 : number_of_continuos_decision_variables

    # indices of coordinates of Lagrangian multipliers (written in the column) correspondent to the y variables
    y_indices = number_of_continuos_decision_variables + 1 : number_of_continuos_decision_variables + number_of_integer_decision_variables

    # creating the array of subproblems
    subproblem = Array{Any}(undef, 1, number_of_scenarios)

    # formulating the subproblems
    for s = 1 : number_of_scenarios

        subproblem[s] = Model(with_optimizer(Gurobi.Optimizer,  OutputFlag = 0))

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
            (1/number_of_scenarios) *

            ( sum( objective_Qs[s][i, i] * w_RNMDT[i, i]
                for i = 1 : number_of_continuos_decision_variables )
            + 2 * sum(objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuos_decision_variables,
                    j = i+1 : number_of_continuos_decision_variables)
            + sum( x[i] * objective_fs[s][1, i]  for i = 1:number_of_continuos_decision_variables)
            + sum( y[j] * objective_fs[s][2, j]  for j = 1:number_of_integer_decision_variables)
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s )[x_indices] .* x )
            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s )[y_indices] .* y )
        )

        @constraint( subproblem[s], [ r = 1 : number_of_constrains ],
            sum( constraint_Qs[s][r][i, i] * w_RNMDT[i, i]
                for i = 1 : number_of_continuos_decision_variables )
            + 2 * sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuos_decision_variables,
                j = i+1 : number_of_continuos_decision_variables)
            + sum( x[i] * constraint_fs[s][r][1, i] for i = 1:number_of_continuos_decision_variables)
            + sum( y[j] * constraint_fs[s][r][2, j] for j = 1:number_of_integer_decision_variables )
            + constraint_fs[s][r][3, 1] <= 0
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

    return subproblem

end

#--------------generating MPS file ---------------------------------------------
#using MathOptFormat
#mps_model = MathOptFormat.MPS.Model()
#original_problem_backend = JuMP.backend(original_problem)
#MOI.copy_to(mps_model, original_problem_backend)
#MOI.write_to_file(mps_model, "my_model.mps")
