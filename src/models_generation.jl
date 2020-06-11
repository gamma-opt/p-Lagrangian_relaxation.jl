using JuMP, Gurobi, Ipopt


#--------------generating original MIQCQP JuMP problem---------------------------------

function original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, is_fixed_int, time_limit, seed)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    # depending on whether the variables are fixed or not we use Ipopt and Gurobi as a solver respectively
    original_problem = Model( (is_fixed_int == "yes") ? ( with_optimizer(Ipopt.Optimizer, print_level=0) ) : ( with_optimizer(Gurobi.Optimizer, NonConvex = 2, MIPGap =  0, Method = 4, OutputFlag=0, TimeLimit = time_limit, Threads = 1) ) )

    # integer decision variables
    if (is_fixed_int == "yes")
        @variable(original_problem, x[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ] )
    else
        @variable(original_problem, x[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int )
    end

    # continuous decision variables
    @variable(original_problem, y[ 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ])

    # quadratic objective
    @objective(original_problem, Max,
        sum( round( (1/number_of_scenarios), digits = 3) *
            (
                sum( y[i, s] * objective_Qs[s][i, j] * y[j, s] for i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables)
                + sum( x[i, s] * objective_c[i] for i = 1 : number_of_integer_decision_variables)
                + sum( y[i, s] * objective_fs[s][i] for i = 1 : number_of_continuous_decision_variables)
            )
        for s in 1:number_of_scenarios)
        )

    # quadratic constraints
    @constraint(original_problem, [s = 1 : number_of_scenarios, i = 1 : number_of_constraints ],
        sum( y[j, s] * constraint_Qs[s][i][j,k] * y[k, s] for j = 1 : number_of_continuous_decision_variables, k = 1: number_of_continuous_decision_variables)
        + sum( x[j, s] * constraint_fs[s][i][1, j] for j = 1 : number_of_integer_decision_variables)
        + sum( y[j, s] * constraint_fs[s][i][2, j] for j = 1:  number_of_continuous_decision_variables)
        + constraint_fs[s][i][3, 1] <= 0 )

    # box constraints for integer variables
    @constraint(original_problem, [ s = 1 : number_of_scenarios ],
        x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

    # box constraints for continuous variables
    @constraint(original_problem, [s = 1 : number_of_scenarios],
        y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

    # non-anticipativity conditions
    @constraint( original_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_integer_decision_variables],
        x[i, s] - x[i, 1] == 0 )

    return original_problem, objective_Qs, objective_fs, objective_c, constraint_Qs, constraint_fs

end

#--------------generating mixed-integer based relaxations-----------------------

function dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, time_limit, seed)

    # all the commments for the auxiliary constraints in this section refer to the paper
    # Enhancing the normalized multiparametric disaggregation
    # technique for mixed-integer quadratic programming

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    RNMDT_problem = Model(with_optimizer(Gurobi.Optimizer, OutputFlag=0, MIPGap =  0, TimeLimit = time_limit, Threads = 1, Method = 4))#, MIPGap =  0, IntFeasTol = 1e-9, FeasibilityTol = 1e-9 ))

    # integer decision variables
    @variable(RNMDT_problem, x[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int)

    # continuous decision variables
    @variable(RNMDT_problem, y[ 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ])

    #auxiliary RNMDT variables
    @variables RNMDT_problem begin
        y_heat_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1, 1 : number_of_scenarios ] # x_heat
        delta_y_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ] # delta_x_heat
        z_RNMDT[ 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1, 1 : number_of_scenarios ], Bin
        w_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ]
        delta_w_RNMDT[1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables, 1 : number_of_scenarios ]
    end

    #quadratic objective
    @objective( RNMDT_problem, Max,
        sum( round( (1/number_of_scenarios), digits = 3) *
            (
            sum(objective_Qs[s][i, j] * w_RNMDT[i, j, s]
                    for i = 1 : number_of_continuous_decision_variables,
                        j = 1 : number_of_continuous_decision_variables)
                + sum( x[i, s] * objective_c[i] for i = 1:number_of_integer_decision_variables)
                + sum( y[j, s] * objective_fs[s][j] for j = 1:number_of_continuous_decision_variables)
            )
        for s = 1 : number_of_scenarios)
    ) # 26

    #quadratic constraints
    @constraint( RNMDT_problem, [ s = 1 : number_of_scenarios, r = 1 : number_of_constraints ],
        sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j, s]
            for i = 1 : number_of_continuous_decision_variables,
                j = 1 : number_of_continuous_decision_variables)
        + sum( x[i, s] * constraint_fs[s][r][1, i] for i = 1:number_of_integer_decision_variables)
        + sum( y[j, s] * constraint_fs[s][r][2, j] for j = 1:number_of_continuous_decision_variables )
        + constraint_fs[s][r][3, 1] <= 0
    ) # 27

## auxiliary RNMDT constraints

    @constraint( RNMDT_problem, [ j  = 1 : number_of_continuous_decision_variables, s  = 1 : number_of_scenarios],
        y[j, s] == (y_boundaries[j, 2] - y_boundaries[j, 1])  *
        ( sum( 2.0^l * z_RNMDT[j, l, s] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j, s] ) ) # 28

    @constraint( RNMDT_problem, [ i  = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios ],
        w_RNMDT[i, j, s] == (y_boundaries[j, 2] - y_boundaries[j, 1]) *
        ( sum( 2.0^l * y_heat_RNMDT[i, j, l, s] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j, s] ) ) # 29

    @constraint( RNMDT_problem, [ j  = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios ],
        0 <= delta_y_RNMDT[j, s] <= 2.0^precision_p[j, s] ) # 30

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios],
        2.0^precision_p[i, s] * ( y[i, s] - y_boundaries[i, 2] ) + y_boundaries[i, 2] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) # 31

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios],
        delta_w_RNMDT[i, j, s] <= 2.0^precision_p[i, s] * ( y[i, s] - y_boundaries[i, 1] ) + y_boundaries[i, 1] * delta_y_RNMDT[j, s]) # 31

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios ],
        y_boundaries[i, 1] * delta_y_RNMDT[j, s] <= delta_w_RNMDT[i, j, s]) #32

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios],
        delta_w_RNMDT[i, j, s] <= y_boundaries[i, 2] * delta_y_RNMDT[j, s]) # 32

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios, l = precision_p[j, s] : -1 ],
        y_boundaries[i, 1] * z_RNMDT[j, l, s]  <= y_heat_RNMDT[i, j, l, s]) # 33

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios, l = precision_p[j, s] : -1 ],
        y_heat_RNMDT[i, j, l, s] <= y_boundaries[i, 2] * z_RNMDT[j, l, s]) # 33

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios, l = precision_p[j, s] : -1 ],
        y_boundaries[i, 1] * (1 - z_RNMDT[j, l, s]) <= y[i,s] - y_heat_RNMDT[i, j, l, s] ) # 34

    @constraint( RNMDT_problem, [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, s = 1 : number_of_scenarios, l = precision_p[j, s] : -1 ],
        y[i,s] - y_heat_RNMDT[i, j, l, s] <= y_boundaries[i, 2] * (1 - z_RNMDT[j, l, s]) ) # 34

##
    # box constraints for integer variables
    @constraint( RNMDT_problem, [ s = 1 : number_of_scenarios ],
        x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])

    # box constraints for continuous variables
    @constraint( RNMDT_problem, [s = 1 : number_of_scenarios],
        y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])

    # non-anticipativity conditions
    @constraint( RNMDT_problem, [s in 2 : number_of_scenarios, i = 1 : number_of_integer_decision_variables],
        x[i, s] - x[i, 1] == 0 )

    return RNMDT_problem

end

#--------------generating subproblems of p_Lagrangian relaxation----------------

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

function dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    # generating the parameters
    constraint_Qs, constraint_fs, objective_Qs, objective_fs, objective_c, x_boundaries, y_boundaries = parameters_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)

    # randomized lagrangian relaxation variables for the x and y non anticipativity conditions written in the column
    Random.seed!(seed)
    vector_of_lambda_lagrangian = Array{Any}(undef, number_of_scenarios - 1)
    [ vector_of_lambda_lagrangian[i] = 0 .+ 0.0 .* rand(1, number_of_integer_decision_variables)
            for i = 1 : number_of_scenarios - 1 ]

    # creating the array of subproblems
    subproblem = Array{Any}(undef, 1, number_of_scenarios)

    # formulating the subproblems based on the scenarios
    for s = 1 : number_of_scenarios

        subproblem[s] = Model(with_optimizer(Gurobi.Optimizer, Method = 4, OutputFlag=0, MIPGap =  0, Threads = 1))

        # integer decision variables
        @variable( subproblem[s], x[1 : number_of_integer_decision_variables], Int )

        # continuous decision variables
        @variable( subproblem[s], y[1 : number_of_continuous_decision_variables] )

        # RNMDT variables
        @variables subproblem[s] begin
            y_heat_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1 ] # x_heat
            delta_y_RNMDT[ 1 : number_of_continuous_decision_variables ] # delta_y_heat
            z_RNMDT[ 1 : number_of_continuous_decision_variables, minimum(precision_p) : -1 ], Bin
            w_RNMDT[ 1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables ]
            delta_w_RNMDT[1 : number_of_continuous_decision_variables, 1 : number_of_continuous_decision_variables ]

        end

        # quadratic objective
        @objective( subproblem[s], Max,
            round( (1/number_of_scenarios), digits = 3) *

            ( sum(objective_Qs[s][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuous_decision_variables,
                    j = 1 : number_of_continuous_decision_variables)
            + sum( x[i] * objective_c[i]  for i = 1:number_of_integer_decision_variables)
            + sum( y[j] * objective_fs[s][j]  for j = 1:number_of_continuous_decision_variables)
            )

            +  sum( f_lambda_lagrangian( vector_of_lambda_lagrangian[:], s ) .* x )

        )

        #quadratic constraint
        @constraint( subproblem[s], [ r = 1 : number_of_constraints ],
            sum(constraint_Qs[s][r][i, j] * w_RNMDT[i, j]
                for i = 1 : number_of_continuous_decision_variables,
                j = 1 : number_of_continuous_decision_variables)
            + sum( x[i] * constraint_fs[s][r][1, i] for i = 1:number_of_integer_decision_variables)
            + sum( y[j] * constraint_fs[s][r][2, j] for j = 1:number_of_continuous_decision_variables )
            + constraint_fs[s][r][3, 1] <= 0
        ) # 27

## auxiliary RNMDT constraints

        @constraint( subproblem[s], [ j  = 1 : number_of_continuous_decision_variables],
            y[j] == (y_boundaries[j, 2] - y_boundaries[j, 1])  *
            ( sum( 2.0^l * z_RNMDT[j, l] for l = precision_p[j, s] : -1) + delta_y_RNMDT[j] ) ) # 28

        @constraint( subproblem[s], [ i  = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables ],
            w_RNMDT[i, j] == (y_boundaries[j, 2] - y_boundaries[j, 1]) *
            ( sum( 2.0^l * y_heat_RNMDT[i, j, l] for l = precision_p[j, s] : -1) + delta_w_RNMDT[i, j] ) ) # 29

        @constraint( subproblem[s], [ j  = 1 : number_of_continuous_decision_variables ],
            0 <= delta_y_RNMDT[j] <= 2.0^precision_p[j, s] ) # 30

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            2.0^precision_p[i, s] * ( y[i] - y_boundaries[i, 2] ) + y_boundaries[i, 2] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) # 31

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            delta_w_RNMDT[i, j] <= 2.0^precision_p[i, s] * ( y[i] - y_boundaries[i, 1] ) + y_boundaries[i, 1] * delta_y_RNMDT[j]) # 31

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            y_boundaries[i, 1] * delta_y_RNMDT[j] <= delta_w_RNMDT[i, j]) #32

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables],
            delta_w_RNMDT[i, j] <= y_boundaries[i, 2] * delta_y_RNMDT[j]) # 32

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_boundaries[i, 1] * z_RNMDT[j, l]  <= y_heat_RNMDT[i, j, l]) # 33

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_heat_RNMDT[i, j, l] <= y_boundaries[i, 2] * z_RNMDT[j, l]) # 33

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y_boundaries[i, 1] * (1 - z_RNMDT[j, l]) <= y[i] - y_heat_RNMDT[i, j, l] ) # 34

        @constraint( subproblem[s], [ i = 1 : number_of_continuous_decision_variables, j = 1 : number_of_continuous_decision_variables, l = precision_p[j, s] : -1],
            y[i] - y_heat_RNMDT[i, j, l] <= y_boundaries[i, 2] * (1 - z_RNMDT[j, l]) ) # 34
##
        # box constraints for integer variables
        @constraint( subproblem[s],
            x_boundaries[:, 1] .<= x .<= x_boundaries[:, 2])

        # box constraints for continuous variables
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
