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
for s = 1 : number_of_scenarios
    for i = 1 : number_of_constrains
        @constraint(original_problem, x[:, s]' * constraint_Qs[s][i] * x[:, s] + sum( ( x[:, s] .* constraint_fs[s][i][1, :] ) .+ ( y[:, s] .* constraint_fs[s][i][2, :] ) ) - auxilary_constant_for_affine_part <= 0 )
    end
end

# constraint 4-3
for s = 1 : number_of_scenarios
    @constraint(original_problem, x_boundaries[:, 1] .<= x[:, s] .<= x_boundaries[:, 2])
end

# constraint 4-4
for s = 1 : number_of_scenarios
    @constraint(original_problem, y_boundaries[:, 1] .<= y[:, s] .<= y_boundaries[:, 2])
end

# constraint 4-5
    @constraint(original_problem, sum( constraint_A1[s] * x[:, s] + constraint_B1[s] * y[:, s] for s = 1 : number_of_scenarios) .== constraint_b1)

#--------------generating mixed-integer based relaxations-----------------------
p = -4

RNMDT_problem = Model(with_optimizer(Gurobi.Optimizer))

# continuous decision variables
@variable(RNMDT_problem, x[ 1 : number_of_continuos_decision_variables, 1 : number_of_scenarios ])

# integer decision variables
@variable(RNMDT_problem, y[ 1 : number_of_integer_decision_variables, 1 : number_of_scenarios ], Int)

@variables RNMDT_problem begin
    #RNMDT variables
    q_aux_RNMDT[1 : number_of_integer_decision_variables, p:-1] #x_heat
    delta_q_RNMDT[1 : number_of_integer_decision_variables] >= 0 #delta_x_heat
    z_RNMDT[1 : number_of_integer_decision_variables, p:-1], Bin
    w_RNMDT[1 : number_of_integer_decision_variables]
    delta_w_RNMDT[1 : number_of_integer_decision_variables]
end
    #@objective(RNMDT_problem, Max,
    #    sum(objective_Qs[s][i, i] * w[i, i] for i = 1 : number_of_constraints)
    #    + sum(objective_Qs[s][i, i] * w[i, i] for i = 1 : )                            )
    #end

    #--------------generating JuMP subproblems--------------------------------------

# auxiliary function for Lagrangian multipliers update
f_lambda_lagrangian(lambda_lagrangian, dec_index) = (dec_index == 1 ? sum(lambda_lagrangian[1:end]) : - lambda_lagrangian[dec_index-1])

vector_of_lambda_lagrangian = Array{Any}(undef, 1, number_of_scenarios-1)
[ vector_of_lambda_lagrangian[i] = 0.5 .+ 1.0 .* rand(1, number_of_continuos_decision_variables) for i = 1:number_of_scenarios-1 ]

# creating the array of subproblems
subproblem  = Array{Any}(undef, 1, number_of_scenarios)

# formulating the subproblems
for s = 1:number_of_scenarios

    global subproblem[s] = Model(with_optimizer(Gurobi.Optimizer))

    # continuous decision variables
    @variable(subproblem[s], x[1 : number_of_continuos_decision_variables])

    # integer decision variables
    @variable(subproblem[s], y[1 : number_of_integer_decision_variables], Int )

    # objective with relaxed
    @objective(subproblem[s], Max,  x' * objective_Qs[s] * x + sum( ( x .* objective_fs[s][1, :] ) .+ ( y .* objective_fs[s][2, :] ) ) + objective_fs[s][3, 1]
        +  sum( f_lambda_lagrangian(vector_of_lambda_lagrangian, s ) .* y ) )

    # constraints 4-2
    for i = 1:number_of_constrains
        @constraint(subproblem[s], x' * constraint_Qs[s][i] * x + sum( ( x .* constraint_fs[s][i][1, :] ) .+ ( y .* constraint_fs[s][i][2, :] ) ) <= 0 )
    end

    # constraint 4-3
    @constraint(subproblem[s], x_boundaries[:, 1] .<= x .<= x_boundaries[:, 2])

    # constraint 4-4
    @constraint(subproblem[s], y_boundaries[:, 1] .<= y .<= y_boundaries[:, 2])
end

optimize!(original_problem)

value.(x)

x = zeros(3,3)
y = 0.5 .* ones(3,3)

for s = 1 : number_of_scenarios
    for i = 1 : number_of_constrains
        print( x[:, s]' * constraint_Qs[s][i] * x[:, s] + sum( ( x[:, s] .* constraint_fs[s][i][1, :] ) .+ ( y[:, s] .* constraint_fs[s][i][2, :] ) ) <= 0 )
    end
end

print(sum( x[:, s]' * objective_Qs[s] * x[:, s] + sum( ( x[:, s] .* objective_fs[s][1, :] ) .+ ( y[:, s] .* objective_fs[s][2, :] ) )  + objective_fs[s][3, 1] for s in 1:number_of_scenarios))

print( sum( constraint_A1[s] * x[:, s] + constraint_B1[s] * y[:, s] for s = 1 : number_of_scenarios) .== constraint_b1 )


for s = 1: number_of_scenarios
    for i = 1:number_of_constrains
        print("scenario $s, constraint $i, PSD: $(isposdef(constraint_Qs[s][i])) \n")
        print("quadmatrix = $(constraint_Qs[s][i]) \n \n")
        print()
    end
end

for s = 1: number_of_scenarios
        print("scenario $s,  PSD: $(isposdef(objective_Qs[s])) \n")
        print("matrix = $(constraint_Qs[s]) \n \n")
    end
end
