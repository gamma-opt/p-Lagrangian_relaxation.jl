# __p_Lagrangian_Decomposition-__

## The code in this folder contains the implementation of the following items:

1. generation the parameters for MIQCQP problems formulation
2. generation of the original MIQCQP formulation, mixed-integer relaxation of it based on RNMDT technique and formulation of Lagrangian dual subproblems after applying Lagrangian decomposition
3. the realisation of the Bundle method inspired Lagrangian decomposition

The code was implemented in Julia 1.2.0 using JuMP package of 0.20 version.

The folder containts the following files:

---

## __parameters_generation.jl__

---

 <span style="color:blue">function</span>
__quadratic_matrix_generation(density, dimention, min_range, max_range, PSD)__

_input_:
*  __density__ - the density of the matrix
*  __dimention__ - the dimention of the matrix
*  __min_range__ - minimum value for the matrix elements generation
*  __max_range__ - maximum value for the matrix elements generation
* __PSD__ - boolean variable with the values "true" or "false" depending on whether the matrix should be PSD or not

_output_:
* returns __dimention__ by __dimention__ matrix with the __density__ density which is PSD or not depending on the value of the __PSD__ variable

---

<span style="color:blue">function</span>
__parameters_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)__

_input_:

* __number_of_scenarios__ - number of the scenarios in the problem
* __number_of_continuos_decision_variables__ - number of the continuous decision variables per scenario
* __number_of_integer_decision_variables__ - number of the integer decision variables per scenario
* __number_of_constrains__ - number of the constraints per scenario
* __Qdensity__ - the density of the matrices for the quadratic constraints

_output_:

* returns the array containing the following structures:
 * __constraint_Qs__ - __number_of_scenarios__ * __number_of_constrains__ matrices for the constraints quadratic part formulation
 * __constraint_fs__ - __number_of_scenarios__ * __number_of_constrains__ matrices containing coefficients for the continuous and integer variables and constants utilizing for the formulating the affine part of the constraints
 * __objective_Qs__ - __number_of_scenarios__ matrices for the objective function quadratic part formulation
 * __objective_fs__ - __number_of_scenarios__ matrices containing coefficients for the continuous and integer variables utilizing for the formulating the linear part of the objective function
 * __x_boundaries__ - the boundaries for the continuous decision variables for all of the scenarios
 * __y_boundaries__ - the boundaries for the integer decision variables for all of the scenarios

---

## __models_generation.jl__

---

<span style="color:blue">function</span>
__original_problem_generation(number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)__

_input_

* __number_of_scenarios__ - number of the scenarios in the problem
* __number_of_continuos_decision_variables__ - number of the continuous decision variables per scenario
* __number_of_integer_decision_variables__ - number of the integer decision variables per scenario
* __number_of_constrains__ - number of the constraints per scenario
* __Qdensity__ - the density of the matrices for the quadratic constraints

_output_

* returns JuMP model representing the original MIQCQP model formulation

---

<span style="color:blue">function</span>
__RNMDT_problem_generation(p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)__

_input_
* __p__ - precision factor
* __number_of_scenarios__ - number of the scenarios in the problem
* __number_of_continuos_decision_variables__ - number of the continuous decision variables per scenario
* __number_of_integer_decision_variables__ - number of the integer decision variables per scenario
* __number_of_constrains__ - number of the constraints per scenario
* __Qdensity__ - the density of the matrices for the quadratic constraints

_output_

* returns JuMP model representing the mixed-integer based relaxation of the original MIQCQP based on applying the RNDMT technique with precision factor __p__

---

<span style="color:blue">function</span> __f_lambda_lagrangian(lambda_lagrangian, dec_index)__
auxiliary function for Lagrangian multipliers update

_input_
* __lambda_lagrangian__ - the vector that containts lagrangian multipliers values
* __dec_index__ - the index that corresponds to the consdering subproblem

_output_
* the relaxation part for the Lagrangian dual problem objective correspondent to the considering subproblem

---

<span style="color:blue">function</span>
__LD_RNDMT_problem_generation(p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity)__

_input_

* __p__ - precision factor
* __number_of_scenarios__ - number of the scenarios in the problem
* __number_of_continuos_decision_variables__ - number of the continuous decision variables per scenario
* __number_of_integer_decision_variables__ - number of the integer decision variables per scenario
* __number_of_constrains__ - number of the constraints per scenario
* __Qdensity__ - the density of the matrices for the quadratic constraints

_output_

* returns the array of the Lagrangian dual subproblems resulting in applying Lagrangian decomposition to the mixed-integer based relaxation of the original problem

---
## p_lagrangian_decomposition_Budnle.jl

---
<span style="color:blue">function</span>
__Lagrangian_decomposition_bundle(p, number_of_scenarios, number_of_continuos_decision_variables, number_of_integer_decision_variables, number_of_constrains, Qdensity, number_of_iterations)__

 The function containing the implementation of the Lagrangian decomposition with the bundle method multipliers update applied to the mixed integer based relaxation of the original problem

As the stopping criteria of the algorithms the fact that the difference between the values of the objective function during the last 5 iterations is small enough comparing to predefined tolerance is utilised __eps_stop__

_input_

* __p__ - precision factor
* __number_of_scenarios__ - number of the scenarios in the problem
* __number_of_continuos_decision_variables__ - number of the continuous decision variables per scenario
* __number_of_integer_decision_variables__ - number of the integer decision variables per scenario
* __number_of_constrains__ - number of the constraints per scenario
* __Qdensity__ - the density of the matrices for the quadratic constraints
* __number_of_iterations__ - limit in the number of iterations for the Budnle method

_output_

* the array containing the following items
 *  the array containing the values of the dual lagrangian function at all the interations
 * the performance time
 * the array containing the values of the centre of mass at all of the iterations
---
