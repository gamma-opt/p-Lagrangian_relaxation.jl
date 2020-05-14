# *p*-Lagrnaigan relaxation


This package contains an implementation *p*-Lagrnagian relaxation method for the noncovnex MIQCQP problems. It includes the primal model formulation, it's rnmdt and p-LR relaxations and dynamic precision-based algorithm realisation in the context of RNMDT and p-Lagranion relaxation methods. The file [experiments]() contains the application of the p-Lagrangian relaxation, RNMDT dynamc-based algorithm and Gurobi solver to solve a set of MIQCOP problems discussed in the section [Experiments](#experiments)

This package is authored by *Nikita Belyak* and *Fabricio Oliveira* and is a part of the study presented in the scientific paper *p*-Lagrnagian relaxation

**Contents**

<!-- TOC -->

- [Problems formulations](#problems_formulation)
- [Solution methods](#solution_methods)
- [Experiments](#experiments)

<!-- /TOC -->

## problems formulation
The noncovex MIQCQP problem can be created by calling the function

```julia

function original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, is_fixed_int, time_limit, seed)

```
where

* number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

* Qdensity is a density of qudratic matrices,

* is_fixed_int a string parameter and should be set to "no",

* time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

* seed is a necesarry parameter for a randomisation and can be set to any integer number.

The RNMDT relaxation of the primal MIQCQP problem can be created by calling the function

```julia
function dynamic_precision_RNMDT_problem_generation(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, time_limit, seed)
```
where

* precision_p is a vector containg the values of the precision factor of each of the continuous variables

* number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

* Qdensity is a density of qudratic matrices,

* is_fixed_int a string parameter and should be set to "no",

* time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

* seed is a necesarry parameter for a randomisation and can be set to any integer number.

the vector of containg subprobles resulting from applying p-Lagrangian relxation can be created by calling the fucntion

```julia
function dynamic_precision_based_LD_RNDMT_problem_generation(precision_p, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed)
```
where

* precision_p is a vector containg the values of the precision factor of each of the continuous variables

* number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

* Qdensity is a density of qudratic matrices,

* is_fixed_int a string parameter and should be set to "no",

* time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

* seed is a necesarry parameter for a randomisation and can be set to any integer number.

the vector of containg subprobles resulting from applying p-Lagrangian relxation can be created by calling the fucntion

## solution_methods

To solve the primal MIQCQP problem one can directly call the solver

To solve the RNMDT relxation of the primal problem one can call the function

```julia
function dynamic_precision_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, seed )
```
where

* N1 is the number of variables for which the precision will increased every iteration that is not a multiple of N2

* N2 is the number of iterations at wich the precision factors of the variables that have the hightest values will be decreased to ensure the general convergence

* number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

* Qdensity is a density of qudratic matrices,

* is_fixed_int a string parameter and should be set to "no",

* time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

* seed is a necesarry parameter for a randomisation and can be set to any integer number.

To solve p-Lagrangian relaxation of the primal MIQCQP problem usin dynamic precision-based algorithm one can call the function

```julia

dynamic_precision_LD_RNMDT_algorithm(N1, N2, tolerance, time_limit, max_number_of_iterations, number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, parallelised, seed)

```

where
* N1 is the number of variables for which the precision will increased every iteration that is not a multiple of N2

* N2 is the number of iterations at wich the precision factors of the variables that have the hightest values will be decreased to ensure the general convergence

* number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

* Qdensity is a density of qudratic matrices,

* is_fixed_int a string parameter and should be set to "no",

* time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

* parallelised is a string papramter that is equal to  if   "parallelised" if the computations should be performed in parpallel scnarios based and to "non_parallelised" otherwise 

* seed is a necesarry parameter for a randomisation and can be set to any integer number.
