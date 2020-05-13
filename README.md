# *p*-Lagrnaigan relaxation


This package contains an implementation *p*-Lagrnagian relaxation method for the noncovnex MIQCQP problems. It includes the primal model formulation, it's rnmdt and p-LR relaxations and dynamic precision-based algorithm realisation in the context of RNMDT and p-Lagranion relaxation methods. The file [experiments]() contains the application of the p-Lagrangian relaxation, RNMDT dynamc-based algorithm and Gurobi solver to solve a set of MIQCOP problems discussed in the section [Experiments](#experiments)

This package is authored by *Nikita Belyak* and *Fabricio Oliveira* and is a part of the study presented in the scientific paper *p*-Lagrnagian relaxation

**Contents**

<!-- TOC -->

- [problems formulations](#problems_formulation)
- [RNMDT relaxation](#rnmdt)
- [Budnle method inspired p-lagrangian relaxation](#p-lr)
- [Experiments](#experiments)

<!-- /TOC -->

## problems formulation
The noncovex MIQCQP problem of the form
```\LaTeX

can be created by calling the function

```julia
function original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, is_fixed_int, time_limit, seed)
```
where
number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

Qdensity is a density of qudratic matrices,

is_fixed_int parameter should be set to "no",

time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

seed is a necesarry parameter for a randomisation and can be set to any integer number. s
