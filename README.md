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
 ``\LaTeX
\begin{align}
\text{RDE}: \text{max. } \quad  &  \sum_{s \in {\mathcal{S}}} P^s \left( \sum_{j \in \mathcal{VI}} I^{0}_j x_{j}^{s}  +  \sum_{i \in {\mathcal{VC}}}\sum_{j \in {\mathcal{VC}}} Q^{s, 0}_{i, j} y_{i}^{s} y_{j}^{s}  + \sum_{i \in {\mathcal{VC}}} C^{s,0}_{i} y_{i}^{s}\right) \label{decomposable_discretised_original_problem_objective} \\ \nonumber \\
\text{s.t.: } \quad & \sum_{i \in {\mathcal{VC}}}\sum_{j \in {\mathcal{VC}}} Q^{s, r}_ {i, j} y_{i}^{s} y_{j}^{s}   + \sum_{i \in {\mathcal{VC}}} C^{s,r}_{i} y_{i}^{s} +  \sum_{j \in {\mathcal{VI}}} I^{s, r}_{j} x_{j}^{s}  + K^{s,r}  \le 0,  \nonumber \\ & \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad \quad  \forall s \in \mathcal{S}, \forall r \in \mathcal{R} \label{decomposable_discretised_original_problem_quadratic_constraint} \\
 & x_{j}^{s} \in \{X^L_{j}, \dots, X^U_{j}\},    \forall s \in \mathcal{S}  , \forall j \in \mathcal{VI} \label{decomposable_discretised_original_problem_box_constraint_integer} \\
 &  x_{j}^{s^\prime} - x_{j}^{s} = 0 ,   \forall s \in \mathcal{S} \setminus \{s^\prime\},  \forall j \in \mathcal{VI}  \label{decomposable_discretised_original_problem_non-anticipativity_constraint_integer} \\
 & \text{ and (\ref{discretised_original_problem_box_constraint_continuous}).} \nonumber&
\end{align}

``

can be created by calling the function

```julia

function original_problem_generation(number_of_scenarios, number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints, Qdensity, is_fixed_int, time_limit, seed)

```
where

number_of_integer_decision_variables, number_of_continuous_decision_variables, number_of_constraints are defined per scenario,

Qdensity is a density of qudratic matrices,

is_fixed_int a string parameter and should be set to "no",

time_limit corresposnds to the time limit for gurobi solver allowed for solving the problem

seed is a necesarry parameter for a randomisation and can be set to any integer number. s
