
Variable x(i);
integer Variable y(e);
Variable y_continuous(e);
Variable w(i,j);
Binary Variable z(j,l);
Variable hat_x(i,j,l);
Positive Variable delta_x(j);
Variable v(i,j);
variable FO;

parameter
x_relaxado(i)
w_relaxado(i,j)
v_relaxado(i,j)
delta_x_relaxado(i)
pp_ws(j,sc)
maior_j_ws(j,sc)
;



Equations
R_FO
Main_L
Main_E
Main_G
Main_L_ori
Main_E_ori
Main_G_ori
*R_FO_fixed
*Main_fixed
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod
MC_1
MC_2
MC_3
MC_4
FO_LB
FO_UB
outer_approximation

Main_continuous
R_FO_Rel_continuous

MC_w_1
MC_w_2
MC_w_3
MC_w_4
;

Equations
R_FO_Or
R_FO_rel
Main_Or
R_FO_rel_lagrange
;

Variable extra_multiplicador(i);
equation eq_extra_multiplicador;
parameter par_extra_multiplicador(i,sc);
parameter coeficiente_multiplicador(i);

eq_extra_multiplicador(i)$iVc(i).. extra_multiplicador(i) =e= coeficiente_multiplicador(i)*x(i) ;

R_FO_rel(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) =E= FO;

R_FO_rel_lagrange(r)$(re_FO(r))..
         (sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
)
/card(sc)
*         + sum((i,sc)$(iVc(i) and ord(sc)<card(sc)),(multiplicadores(i,sc)*x(i))) -  sum((i,sc)$(iVc(i) and ord(sc)>1), multiplicadores(i,sc-1)*x(i))
         + sum(i$iVc(i), coeficiente_multiplicador(i)*x(i))
          =E= FO;

R_FO_rel_continuous(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_continuous(e)) =E= FO;


*R_FO_fixed(r)$(re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_fixed(e)) =E= FO;

*Main_fixed(r)$(not re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_fixed(e)) + c(r) =L= 0;

Main_Or(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) + c(r) =L= 0;

FO_LB..
         FO =G= LB;

FO_UB..
         FO =L= UB;

parameter mult(r);
mult(r) = 0;

R_FO(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(r1$(not re_FO(r1)), mult(r1)*( sum((i,j)$BL(r1,i,j),a(r1,i,j)*w(i,j)) + sum(i,b(r1,i)*x(i)) + c(r1))) =E= FO;

Main_L(r)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =L= 0;
Main_E(r)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =E= 0;
Main_G(r)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =G= 0;

****

R_FO_Or(r)$(re_FO(r))..
         (sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
         + sum(e,b_int(r,e)*y(e))
)
/card(sc)
          =E= FO;

Main_L_ori(r)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =L= 0;
Main_E_ori(r)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =E= 0;
Main_G_ori(r)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =G= 0;

****

Main_continuous(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_continuous(e)) + c(r) =L= 0;

Def_w(i,j)$QT(i,j)..
         w(i,j) =E= (x.up(j)-x.lo(j))*sum(l$(ord(l)<=-pp(j)),2**(-ord(l))*hat_x(i,j,l)) + (x.up(j) - x.lo(j))*v(i,j) + x(i)*x.lo(j);

Def_x(j)$DS(j)..
         x(j) =E= (x.up(j) - x.lo(j))*sum(l$(ord(l)<=-pp(j)),2**(-ord(l))*z(j,l)) + (x.up(j) - x.lo(j))*delta_x(j) + x.lo(j);

UB_hat_x(i,j,l)$(QT(i,j) and ord(l)<=-pp(j))..
         hat_x(i,j,l) =L= x.up(i)*z(j,l);

LB_hat_x(i,j,l)$(QT(i,j) and ord(l)<=-pp(j))..
         hat_x(i,j,l) =G= x.lo(i)*z(j,l);

UB_prod(i,j,l)$(QT(i,j) and ord(l)<=-pp(j))..
         x(i) - hat_x(i,j,l) =L= x.up(i)*(1-z(j,l));

LB_prod(i,j,l)$(QT(i,j) and ord(l)<=-pp(j))..
         x(i) - hat_x(i,j,l) =G= x.lo(i)*(1-z(j,l));

MC_1(i,j)$QT(i,j)..
         v(i,j) =G= delta_x(j)*x.lo(i);

MC_2(i,j)$QT(i,j)..
         v(i,j) =L= delta_x(j)*x.up(i);

MC_3(i,j)$QT(i,j)..
         v(i,j) =G= (2**pp(j))*(x(i)-x.up(i)) + delta_x(j)*x.up(i);

MC_4(i,j)$QT(i,j)..
         v(i,j) =L= (2**pp(j))*(x(i)-x.lo(i)) + delta_x(j)*x.lo(i);

MC_w_1(i,j)$QT(i,j)..
         w(i,j) =G= x(i)*x.up(j) + x.up(i)*x(j) - x.up(i)*x.up(j);
MC_w_2(i,j)$QT(i,j)..
         w(i,j) =G= x(i)*x.lo(j) + x.lo(i)*x(j) - x.lo(i)*x.lo(j);
MC_w_3(i,j)$QT(i,j)..
         w(i,j) =L= x(i)*x.lo(j) + x.up(i)*x(j) - x.up(i)*x.lo(j);
MC_w_4(i,j)$QT(i,j)..
         w(i,j) =L= x(i)*x.up(j) + x.lo(i)*x(j) - x.lo(i)*x.up(j);

parameter x_value(i,l);
*scalar pp /-50/;
*x_value(i,l) := LB_x(i) + (UB_x(i) - LB_x(i))*sum(ll$(ord(ll)<ord(l)), 2**(pp));
*outer_approximation(i,l)$sum(r, BL(r,i,i) and ord(l)<=-pp+2)..
*         w(i,i) =G= 2*x_value(i,l)*x(i) - power(x_value(i,l), 2);

* Bounds das variáveis
delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**pp(j);
x.lo(i)$(LB_x(i) > -power(10,10)) = LB_x(i);
x.up(i)$(UB_x(i) and UB_x(i) < power(10,10)) = UB_x(i);
y.lo(e) = LB_y(e);
y.up(e) = UB_y(e);
y_continuous.lo(e) = LB_y(e);
y_continuous.up(e) = UB_y(e);

equation eq_nao_linear;
eq_nao_linear(i,j)$QT(i,j)..
         w(i,j) =e= x(i)*x(j);


equation fixa_var_primeiro_estagio;
parameter valor_fixado(i);
fixa_var_primeiro_estagio(i)$iVc(i)..
         x(i) =e= valor_fixado(i);

parameter direcao_valor_fixado(i,sc);

model m /
R_FO
*Main
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod
MC_1
MC_2
MC_3
MC_4
/;
*solve m min FO us MIP;

model Ori /
R_FO_Or
Main_Or
/;

*model Fixed /
*R_FO_fixed
*Main_fixed
*/;

Equation
FO_bound_UB
FO_bound_LB
FO_Bound
;

model Rel /
R_FO_rel
Main_L
Main_E
Main_G
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod
MC_1
MC_2
MC_3
MC_4

*MC_w_1
*MC_w_2
*MC_w_3
*MC_w_4
*FO_bound_UB
*FO_bound_LB
*outer_approximation
/;

model Rel_lagrange /
eq_extra_multiplicador
R_FO_rel_lagrange
Main_L
Main_E
Main_G
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod
MC_1
MC_2
MC_3
MC_4
/;

model m_fixa_primeiro_estagio /
R_FO_or
Main_L_ori
Main_E_ori
Main_G_ori
fixa_var_primeiro_estagio
*eq_nao_linear
/;

equations
R_FO_fixed
Main_L_fixed
Main_E_fixed
Main_G_fixed
;


R_FO_fixed(r)$(re_FO(r))..
         (sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
         + sum(e,b_int(r,e)*y_fixed(e))
)
/card(sc)
          =E= FO;

Main_L_fixed(r)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y_fixed(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =L= 0;
Main_E_fixed(r)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y_fixed(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =E= 0;
Main_G_fixed(r)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i))
+ sum(e,b_int(r,e)*y_fixed(e))
+ c(r) + sum(sc$scElemento(sc), cStochastic(r, sc)) =G= 0;

model m_fixa_primeiro_estagio_fixed /
R_FO_fixed
Main_L_fixed
Main_E_fixed
Main_G_fixed
fixa_var_primeiro_estagio
*eq_nao_linear
/;


Variable FO_Bound_variavel;

Parameter Bound_UB;
Parameter Escolhe_Variavel(i);
Parameter novo_LB(i), novo_UB(i);

FO_bound_UB..
         FO =L= UB;

FO_bound_LB..
         FO =G= LB;

FO_Bound..
         FO_Bound_variavel =E= sum(i, Escolhe_Variavel(i)*x(i));

model Rel_bound /
R_FO_rel_continuous
Main_continuous
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod
MC_1
MC_2
MC_3
MC_4
FO_bound_UB
FO_bound_LB
FO_Bound
*outer_approximation
/
;



** metododo plano de cortes
* metodo do plano de cortes
scalar tamanho_passo / 10/;
set cuts /1*1000/;
set cuts_used(cuts);
scalar ultimo_corte /0/;
variable new_multiplier(i,sc);
parameter cut_a(cuts,i,sc);
parameter cut_b(cuts);
*define modello de plano de cortes
variable cuts_Fo;
Equations cuts_Cons(cuts);
cuts_Cons(cuts)$cuts_used(cuts).. cuts_Fo =G= sum((i,sc)$dominio_multip(i,sc), cut_a(cuts,i,sc)*new_multiplier(i,sc)) + cut_b(cuts);
model planning_cuts /cuts_Cons/;



******* bundle_method (metodo dos feixes

scalar bundle_tamanho_maximo /100/;
scalar bundle_ultimo_feixe_usado /0/;
scalar bundle_sigma;
bundle_sigma := 0.01;
*(1/10)**6;
set bundle_fe /1*1000/;
alias(bundle_fe, bundle_fe_linha);
set bundle_fe_usado(bundle_fe);
scalar bundle_fe_ultima /0/;
bundle_fe_usado(bundle_fe) := no;
parameter bundle_x(i,sc), bundle_x_novo(i,sc), bundle_y(i,sc), bundle_y_novo(i,sc), bundle_z(i,sc);
* x e z multiplicadores
* y subgradient
scalar bundle_gamma /1/;
scalar bundle_fx, bundle_fx_novo, bundle_fz, bundle_e;
scalar bundle_delta
scalar passo_serio
scalar bundle_i1, bundle_i2;
scalar bundle_iterations /0/;

parameter bundle_v(bundle_fe), bundle_epsilon, bundle_d(i,sc);

variable bundle_new_multiplier(i,sc);
variable bundle_FO;
parameter feixe_1(bundle_fe,i,sc), feixe_2(bundle_fe);
scalar bundle_erro;

Equation bundle_con, bundle_con_fo;

* modelo de minimizacao do metodo de feixe

variable bundle_new_multiplier(i,sc);
variable bundle_FO, bundle_psi;
Equation bundle_con, bundle_con_fo;

bundle_con(bundle_fe)$bundle_fe_usado(bundle_fe)..
         bundle_psi =G= bundle_fx - feixe_2(bundle_fe) + sum((i,sc)$dominio_multip(i,sc), feixe_1(bundle_fe,i,sc)*(bundle_new_multiplier(i,sc)-bundle_x(i,sc)));

bundle_con_fo..
         bundle_FO =E= bundle_psi + bundle_gamma/2 * sum((i,sc)$dominio_multip(i,sc), power((bundle_new_multiplier(i,sc)-bundle_x(i,sc)),2));

model m_bundle /bundle_con, bundle_con_fo/;



*** bunsca local
set busca_local(*) / 1*1 /;
scalar local_viavel /1/;
scalar local_viavel_busca /1/;
scalar LB_iteracao_local, LB_busca;
Option Seed=3141;



*****
*without decomposing (ds)
***
Variable x_ws(i,sc);
integer Variable y_ws(e,sc);
Variable y_continuous_ws(e,sc);
Variable w_ws(i,j,sc);
Binary Variable z_ws(j,l,sc);
Variable hat_x_ws(i,j,l,sc);
Positive Variable delta_x_ws(j,sc);
Variable v_ws(i,j,sc);
variable FO;

parameter y_fixed_ws
                 x_relaxado_ws(i,sc)
                 w_relaxado_ws(i,j,sc)
                 v_relaxado_ws(i,j,sc)
                 delta_x_relaxado_ws(i,sc);

Equations
R_FO_Or_ws
R_FO_rel_ws
R_FO_rel_lagrange_ws
R_FO_ws
Main_L_ws
Main_E_ws
Main_G_ws
Main_L_ori_ws
Main_E_ori_ws
Main_G_ori_ws
R_FO_fixed_ws
Main_fixed_ws
Main_L_fixed_ws
Main_E_fixed_ws
Main_G_fixed_ws
Def_w_ws
Def_x_ws
UB_hat_x_ws
LB_hat_x_ws
UB_prod_ws
LB_prod_ws
MC_1_ws
MC_2_ws
MC_3_ws
MC_4_ws
Main_continuous_ws
R_FO_Rel_continuous_ws

MC_w_1_ws
MC_w_2_ws
MC_w_3_ws
MC_w_4_ws

linking_ws
;

* define bounds
x_ws.lo(i,sc) = LB_x(i);
x_ws.up(i,sc)$UB_x(i) = UB_x(i);
y_ws.lo(e,sc) = LB_y(e);
y_ws.up(e,sc) = UB_y(e);
y_continuous_ws.lo(e,sc) = LB_y(e);
y_continuous_ws.up(e,sc) = UB_y(e);

****

linking_ws(i,sc)$(iVc(i))..
         x_ws(i,sc) =e= sum(sc1, prob(sc1)*x_ws(i,sc1));

R_FO_rel_ws(r)$(re_FO(r))..
         sum((i,j,sc)$BL(r,i,j),a(r,i,j)*prob(sc)*w_ws(i,j,sc)) + sum((i,sc),b(r,i)*prob(sc)*x_ws(i,sc)) + sum((e,sc),b_int(r,e)*prob(sc)*y_ws(e,sc)) =E= FO;

R_FO_rel_lagrange_ws(r)$(re_FO(r))..
         sum((i,j,sc)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum((i,sc),b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y(e))
         + sum((i,sc)$iVc(i),(multiplicadores(i,sc)$(ord(sc)<card(sc))-multiplicadores(i,sc-1))*x_ws(i,sc))
          =E= FO;

R_FO_rel_continuous_ws(r)$(re_FO(r))..
         sum((i,j,sc)$BL(r,i,j),a(r,i,j)*prob(sc)*w_ws(i,j,sc)) + sum((i,sc),b(r,i)*prob(sc)*x_ws(i,sc)) + sum((e,sc),b_int(r,e)*prob(sc)*y_continuous_ws(e,sc)) =E= FO;


R_FO_fixed_ws(r)$(re_FO(r))..
         sum((i,j,sc)$BL(r,i,j),a(r,i,j)*prob(sc)*x_ws(i,sc)*x_ws(j,sc)) + sum((i,sc),b(r,i)*prob(sc)*x_ws(i,sc)) + sum((e,sc),b_int(r,e)*prob(sc)*y_fixed_ws(e,sc)) =E= FO;

Main_L_fixed_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc)) + sum(e,b_int(r,e)*y_fixed_ws(e,sc))
         + c(r) + cStochastic(r, sc)=L= 0;

Main_E_fixed_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc)) + sum(e,b_int(r,e)*y_fixed_ws(e,sc))
         + c(r) + cStochastic(r, sc) =E= 0;


Main_G_fixed_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc)) + sum(e,b_int(r,e)*y_fixed_ws(e,sc))
         + c(r) + cStochastic(r, sc) =G= 0;

*Main_Or_ws(r,sc)$(not re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc)) + sum(e,b_int(r,e)*y_ws(e,sc)) + c(r) + cStochastic(r, sc) =L= 0;

*FO_LB_ws..
*         FO =G= LB;

*FO_UB_ws..
*         FO =L= UB;

*parameter mult(r);
mult(r) = 0;

R_FO_ws(r)$(re_FO(r))..
         FO =e= sum((i,j,sc)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum((i,sc),b(r,i)*x_ws(i,sc));
* + sum(r1$(not re_FO(r1)), mult(r1)*( sum((i,j)$BL(r1,i,j),a(r1,i,j)*w_ws(i,j,sc)) + sum(i,b(r1,i)*x(i,sc)) + c(r1))) =E= FO;

Main_L_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + cStochastic(r, sc) =L= 0;
Main_E_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + cStochastic(r, sc) =E= 0;
Main_G_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y(e))
+ c(r) + cStochastic(r, sc) =G= 0;

****

R_FO_Or_ws(r)$(re_FO(r))..
         sum(sc, prob(sc)*[sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum((i),b(r,i)*x_ws(i,sc))
         + sum(e,b_int(r,e)*y_ws(e,sc))])
          =E= FO;

Main_L_ori_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=-1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y_ws(e,sc))
+ c(r) + cStochastic(r, sc) =L= 0;

Main_E_ori_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=0)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y_ws(e,sc))
+ c(r) + cStochastic(r, sc) =E= 0;

Main_G_ori_ws(r,sc)$(not re_FO(r) and tipo_restricao(r)=1)..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x_ws(i,sc)*x_ws(j,sc)) + sum(i,b(r,i)*x_ws(i,sc))
+ sum(e,b_int(r,e)*y_ws(e,sc))
+ c(r) + cStochastic(r, sc) =G= 0;

****

Main_continuous_ws(r,sc)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w_ws(i,j,sc)) + sum(i,b(r,i)*x_ws(i,sc)) + sum(e,b_int(r,e)*y_continuous_ws(e,sc)) + c(r) =L= 0;

Def_w_ws(i,j,sc)$QT(i,j)..
         w_ws(i,j,sc) =E= (x_ws.up(j,sc)-x_ws.lo(j,sc))*sum(l$(ord(l)<=-pp_ws(j,sc)),2**(-ord(l))*hat_x_ws(i,j,l,sc)) + (x_ws.up(j,sc) - x_ws.lo(j,sc))*v_ws(i,j,sc) + x_ws(i,sc)*x_ws.lo(j,sc);

Def_x_ws(j,sc)$DS(j)..
         x_ws(j,sc) =E= (x_ws.up(j,sc) - x_ws.lo(j,sc))*sum(l$(ord(l)<=-pp_ws(j,sc)),2**(-ord(l))*z_ws(j,l,sc)) + (x_ws.up(j,sc) - x_ws.lo(j,sc))*delta_x_ws(j,sc) + x_ws.lo(j,sc);

UB_hat_x_ws(i,j,l,sc)$(QT(i,j) and ord(l)<=-pp_ws(j,sc))..
         hat_x_ws(i,j,l,sc) =L= x_ws.up(i,sc)*z_ws(j,l,sc);

LB_hat_x_ws(i,j,l,sc)$(QT(i,j) and ord(l)<=-pp_ws(j,sc))..
         hat_x_ws(i,j,l,sc) =G= x_ws.lo(i,sc)*z_ws(j,l,sc);

UB_prod_ws(i,j,l,sc)$(QT(i,j) and ord(l)<=-pp_ws(j,sc))..
         x_ws(i,sc) - hat_x_ws(i,j,l,sc) =L= x_ws.up(i,sc)*(1-z_ws(j,l,sc));

LB_prod_ws(i,j,l,sc)$(QT(i,j) and ord(l)<=-pp_ws(j,sc))..
         x_ws(i,sc) - hat_x_ws(i,j,l,sc) =G= x_ws.lo(i,sc)*(1-z_ws(j,l,sc));

MC_1_ws(i,j,sc)$QT(i,j)..
         v_ws(i,j,sc) =G= delta_x_ws(j,sc)*x_ws.lo(i,sc);

MC_2_ws(i,j,sc)$QT(i,j)..
         v_ws(i,j,sc) =L= delta_x_ws(j,sc)*x_ws.up(i,sc);

MC_3_ws(i,j,sc)$QT(i,j)..
         v_ws(i,j,sc) =G= (2**pp_ws(j,sc))*(x_ws(i,sc)-x_ws.up(i,sc)) + delta_x_ws(j,sc)*x_ws.up(i,sc);

MC_4_ws(i,j,sc)$QT(i,j)..
         v_ws(i,j,sc) =L= (2**pp_ws(j,sc))*(x_ws(i,sc)-x_ws.lo(i,sc)) + delta_x_ws(j,sc)*x_ws.lo(i,sc);

MC_w_1_ws(i,j,sc)$QT(i,j)..
         w_ws(i,j,sc) =G= x_ws(i,sc)*x_ws.up(j,sc) + x_ws.up(i,sc)*x_ws(j,sc) - x_ws.up(i,sc)*x_ws.up(j,sc);
MC_w_2_ws(i,j,sc)$QT(i,j)..
         w_ws(i,j,sc) =G= x_ws(i,sc)*x_ws.lo(j,sc) + x_ws.lo(i,sc)*x_ws(j,sc) - x_ws.lo(i,sc)*x_ws.lo(j,sc);
MC_w_3_ws(i,j,sc)$QT(i,j)..
         w_ws(i,j,sc) =L= x_ws(i,sc)*x_ws.lo(j,sc) + x_ws.up(i,sc)*x_ws(j,sc) - x_ws.up(i,sc)*x_ws.lo(j,sc);
MC_w_4_ws(i,j,sc)$QT(i,j)..
         w_ws(i,j,sc) =L= x_ws(i,sc)*x_ws.up(j,sc) + x_ws.lo(i,sc)*x_ws(j,sc) - x_ws.lo(i,sc)*x_ws.up(j,sc);

model Rel_ws /
R_FO_rel_ws
Main_L_ws
Main_E_ws
Main_G_ws
Def_w_ws
Def_x_ws
UB_hat_x_ws
LB_hat_x_ws
UB_prod_ws
LB_prod_ws
MC_1_ws
MC_2_ws
MC_3_ws
MC_4_ws
linking_ws

*MC_w_1
*MC_w_2
*MC_w_3
*MC_w_4
*FO_bound_UB
*FO_bound_LB
*outer_approximation
/;

model Fixed_ws /
R_FO_fixed_ws
Main_L_fixed_ws
Main_E_fixed_ws
Main_G_fixed_ws
linking_ws
/;

model m_ws /
R_FO_Or_ws
Main_L_ori_ws
Main_E_ori_ws
Main_G_ori_ws
linking_ws
/;

****
