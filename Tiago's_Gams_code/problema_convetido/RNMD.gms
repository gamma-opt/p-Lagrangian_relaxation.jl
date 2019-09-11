$offorder
$offdigit
$onempty
option limrow = 1000;
option limcol = 1000;
*set e(*) empty set //;

* scalares de tempo
Scalar tcomp, texec, telapsed;


* Indices das variáveis x
*Sets i;
$include "C:\Artigos\p-LD\numerical experiments\v7\problema_convetido\teste_sets.inc";
alias(i,j);
alias (i,ii);

* Indice de precisão
Set l /1*100/;
alias(l,ll);

* Indice de funções
*Set r;
alias(r,r1,r2);
SET re_FO(r) /0/;

* Incide representando os pares i,j de variáveis que aparecem nas funções r
Set BL(r,i,j);


scalar mudanca_bound_variavel;

Scalar p;
p = -20;

Parameter LB_x(i), UB_x(i);
Parameter LB_y(e), UB_y(e);


Parameter a(r,i,j);




Parameter b(r,i);

Parameter b_int(r, e);

Parameter y_fixed(e);

Parameter c(r);


display r, l;

b_int(r, e) = 0;
LB_y(e) = 0;
UB_y(e) = 0;
$include "C:\Artigos\p-LD\numerical experiments\v7\problema_convetido\teste.inc";
*c(r)$(not re_FO(r)) := c(r) - 10**(-7)

* Selecionar quais variáveis serão expandidas
***
Binary Variable seleciona_variavel_expandida(i);
Variable FO_seleciona_variavel_expandida;

Equation
R_FO_seleciona_variavel_expandida
R_seleciona_variavel_expandida
;

R_FO_seleciona_variavel_expandida..
         FO_seleciona_variavel_expandida =E= sum(i, seleciona_variavel_expandida(i));
;

R_seleciona_variavel_expandida(i,j)$sum(r, BL(r,i,j))..
         1 =L= seleciona_variavel_expandida(i) + seleciona_variavel_expandida(j);

model m_seleciona_variavel_expandida
/
R_FO_seleciona_variavel_expandida
R_seleciona_variavel_expandida
/;

option optca = 0;
option optcr = 0;
*solve m_seleciona_variavel_expandida min FO_seleciona_variavel_expandida us MIP;
scalar default_selectiona_variavel_expandida;
default_selectiona_variavel_expandida = sum(j$sum((r,i), BL(r,i,j)), 1);
*display default_selectiona_variavel_expandida, BL;
***

Variable x(i);
integer Variable y(e);
Variable y_continuous(e);
Variable w(i,j);
Binary Variable z(j,l);
Variable hat_x(i,j,l);
Positive Variable delta_x(j);
Variable v(i,j);
variable FO;

Equations
R_FO
Main
R_FO_fixed
Main_fixed
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
;

Equations
R_FO_Or
R_FO_rel
Main_Or
;

R_FO_rel(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) =E= FO;

R_FO_rel_continuous(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_continuous(e)) =E= FO;

R_FO_Or(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) =E= FO;

R_FO_fixed(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_fixed(e)) =E= FO;

Main_fixed(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_fixed(e)) + c(r) =L= 0;

Main_Or(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) + c(r) =L= 0;

scalar UB;
scalar LB;

FO_LB..
         FO =G= LB;

FO_UB..
         FO =L= UB;

parameter mult(r);
mult(r) = 0;

R_FO(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(r1$(not re_FO(r1)), mult(r1)*( sum((i,j)$BL(r1,i,j),a(r1,i,j)*w(i,j)) + sum(i,b(r1,i)*x(i)) + c(r1))) =E= FO;

Main(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y(e)) + c(r) =L= 0;

Main_continuous(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(e,b_int(r,e)*y_continuous(e)) + c(r) =L= 0;

Def_w(i,j)$sum(r,BL(r,i,j))..
         w(i,j) =E= (x.up(j)-x.lo(j))*sum(l$(ord(l)<=-p),2**(-ord(l))*hat_x(i,j,l)) + (x.up(j) - x.lo(j))*v(i,j) + x(i)*x.lo(j);

Def_x(j)$sum((r,i),BL(r,i,j))..
         x(j) =E= (x.up(j) - x.lo(j))*sum(l$(ord(l)<=-p),2**(-ord(l))*z(j,l)) + (x.up(j) - x.lo(j))*delta_x(j) + x.lo(j);

UB_hat_x(i,j,l)$(sum(r,BL(r,i,j)) and ord(l)<=-p)..
         hat_x(i,j,l) =L= x.up(i)*z(j,l);

LB_hat_x(i,j,l)$(sum(r,BL(r,i,j)) and ord(l)<=-p)..
         hat_x(i,j,l) =G= x.lo(i)*z(j,l);

UB_prod(i,j,l)$(sum(r,BL(r,i,j)) and ord(l)<=-p)..
         x(i) - hat_x(i,j,l) =L= x.up(i)*(1-z(j,l));

LB_prod(i,j,l)$(sum(r,BL(r,i,j)) and ord(l)<=-p)..
         x(i) - hat_x(i,j,l) =G= x.lo(i)*(1-z(j,l));

MC_1(i,j)$sum(r,BL(r,i,j))..
         v(i,j) =G= delta_x(j)*x.lo(i);

MC_2(i,j)$sum(r,BL(r,i,j))..
         v(i,j) =L= delta_x(j)*x.up(i);

MC_3(i,j)$sum(r,BL(r,i,j))..
         v(i,j) =G= (2**p)*(x(i)-x.up(i)) + delta_x(j)*x.up(i);

MC_4(i,j)$sum(r,BL(r,i,j))..
         v(i,j) =L= (2**p)*(x(i)-x.lo(i)) + delta_x(j)*x.lo(i);

parameter x_value(i,l);
scalar pp /-50/;
x_value(i,l) := LB_x(i) + (UB_x(i) - LB_x(i))*sum(ll$(ord(ll)<ord(l)), 2**(pp));
outer_approximation(i,l)$sum(r, BL(r,i,i) and ord(l)<=-pp+2)..
         w(i,i) =G= 2*x_value(i,l)*x(i) - power(x_value(i,l), 2);

* Bounds das variáveis
delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
x.lo(i) = LB_x(i);
x.up(i)$UB_x(i) = UB_x(i);
y.lo(e) = LB_y(e);
y.up(e) = UB_y(e);
y_continuous.lo(e) = LB_y(e);
y_continuous.up(e) = UB_y(e);


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

model Fixed /
R_FO_fixed
Main_fixed
/;

Equation
FO_bound_UB
FO_bound_LB
FO_Bound
;

model Rel /
R_FO_rel
Main
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
*FO_bound_UB
*FO_bound_LB
*outer_approximation
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
outer_approximation
/

scalar best /-1000/;
scalar iter /0/;
scalar step /1/;
Parameter melhor_multiplicador(r);
*7,28
scalar stop /0/;
scalar mod /0/;

parameter subg(r);
parameter mult_old(r);
Parameter s_subg(r);
s_subg(r) = 0;
*subg(r)$(ord(r)<card(r)) = 10;

scalar p_aux;
scalar Iter_O /0/;
scalar Iter_i /0/;
scalar stop_2 /0/;
scalar iter_tot /0/;
UB = 10**8;
LB = -UB;
scalar FO_atual;
*option MIP = OSICPLEX;
*option MIP = OSIGUROBI;
option MIP = CPLEX;
*option NLP = minos;
option NLP = conopt;
*option nlp = knitro;
option optca = 0.0001;
option optcr = 0.000001;
stop = 0;


* procedimento de bound contraction
$ontext
option RMIP = CPLEX;
p = 0;
delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
solve Rel us RMIP min FO;
LB = FO.l;
y_fixed(e) = max(min(round(y.l(e)),UB_y(e)),LB_y(e));
solve Fixed us NLP min FO;
Bound_UB = FO.l;
UB = FO.l;
p = -4;
delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
$include "C:\Artigos\p-LD\numerical experiments\v7\problema_convetido\bound_contraction.gms";
*LB_x(i) =  max(novo_LB(i)-0.001, x.lo(i));
*UB_x(i) =  min(novo_UB(i)+0.001, x.up(i));
*x.lo(i) = LB_x(i);
*x.up(i) = UB_x(i);
$offtext

scalar time_Start, time_current, time_limit;
time_limit = 7200.0;
time_Start = TimeElapsed;

*z.prior(i,l) = ord(l);
*Rel.prioropt = 1;


$onecho > Cplex.opt
mipstart 1
$offecho
rel.OptFile = 0;

p = 0;
scalar aux_x;
file results_parcial  /C:\Users\tiago.andrade\Desktop\RNMDT_Julia\WasterWater result\RNMDT.txt/;

put results_parcial;
         put 'RNMDT' /;
putclose;
results_parcial.ap = 0;

while(not stop,
         delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
*         x_value(i,l) := LB_x(i) + (UB_x(i) - LB_x(i))*sum(ll$(ord(ll)<ord(l)), 2**(p));
         iter_O = iter_O + 1;
$ontext
         loop(i,
                 aux_x = (x.l(i)-x.lo(i))/(x.up(i)-x.lo(i));
                 loop(l$(ord(l)<=-p),
                         if(aux_x >= 2**(-ord(l)),
                                 z.l(i,l) = 1;
                                 aux_x = aux_x - 2**(-ord(l));
                         else
                                 z.l(i,l) = 0;
                         );
                 );
                 delta_x.l(i) = aux_x;
         );
         w.l(i,j) = x.l(i)*x.l(j);
         hat_x.l(i,j,l) = x.l(j)*z.l(i,l);
         v.l(i,j) = 0;
         display x.l, z.l, w.l, y.l, delta_x.l;
$offtext

         time_current = TimeElapsed - time_Start;
                 if(time_current < time_limit,
                         Rel.reslim = time_limit - time_current;
                         Fixed.reslim =  time_limit - time_current;
                 else
                         Rel.reslim = 0;
                         Fixed.reslim = 0;
                 );

                 iter_tot = iter_tot + 1;
*                 Rel.cutoff = UB;
                 solve Rel min FO us MIP;
                 if(p=0 and card(e) = 0,
                         FO_atual = FO.l;
                 else
                         FO_atual = Rel.objest;
                 );
                 LB = max(LB,FO_atual);
                 display FO_atual;
* fixa o valor das variáveis inteiras originais do problema
                 y_fixed(e) = y.l(e);

                 solve Fixed min FO us NLP;
                 if(Fixed.modelstat <= 2,
                         UB = min(UB,FO.l);
                 );
                 display UB, LB, p, y_fixed;

                 Bound_UB = UB;
*$include "C:\Artigos\p-LD\numerical experiments\v7\problema_convetido\bound_contraction.gms";
         telapsed = TimeElapsed;
         results_parcial.ap = 1;
         put results_parcial;
                 put p /;
                 put telapsed /;
                 put LB:10:4 /;
                 put UB:10:4 / /;
         putclose;



         if((UB - LB < 10**(-3)) or (UB-LB)/abs(LB) < 10**(-6) or TimeElapsed > time_limit,
                 stop = 1;
         else
                 p = p - 1;
         );
);
*display mult, subg, mod, LB, UB, p, iter_O, iter_i;
display UB, LB, p, iter_i, iter_tot, iter_O, iter_i, telapsed;

file results /C:\Artigos\p-LD\numerical experiments\v7\problema_convetido\results.csv/;
put results;
         put telapsed /;
         put LB:10:4 /;
         put UB:10:4 /;
         put iter_i /;
         put iter_o /;
         put iter_tot /;
         put p;
putclose;

display e, b_int, y.l;
