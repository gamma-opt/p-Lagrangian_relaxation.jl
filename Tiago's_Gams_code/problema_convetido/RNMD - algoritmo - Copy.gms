$offorder
$offdigit
$onempty
option limrow = 1000;
option limcol = 1000;
*set e(*) empty set //;

File Results /'D:\p-LD\numerical experiments\Results 20180118\Results.csv'/;
File Results_iterations /'D:\p-LD\numerical experiments\Results 20180118\Results_iterations.csv'/;
File Results_feixes /'D:\p-LD\numerical experiments\Results 20180118\Results_feixes.csv'/;

* scalares de tempo
Scalar tcomp, texec, telapsed;


* Indices das variáveis x
*Sets i;
$include "D:\p-LD\numerical experiments\v9\problema_convetido\teste_sets.inc";
*base
*$include "D:\p-LD\\numerical experiments\\Results 20180118\\model_size.inc";
*Set i / 1* 10/;
*Set e(*) / /;
*Set r / 1* 20/;
*Set sc /1 * 20/;

alias(i,j);
alias (i,ii);
alias(sc,sc1);
*variavel de primeiro estagio
Set iVc(i);

* Indice de precisão
Set l /1*1000/;
alias(l,ll);

* Indice de funções
*Set r;
alias(r,r1,r2);
Set re_FO(r);
re_FO('1') := yes;
* tipo resticao (-1 <=, 0 ==, 1 >=)
parameter tipo_restricao(r);





* Incide representando os pares i,j de variáveis que aparecem nas funções r
Set BL(r,i,j);
Set QT(i,j);
Set DS(j);


scalar UB_iteracao, LB_iteracao;
parameter x_stoc_rel(i,sc);
parameter y_stoc_rel(e,sc);
parameter x_stoc_ori(i,sc);


Scalar p;
p = 0;
parameter pp(i);
pp(i) := p;

Parameter LB_x(i), UB_x(i);
Parameter LB_y(e), UB_y(e);


Parameter a(r,i,j);




Parameter b(r,i);

Parameter b_int(r, e);
b_int(r, e) := 0;

Parameter y_fixed(e);

Parameter c(r);

*set sc /1*29/;
set scElemento(sc);
parameter cStochastic(r,sc);
cStochastic(r,sc) := 0;
parameter prob(sc);


display r, l;

b_int(r, e) = 0;
LB_y(e) = 0;
UB_y(e) = 0;
$include "D:\p-LD\numerical experiments\v9\problema_convetido\test.inc";
$ontext
*sorteio
re_FO('1') := yes;
UB_x(i) := 1;
LB_x(i) := -1;

b(r,i) := uniform(-1,1);
c(r)$(not re_FO(r)) := uniform(-3,0);
iVc(i)$(ord(i) <= card(i)/10) := yes;
*iVc('1') := yes;
*iVc('2') := yes;
*iVc('3') := yes;
a(r,i,j)$(ord(i)<ord(j) and not iVc(i) and not iVc(j)) := uniform(0,1);
UB_y(e) := 1;
LB_y(e) := -1;
b_int(r,e) := uniform(-1,1);
prob(sc) := 1/card(sc);
cStochastic(r,sc) := uniform(0.,1.);
tipo_restricao(r)$(not re_FO(r)) := -1;
$offtext


****
* parametros referentes ao dual (lagrange)
* only for iVc(i) and sc < card(sc)
parameter subgradient(i,sc);
parameter multiplicadores(i,sc);
multiplicadores(i,sc) := 0;
set dominio_multip(i,sc);
dominio_multip(i,sc)$iVc(i) := yes;
****

BL(r,i,j) := (a(r,i,j) <> 0 and ord(i)<=ord(j));
QT(i,j) = sum(r,BL(r,i,j) and ord(i) <= ord(j));
DS(j) = sum(i, QT(i,j));
*c(r)$(not re_FO(r)) := c(r) - 10**(-7)

scalar UB, LB;

$include "D:\p-LD\numerical experiments\v9\problema_convetido\DefineModelos.gms"


scalar best /-1000/;
scalar iter /0/;
scalar stop /0/;

scalar Iter_O /0/;
scalar Iter_i /0/;
scalar stop_2 /0/;
scalar iter_tot /0/;
UB = 10**13;
LB = -UB;
scalar FO_atual;

*option MIP = OSICPLEX;
*option MIP = OSIGUROBI;
*option MIP = CPLEX;
option NLP = minos;
*option NLP = conopt;
*option nlp = knitro;
*option QCP = knitro;
option QCP = conopt;
option optca = 0.000001;
option optcr = 0.00000001;

*option optca = 0.0001;
*option optcr = 0.000001;
stop = 0;
p = 0;
pp(j) = p;
*$include "D:\p-LD\numerical experiments\v9\problema_convetido\bound_contraction.gms"

* dfine limite de tempo
scalar time_Start, time_current, time_limit;
time_limit = 2*3600.0;
time_Start = TimeElapsed;

time_current = TimeElapsed - time_Start;
if(time_current < time_limit,
         Rel_lagrange.reslim = time_limit - time_current;
else
         Rel_lagrange.reslim = 0;
);


scalar subdradient_modulo;
scalar iteracao_grande /0/;
put Results;
put 'Resulsts' / /;
putclose;
Results.ap := 1;
put Results_iterations;
put Results_iterations;
         put 'Iteration ', ' ',
         put 'Time_elapsed ', ' ',
         put 'UB_iteracao ', ' ',
         put 'LB_iteracao ', ' ',
         put 'p ';
putclose;
Results_iterations.ap := 1;
subdradient_modulo := 1;
put Results_feixes;
         put 'Grande_iteracao ', ' ',
         put 'Bundle_iteracao ', ' ',
         put 'feixe_usado ', ' ',
         put 'feixe_1 ', ' ',
         put 'feixeis_2 ', ' ',
         put 'v ', ' ',
         put 'd ', ' ',
         put 'epsilon ', ' ',
         put 'bundle_erro ',
;
putclose;
Results_feixes.ap := 1;



scalar UB_externo, LB_externo;
UB_externo = UB;
LB_externo = LB;
p = 0;
pp(j) = p;
pp_ws(j,sc) = p;



$ontext
while((time_current < time_limit and iteracao_grande < 100 and UB_externo - LB_externo > 10**(-2) and UB_externo - LB_externo > 10**(-4)*LB_externo),

*$include "D:\\p-LD\\numerical experiments\\v9\\problema_convetido\\subgradient_metod.gms"

*$include "D:\\p-LD\\numerical experiments\\v9\\problema_convetido\\cutting_planes_method.gms"

$include "D:\\p-LD\\numerical experiments\\v9\\problema_convetido\\bundle_method.gms"
         iteracao_grande =  iteracao_grande + 1;
         UB_externo = UB;
         LB_externo = LB;

         time_current = TimeElapsed - time_Start;
         put Results_iterations;
                 put iteracao_grande, ' ',
                  time_current, ' ',
                 UB_externo:10:4, ' ',
                 LB_externo:10:4, ' ',
                 p;
         putclose;

         p = p - 1;
         pp(j) = p;
);
$offtext

b_int(r,e) = 0;
parameter stop, stop_2, maior_valor;
set maior_j(j);
$include "D:\\p-LD\\numerical experiments\\v9\\problema_convetido\\algRNMDT.gms"


put Results_iterations;
         put iteracao_grande, ' ',
         time_current, ' ',
         UB_externo:10:4, ' ',
         LB_externo:10:4, ' ',
         p;
putclose;

display time_current;

display p, UB_externo, LB_externo;

*delta_x_ws.up(j,sc)$sum((r,i),BL(r,i,j)) := 2**pp(j);
*solve Rel_ws us MIP max FO;
*y_fixed_ws(e,sc) = y_ws.l(e,sc);
*solve Fixed_ws max FO us NLP;
*display cStochastic, iVc;

display
re_FO
BL
a
b
b_int;


display
dominio_multip
iVc;
