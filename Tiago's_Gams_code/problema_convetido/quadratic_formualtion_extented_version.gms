$offorder

* scalares de tempo
Scalar tcomp, texec, telapsed;


* Indices das variáveis x
*Sets i;
$include "C:\Users\tiago.andrade\Desktop\problems QCQP\problems QCQP\problema_convetido\teste_sets.inc";
alias(i,j);

* Indice de precisão
Set l /1*100/;
alias(l,ll);

* Indice de funções
*Set r;
alias(r,r1,r2);
SET re_FO(r) /0/;

* Incide representando os pares i,j de variáveis que aparecem nas funções r
Set BL(r,i,j);




Scalar p;
p = -20;

Parameter LB_x(i), UB_x(i);


Parameter a(r,i,j);




Parameter b(r,i);


Parameter c(r);


Parameter r_dsdp(r,i);

display r, l;

$include "C:\Users\tiago.andrade\Desktop\problems QCQP\problems QCQP\problema_convetido\teste.inc";

$include "C:\Users\tiago.andrade\Desktop\problems QCQP\problems QCQP\problema_convetido\r_DSDP.inc";

Positive Variable x(i);
Positive Variable w(i);
Binary Variable z(i,l);
*z.up(i,l) = 1;
Positive Variable hat_x(i,l);
Positive Variable delta_x(i);
Positive Variable v(i);
positive Variable hat_z(i,l,ll);
Positive Variable hat_delta_x(i,l);

variable FO;

Equations
R_FO
Main
Def_w
Def_x
UB_hat_x
LB_hat_x
UB_prod
LB_prod

UB_delta_hat_x
UB_delta_prod
LB_delta_prod

MC_3
MC_4

MC_1
MC_2

MC_w_1
MC_w_2
MC_w_3

logic_1
logic_2
logic_3
;

Equations
R_FO_Or
R_FO_rel
Main_Or
R_FO_decom
Main_decom
;

*R_FO_rel(r)$(re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i)) + sum(i,b(r,i)*x(i)) =E= FO;

R_FO_Or(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) =E= FO;

R_FO_decom(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + (10)*sum(i, power(x(i),2)) - (10)*sum(i, w(i)) =E= FO;

Main_Or(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + c(r) =L= 0;

Main_decom(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + c(r) + (10)*sum(i, power(x(i),2)) - (10)*sum(i, w(i)) =L= 0;

parameter mult(r);
mult(r) = 0;

*R_FO(r)$(re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(r1$(not re_FO(r1)), mult(r1)*( sum((i,j)$BL(r1,i,j),a(r1,i,j)*w(i,j)) + sum(i,b(r1,i)*x(i)) + c(r1))) =E= FO;

*Main(r)$(not re_FO(r))..
*         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + c(r) =L= 0;

Def_w(i)..
*         w(i) =E= (x.up(i)-x.lo(i))*sum(l$(ord(l)<=-p),2**(-ord(l))*hat_x(i,l)) + (x.up(i)-x.lo(i))*v(i) + x(i)*x.lo(i);
         w(i) =E= 2*(x.up(i)-x.lo(i))*(sum(l$(ord(l)<=-p),2**(-ord(l))*z(i,l)) + delta_x(i))*x.lo(i)
                 + x.lo(i)**2
                 + ((x.up(i)-x.lo(i))**2)*(sum(l$(ord(l)<=-p),4**(-ord(l))*z(i,l)) + v(i) + 2*sum(l$(ord(l)<=-p),2**(-ord(l))*hat_delta_x(i,l)) + 2*sum((l,ll)$(ord(l)<=-p and ord(ll) < ord(l) and ord(l)+ord(ll)<=-p),2**(-ord(l)-ord(ll))*hat_z(i,l,ll)) + 2*sum((l,ll)$(ord(l)<=-p and ord(ll) < ord(l) and ord(l)+ord(ll)>-p),2**(-ord(l)-ord(ll))))
;

Def_x(i)..
         x(i) =E= (x.up(i)-x.lo(i))*sum(l$(ord(l)<=-p),2**(-ord(l))*z(i,l)) + (x.up(i)-x.lo(i))*delta_x(i) + x.lo(i);

MC_1(i)..
         v(i) =L= (2**p)*delta_x(i);

MC_2(i)..
         v(i) =G= 2*(2**p)*delta_x(i) - (4**p);

MC_w_1(i)..
         w(i) =L= x.up(i)*x(i) + x.lo(i)*x(i) - x.up(i)*x.lo(i);

MC_w_2(i)..
         w(i) =G= 2*x.up(i)*x(i) - x.up(i)**2;

MC_w_3(i)..
         w(i) =G= 2*x.lo(i)*x(i) - x.lo(i)**2;

logic_1(i,l,ll)$(ord(l)<=-p and ord(ll) < ord(l) and ord(l)+ord(ll)<=-p)..
         hat_z(i,l,ll) =L= z(i,l);

logic_2(i,l,ll)$(ord(l)<=-p and ord(ll) < ord(l) and ord(l)+ord(ll)<=-p)..
         hat_z(i,l,ll) =L= z(i,ll);

logic_3(i,l,ll)$(ord(l)<=-p and ord(ll) < ord(l) and ord(l)+ord(ll)<=-p)..
         hat_z(i,l,ll) =G= z(i,l) + z(i,ll) - 1;

UB_delta_hat_x(i,l)$(ord(l) <= -p)..
         hat_delta_x(i,l) =L= (2**p)*z(i,l);

UB_delta_prod(i,l)$(ord(l) <= -p)..
         delta_x(i) - hat_delta_x(i,l) =L= (2**p)*(1-z(i,l));

LB_delta_prod(i,l)$(ord(l) <= -p)..
          hat_delta_x(i,l) =L= delta_x(i);




* Bounds das variáveis
delta_x.up(i) = 2**p;
x.lo(i) = LB_x(i);
x.up(i) = UB_x(i);

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
/;

model Decom /
R_FO_decom
Main_decom
Def_w
Def_x


MC_1
MC_2


MC_w_1
*MC_w_2
*MC_w_3

logic_1
logic_2
*logic_3

UB_delta_hat_x
UB_delta_prod
LB_delta_prod
/;

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
subg(r)$(ord(r)<card(r)) = 10;

scalar p_aux;
scalar Iter_O /0/;
scalar Iter_i /0/;
scalar stop_2 /0/;
scalar iter_tot /0/;
scalar UB;
scalar LB;
UB = 10**8;
LB = -UB;
scalar FO_atual;
*option MIP = GUROBI;
option MIP = CPLEX;
option NLP = baron;
option miqcp = cplex;
option qcp = conopt;
option optcr = 0.5;
stop = 0;
p = -5;
delta_x.up(i) = 2**p;

scalar time_Start, time_current, time_limit;
time_limit = 3600.0*24*1;
time_Start = TimeElapsed;
Decom.reslim = time_limit

*solve ori min FO us nlp
solve Decom min FO us miqcp
display FO.l
solve Ori min FO us qcp
display FO.l

*display mult, subg, mod, LB, UB, p, iter_O, iter_i;
*display UB, LB, p, iter_i, iter_tot, iter_O, iter_i, telapsed;

$ontext
file results /C:\Users\TA\Desktop\problems QCQP\problema_convetido\results.csv/;
put results;
         put telapsed /;
         put LB:10:4 /;
         put UB:10:4 /;
         put iter_i /;
         put iter_o /;
         put iter_tot /;
         put p;
putclose;
$offtext

$ontext
while(not stop,
         iter = iter + 1;
*         step = 0.1;
         solve m min FO us MIP;

         subg(r1)$(ord(r1)<card(r1)) = sum((i,j)$BL(r1,i,j),a(r1,i,j)*w.l(i,j)) + sum(i,b(r1,i)*x.l(i)) + c(r1);
         s_subg(r1) = 0.3*s_subg(r1) + 0.7*subg(r1);
         mod = (sum(r1$(ord(r1)<card(r1)), s_subg(r1)*s_subg(r1)))**(1/2);
         display mod, subg;
*         subg(r1) = subg(r1)/mod;
         step = (7049.25 - FO.l)/(mod**2);
         display iter, step, mult, mod, FO.l;
         mult(r1)$(ord(r1)<card(r1)) = max(0,mult(r1) + step*s_subg(r1));

         display subg, mult;

         if(FO.l > best,
                 best = FO.l;
                 melhor_multiplicador(r) = mult(r);
         );

         if(iter > 2000,
                 stop = 1;
         );
);

display best, melhor_multiplicador;
$offtext
