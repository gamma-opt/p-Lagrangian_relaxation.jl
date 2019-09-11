*$offorder

* scalares de tempo
Scalar tcomp, texec, telapsed;


* Indices das variáveis x
*Sets i;
$include "C:\Artigos\p-LD\numerical experiments\v2\problema_convetido\teste_sets.inc";
alias(i,j);

* Indice de precisão
Set l /1*100/;

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


display r, l;

$include "C:\Artigos\p-LD\numerical experiments\v2\problema_convetido\teste.inc";

Positive Variable x(i);
Positive Variable w(i,j);
Binary Variable z(j,l);
Positive Variable hat_x(i,j,l);
Positive Variable delta_x(j);
Positive Variable v(i,j);
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
MC_1
MC_2
MC_3
MC_4
;

Equations
R_FO_Or
R_FO_rel
Main_Or
;

R_FO_rel(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) =E= FO;

R_FO_Or(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) =E= FO;

Main_Or(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*x(i)*x(j)) + sum(i,b(r,i)*x(i)) + c(r) =L= 0;

parameter mult(r);
mult(r) = 0;

R_FO(r)$(re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + sum(r1$(not re_FO(r1)), mult(r1)*( sum((i,j)$BL(r1,i,j),a(r1,i,j)*w(i,j)) + sum(i,b(r1,i)*x(i)) + c(r1))) =E= FO;

Main(r)$(not re_FO(r))..
         sum((i,j)$BL(r,i,j),a(r,i,j)*w(i,j)) + sum(i,b(r,i)*x(i)) + c(r) =L= 0;

Def_w(i,j)$sum(r,BL(r,i,j))..
         w(i,j) =E= (x.up(j)-x.lo(j))*sum(l$(ord(l)<=-p),2**(-ord(l))*hat_x(i,j,l)) + (x.up(j)-x.lo(j))*v(i,j) + x(i)*x.lo(j);

Def_x(j)$sum((r,i),BL(r,i,j))..
         x(j) =E= (x.up(j)-x.lo(j))*sum(l$(ord(l)<=-p),2**(-ord(l))*z(j,l)) + (x.up(j)-x.lo(j))*delta_x(j) + x.lo(j);

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

* Bounds das variáveis
delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
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
scalar UB;
scalar LB;
UB = 10**8;
LB = -UB;
scalar FO_atual;
*option MIP = GUROBI;
option MIP = CPLEX;
option NLP = CONOPT;
option optcr = 0.0005;
stop = 0;
p = 0;

scalar time_Start, time_current, time_limit;
time_limit = 3600.0;
time_Start = TimeElapsed;

$ontext
while(not stop,
         delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
         iter_O = iter_O + 1;
*         mult(r1) = 0;
         iter_i = 0;
         stop_2 = 0;
         while (not stop_2,
                 time_current = TimeElapsed - time_Start;
                 if(time_current < time_limit,
                         m.reslim = time_limit - time_current;
                         Ori.reslim =  time_limit - time_current;
                 else
                         m.reslim = 0;
                         Ori.reslim = 0;
                 );
                 iter_i = iter_i + 1;
                 iter_tot = iter_tot + 1;
                 solve m min FO us MIP;
                 if(p=0,
                         FO_atual = FO.l;
                 else
                         FO_atual = m.objest;
                 );
                 LB = max(LB,FO_atual);
                 display FO_atual;
                 subg(r1)$(not re_FO(r1)) = sum((i,j)$BL(r1,i,j),a(r1,i,j)*w.l(i,j)) + sum(i,b(r1,i)*x.l(i)) + c(r1);
                 mod = (sum(r1$(not re_FO(r1)), subg(r1)*subg(r1)))**(1/2);

*                 p_aux = p;
*                 p = 0;
*                 solve Rel min FO us MIP;

                 solve Ori min FO us NLP;
*                 p = p_aux;
                 UB = min(UB,FO.l);
                 display UB, LB, step, mult, p, subg;
                 step = 2*(UB - FO_atual)/(mod**2);

                 mult_old(r1)$(not re_FO(r1)) = mult(r1);
                 mult(r1)$(not re_FO(r1)) = max(0,mult(r1) + step*subg(r1));

                 if( iter_i >= 5 or (UB - LB < 10**(-3)) or sum(r1,abs(mult(r1)-mult_old(r1)))<=10**(-6) or TimeElapsed > time_limit,
                         stop_2 = 1;
                 );
                 display iter_o, iter_i, stop, stop_2, p, UB, LB;
         );

         telapsed = TimeElapsed;
         if( iter_o >= 100 or (UB - LB < 10**(-3)) or TimeElapsed > time_limit,
                 stop = 1;
         else
                 p = p - 1;
         );
);
*display mult, subg, mod, LB, UB, p, iter_O, iter_i;
display UB, LB, p, iter_i, iter_tot, iter_O, iter_i, telapsed;
$offtext




while(not stop,
         delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**p;
         iter_O = iter_O + 1;

         time_current = TimeElapsed - time_Start;
                 if(time_current < time_limit,
                         Rel.reslim = time_limit - time_current;
                         Ori.reslim =  time_limit - time_current;
                 else
                         Rel.reslim = 0;
                         Ori.reslim = 0;
                 );

                 iter_tot = iter_tot + 1;
                 solve Rel min FO us MIP;
                 if(p=0,
                         FO_atual = FO.l;
                 else
                         FO_atual = Rel.objest;
                 );
                 LB = max(LB,FO_atual);
                 display FO_atual;


                 solve Ori min FO us NLP;
                 UB = min(UB,FO.l);
                 display UB, LB, p;


         telapsed = TimeElapsed;
         if( iter_o >= 100 or (UB - LB < 10**(-3)) or TimeElapsed > time_limit,
                 stop = 1;
         else
                 p = p - 1;
         );
);
*display mult, subg, mod, LB, UB, p, iter_O, iter_i;
display UB, LB, p, iter_i, iter_tot, iter_O, iter_i, telapsed;

file results /C:\Artigos\p-LD\numerical experiments\v2\problema_convetido\results.csv/;
put results;
         put telapsed /;
         put LB:10:4 /;
         put UB:10:4 /;
         put iter_i /;
         put iter_o /;
         put iter_tot /;
         put p;
putclose;

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
