scalar aux_x;
file results_parcial  /C:\Users\tiago.andrade\Desktop\RNMDT_Julia\WasterWater result\RNMDT.txt/;

put results_parcial;
         put 'RNMDT' /;
putclose;
results_parcial.ap = 0;

Parameter maior_valor /0/;
Set maior_j(j);
maior_j(j) := 0;
iter_O = 0;
while(not stop,
         delta_x.up(j)$sum((r,i),BL(r,i,j)) := 2**pp(j);
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
                 x_relaxado(i) = x.l(i);
                 w_relaxado(i,j) = w.l(i,j);
                 v_relaxado(i,j) = v.l(i,j);
                 delta_x_relaxado(i) = delta_x.l(i);
                 if(sum(j, pp(j)) = 0 and card(e) = 0,
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
                 maior_valor = smax(j, sum((r,i)$BL(r,i,j), abs(a(r,i,j)*(w_relaxado(i,j) - x_relaxado(i)*x_relaxado(j)))));
                 maior_j(j) = 0;
                 if((mod(iter_O+1,10) = 0),
                         pp(j)$sum((r,i),BL(r,i,j)) := pp(j) -1;
                 else
                         maior_j(j) = no;
                         while(sum(j$maior_j(j), 1) < 5 and sum(j$(not maior_j(j) and sum((r,i),BL(r,i,j))), 1) > 0,
                                 stop_2 = 0;
                                 maior_valor = smax(j$(not maior_j(j)), sum((r,i)$BL(r,i,j), abs(a(r,i,j)*(w_relaxado(i,j) - x_relaxado(i)*x_relaxado(j)))));
                                 loop(j$(not maior_j(j) and stop_2 = 0),
                                         if(sum((r,i)$BL(r,i,j), abs(a(r,i,j)*(w_relaxado(i,j) - x_relaxado(i)*x_relaxado(j)))) >= maior_valor*.999,        maior_j(j) = yes; stop_2 = 1;);
                                 );
                         );
*                        loop(j$(sum((r,i)$BL(r,i,j), abs(a(r,i,j)*(w.l(i,j) - x.l(i)*x.l(j)))) >= maior_valor*.5 and sum((r,i),BL(r,i,j))),
*                                maior_j(j) = 1;
*                        );
                        pp(j)$maior_j(j) :=  pp(j) - 1;
*min(pp(j),ceil(smin(i$sum(r,BL(r,i,j) and (abs(v_relaxado(i,j)-x_relaxado(i)*delta_x_relaxado(j))>0)), log2(4*abs(v_relaxado(i,j)-x_relaxado(i)*delta_x_relaxado(j))/(UB_x(i)-LB_x(i))))) - 1) ;
                        display maior_j, pp, maior_valor;
                 );
                 p = p - 1;
         );
);
*