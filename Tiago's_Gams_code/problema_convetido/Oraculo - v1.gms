delta_x.up(j)$sum((r,i),BL(r,i,j)) = 2**pp(j);

 UB_iteracao := 0;
 LB_iteracao := 0;
 x_stoc_rel(i,sc) := 0;
 x_stoc_ori(i,sc) := 0;

* obtem solucao dual
loop(sc1,
         display iteracao_grande, sc1;
         scElemento(sc) := no;
         scElemento(sc1) := yes;
display scElemento;
         solve Rel_lagrange us MIP max FO;
         UB_iteracao  := UB_iteracao  + FO.l*prob(sc1);
         x_stoc_rel(i,sc1) := x.l(i);
 );
 UB = min(UB, UB_iteracao);
 subgradient(i,sc)$dominio_multip(i,sc) := x_stoc_rel(i, sc) - x_stoc_rel(i,sc+1);
 valor_fixado(i)$iVc(i) := 0.0 + .0*smin(sc, x_stoc_rel(i, sc)) + 1*sum(sc, prob(sc)*x_stoc_rel(i, sc)) + 0.00*smax(sc, x_stoc_rel(i, sc));
display valor_fixado;
*sum(sc, x_stoc_rel(i, sc))/card(sc);
* obtem solucao primal
local_viavel := 1;
loop(sc1,
         display iteracao_grande, sc1;
         scElemento(sc) := no;
         scElemento(sc1) := yes;
         x.l(i) := x_stoc_rel(i,sc1);
         solve m_fixa_primeiro_estagio us NLP max FO;
         LB_iteracao  := LB_iteracao  + FO.l*prob(sc1);
         x_stoc_ori(i,sc1) := x.l(i);
         direcao_valor_fixado(i,sc1)$iVc(i) := fixa_var_primeiro_estagio.m(i);
         if(m_fixa_primeiro_estagio.modelStat > 2,
                 local_viavel := 0;
         );
 );
if(local_viavel = 1,
         LB := max(LB, LB_iteracao);
);
display valor_fixado,direcao_valor_fixado;
*valor_fixado(i)$iVc(i) := max(0,
*                          valor_fixado(i)
*                         + 0.5*sum(sc, prob(sc)*direcao_valor_fixado(i,sc))
*/sqrt(sum(sc, power(prob(sc)*direcao_valor_fixado(i,sc),2)))
*                         );
display valor_fixado;

loop(busca_local,
         local_viavel_busca := 1;

         LB_iteracao_local := 0;
         local_viavel := 1;
         loop(sc1,
                  display iteracao_grande, sc1;
                  scElemento(sc) := no;
                  scElemento(sc1) := yes;
*                  x.l(i) := x_stoc_rel(i,sc1);
                  x.l(i) := x_stoc_rel(i,sc1) * uniform(0.9,1.1);
                  solve m_fixa_primeiro_estagio us NLP max FO;
                  LB_iteracao_local  := LB_iteracao_local  + FO.l*prob(sc1);
                  x_stoc_ori(i,sc1) := x.l(i);
                  direcao_valor_fixado(i,sc1)$iVc(i) := fixa_var_primeiro_estagio.m(i);
                  if(m_fixa_primeiro_estagio.modelStat > 2,
                         local_viavel_busca := 0;
                 );
          );
         if(local_viavel_busca = 1,
                 LB_iteracao := max(LB_iteracao, LB_iteracao_local);
                 LB := max(LB, LB_iteracao);
                 local_viavel := 1;
         );
         display LB_iteracao_local, LB_iteracao, valor_fixado, direcao_valor_fixado;
         valor_fixado(i)$iVc(i) :=  valor_fixado(i)* uniform(0.9,1.1);
         valor_fixado(i)$iVc(i) := min(UB_x(i),max(LB_x(i),valor_fixado(i)));
*         valor_fixado(i)$iVc(i) := max(0,
*                                 valor_fixado(i)
*                         + 0.5*sum(sc, prob(sc)*direcao_valor_fixado(i,sc))
*./sqrt(sum(sc, power(prob(sc)*direcao_valor_fixado(i,sc),2)))
*                         );
         display valor_fixado;

);
