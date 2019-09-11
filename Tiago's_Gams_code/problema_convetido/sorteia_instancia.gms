*base
Set i / 1* 20/;
Set e(*) / /;
Set r / 1* 30/;
Set sc /1 * 10/;



*sorteio
b(r,i) := uniform(-1,1);
c(r)$(not re_FO(r)) := uniform(-5,0);
iVc('1') := yes;
iVc('2') := yes;
iVc('3') := yes;
Q(r,i,j)$(i<j and not iVc(i)) := uniform(-1,1);
prob(sc) := 1/card(sc);
cStochastic(r,sc) := uniform(-2,1);

