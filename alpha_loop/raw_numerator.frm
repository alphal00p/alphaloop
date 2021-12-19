#-
On statistics;
CF dotp, myv;

#include numerator.frm # setup

* Load the diagrams
#include- input_`SGID'.h

.sort

* Do algebra
#include numerator.frm # feynman-rules

* WARNING THE FEYNMAN RULES NEED TO BE REVAMPED SO AS TO SPLIT IT INTO: Couplings, Colour and Lorentz
* THis is necessary in order to support four-gluon vertices since in that case only couplings factorize completely.

L F1 = F;
.sort
Hide F1;
#call FeynmanRulesGlobal()
id vx(?a) = 1;
id prop(?a) = 1;
* Print +S F;
.sort
#write<out_raw_prefactor_`SGID'.txt> "%E" F
Drop F;

Unhide F1;
#call FeynmanRulesMomentum()
.sort
PolyRatFun rat(expand,ep,50);
id rat(x1?,x2?) = x1 / x2;

* Print +S F1;
#write<out_raw_numerator_`SGID'.txt> "%E" F1
.sort

.end;