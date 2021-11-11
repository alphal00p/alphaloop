#-
On statistics;
CF dotp, myv;

#include numerator.frm # setup

* Load the diagrams
#include- input_`SGID'.h

.sort

* Do algebra
#include numerator.frm # feynman-rules

#call FeynmanRules()

B+ x;
.sort

#write<out_raw_numerator_`SGID'.txt> "%E" F
.sort
.end;
