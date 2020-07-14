#-
#:maxtermsize 1M

Off statistics;
CF dotp;

#include numerator.frm # setup

* Load the diagrams
#include iso_check_`SGID'_`ID0'_`IDn'.frm
* Do algebra
#include numerator.frm # feynman-rules

********************************
* Check multiplicity 
********************************
Multiply replace_(D ,4);
Multiply replace_(gs ,rat(gs,1));
Multiply replace_(mass_t ,rat(mass_t,1));
id k1?.k2? = dotp(k1*k2);
argtoextrasymbol dotp;
id dotp(x?) = rat(x,1);
.sort
L isoF = 
	#if termsin(F1)==0
		#if termsin(F2)==0
			0;
		#endif
	#else
		F2*1/F1;
	#endif
id 1/rat(x1?,x2?) = rat(x2,x1);
print isoF;
.sort
.end;