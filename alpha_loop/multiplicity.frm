#-
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
* Dot product to symbols
id k1?.k2? = dotp(k1*k2);
argtoextrasymbol dotp;
id dotp(x?) = x;
.sort

* Infer multiplicity factor
L sF1 = firstterm_(F1);
L sF2 = firstterm_(F2);
L isoF = 
	#if termsin(sF1)==0
		#if termsin(sF2)==0
			0;
		#endif
	#else
		sF2*1/sF1;
	#endif
id 1/rat(x1?,x2?) = rat(x2,x1);
.sort
drop sF1, sF2;

* Apply multiplicity factor
L isoCHECK = isoF * F1 - F2;
.sort
* sort again to work around a tform bug having to do with normalization of the rat
.sort

#if termsin(isoCHECK) != 0
	print +s isoCHECK;
	.sort;
	exit "ISO CHECK FAIL: multiplicity factor not found";
#endif
print isoF;
.sort
.end;