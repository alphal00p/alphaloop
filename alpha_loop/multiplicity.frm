#-
On statistics;
CF dotp, myv;

#define MULPHASE "1"

#include numerator.frm # setup 

* Load the diagrams
#include iso_check_`SGID'.frm

.sort
Drop F1,...,F`NUMD';
L F = F1 + x * (F2+...+F`NUMD');

* Do algebra
#include numerator.frm # feynman-rules

#call FeynmanRules()

*patch form factors so it doesn't break the check
#ifdef `FORMFACTORS'
#call MultiplictyPatch()
#endif

B+ x;
.sort
Drop F;

L F1 = F[1];
L F2 = F[x];
.sort

********************************
* Check multiplicity 
********************************
* Dot product to symbols
#redefine oldextrasymbols "`extrasymbols_'"
id k1?.k2? = dotp(k1.k2);
argtoextrasymbol dotp;
id dotp(x?) = x;
.sort
Hide F1, F2;

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
print +s isoF;
.sort
drop sF1, sF2;

* Apply multiplicity factor
L isoCHECK = isoF * F1 - F2;
.sort
* sort again to work around a tform bug having to do with normalization of the rat
.sort

#if termsin(isoCHECK) != 0
	Multiply replace_(<Z{`oldextrasymbols'+1}_,extrasymbol_({`oldextrasymbols'+1})>,...,<Z`extrasymbols_'_,extrasymbol_(`extrasymbols_')>);
	.sort
* Test if the expression is really non zero by expanding the scalar products	
	id k1?.k2? = dotp(myv(k1,0))*dotp(myv(k2,0))
		   - dotp(myv(k1,1))*dotp(myv(k2,1))
		   - dotp(myv(k1,2))*dotp(myv(k2,2))
		   - dotp(myv(k1,3))*dotp(myv(k2,3));
	.sort
	
	#if termsin(isoCHECK) != 0
		print +s isoCHECK, isoF;
		.sort;
		exit "ISO CHECK FAIL: multiplicity factor not found";
	#endif
#endif
print isoF;
.sort
.end;