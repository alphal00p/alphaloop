#-
Off statistics;

#include numerator.frm # setup

* Load the diagrams
#include iso_check_`SGID'_`ID0'_`IDn'.frm

L FFF = 1;

L F12p = F1+F2;
L F12m = F1-F2;

#include numerator.frm # feynman-rules

********************************
* Check multiplicity +/-1
********************************
PolyRatFun;
L isoF = 
	#if termsin(F12m) == 0
        #if termsin(F1) == 0
            10;
        #else
		    1;
        #endif
	#elseif termsin(F12p) == 0
		-1;
	#else 
		0;
	#endif
	

print isoF;
.sort
.end;
