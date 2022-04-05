#-
*Off statistics;
S D, ep;
Auto V k,p,q,fmb;
Auto S n,x,m,t,cMi,log,z,ALARM;
CF vx,rat,map,uvprop,g(s),vxs(s),tmps(s);
V energyselector;

Auto I s=4,mu=D;
CF counter,vec,vec1,gamma,NN;
CT opengammastring;

#define MAXPOLE "3"
#define SELECTEDEPSILONORDER "0"
#define PYSECDECCOMPARE "1"

Polyratfun rat;

#include ../../tensorreduce.frm
#include ../../integrateduv.frm

L F = %s;

.sort

#call IntegrateUV()

.sort

#call SubstituteMasters()

.sort

*Format float;
Format Float 32;

id pi^n? = (30246273033735921/9627687726852338)^n;

* for mUV=2
*id mUV2^n? = 2^n;
*id logmUVmu = 18947212500888319/13667524756850476;

* for mUV=1
id mUV2^n? = 1;
id logmUVmu = 0;

argument rat;
    id ep^%d = 0;
endargument;

Print +s;
.end
