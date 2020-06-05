
#define NLOOPS "2"
* use True if LookTables for IBPs should be used


#-
CFunction den, topo1, topo2;
Symbols mUV, d, eps, n1,...,n3, M1, M2;
Vectors V p1,...,p40,k1,k2; 
CFunction rat; 
Polyratfun rat;

*Local expr= topo2(3,2,3)+topo2(-5,4,3);
*(k1.k1*mUV^2+p1.p2*k1.k1)*topo1(3);
#include fullTest.inc
* example of slow ibp reductions: ibpreduced result - original seed 
* Local ibps1=-((-36633600 + 57571584*d - 34603072*d^2 + 10196480*d^3 - 1603480*d^4 + 135716*d^5 - 5758*d^6 + 95*d^7)*topo2(1, 1, 0))/(33592320*mUV^16) - ((-23042880 + 20938896*d - 7321734*d^2 + 1271809*d^3 - 117180*d^4 + 5614*d^5 - 126*d^6 + d^7)*topo2(1, 1, 1))/(11022480*mUV^14) + topo2(1, 1, 8);
*Local expr=topo2(2,2,2);



* replace partial fraction scalar-products 
repeat;
    id k1.k1*topo1(n1?) = topo1(n1-1)+mUV^2*topo1(n1);
endrepeat;


* topo2 =[k1^2, (k1-k2)^2,k2^2]
* k1*k2=-(prop2-k1.k1 -k2.k2 +mUV^2 )/2
repeat;
    id k1.k2*topo2(n1?,n2?,n3?) =-1/2*( 
        topo2(n1, n2-1,n3)
        +(- k1.k1 - k2.k2 + mUV)*topo2(n1,n2,n3) 
        );

    id k1.k1*topo2(n1?,?a) = topo2(n1-1,?a)+mUV^2*topo1(n1,?a);
    id k2.k2*topo2(?a,n3?) = topo2(?a,n3-1)+mUV^2*topo1(?a,n3);
endrepeat;



#if 'NLOOPS'==1
    repeat;
        id topo1(n1?{>1}) = ((2 + d - 2*n1)*topo1(-1 + n1))/(2*mUV^2*(-1 + n1));
    endrepeat;
    id topo1(n1?{<1}) = 0;
    .sort
#endif


#if 'NLOOPS'==2
* parametric reduction
    #call paraIBPSunrise()
    .sort
#endif


.sort
*print;

* export for running fullTest
format mathematica;
#Do i=1,889;
    #write <ibpTest.log>  "res['i']-> %e", ibps'i'
#enddo;  

.end
