#define epsGoal "0"
#-
CFunction  topo1, topo2;

Symbol mUV, d, eps, n1,...,n3, tp2P011, tp2P111, mi1L1 ,mi1L2, mi2L2, lMUV, x1, x2;
Vectors V p1,...,p40,k1,k2; 
CFunction rat, num,den,freeEps;
PolyRatFun rat;

*PolyRatFun expEps(expand,eps,0);
* coefficints of eps expansion: CMi; not computed coeff: alarmMI
AutoDeclare Symbol cMi, alarmMi;

*Local expr= topo2(3,2,3)+topo2(-5,4,3);
*(k1.k1*mUV^2+p1.p2*k1.k1)*topo1(3);
* full test checks IBPs
*#include fullTest.inc
* fullExpansionTest computes expansions 
#include fullExpansionTest.inc

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


* 1loop ibps
id topo1(n1?{<1}) = 0;
  repeat;
      id topo1(n1?{>1}) = topo1(-1 + n1)*rat((2 + d - 2*n1),2*mUV^2*(-1 + n1));
  endrepeat;
* topo1(2) as master
id topo1(1)=mi1L1*rat(2*mUV^2,d - 2);
.sort


* parametric reduction of 2loop
#call paraIBPSunrise()
.sort

* plug in eps expansion
multiply replace_(d,4-2*eps);
#include mi1L1epsEpxansion.inc;
#include mi1L2epsEpxansion.inc;
#include mi2L2epsEpxansion.inc;

.sort
*on nospacesinnumbers;
*format mathematica;
*#write <testExpA0.log>  "formRes=%E", expr1
*.sort
PolyRatFun;
id rat(x1?,x2?) = num(x1)*den(x2);
FactArg den, num;
ChainOut,num;
ChainOut,den;

id num(x1?)=x1;
id den(x1?number_)=1/x1;
* since we are only having one scale, the mass will factorize
id den(mUV) = 1/mUV;
id den(x1?)^n1? = rat(1,x1^n1);
id den(x1?) = rat(1,x1); 
* treat the numerator
id eps^n1?  = rat(eps^n1,1);
id eps = rat(eps,1);
repeat;
    id rat(x1?,x2?)*rat(n1?,n2?)= rat(x1*n1,x2*n2);
endrepeat;
.sort
PolyRatFun rat;
Symbol eps;
.sort
*on nospacesinnumbers;
*format mathematica;
*#write <testExpA.log>  "formRes=%E", expr1
*.sort;
PolyRatFun rat(expand,eps,10);
Symbol eps;
.sort
PolyRatFun;
id rat(x1?)=x1;
id eps^n1?{>'epsGoal'}=0;
*print;
.sort

on nospacesinnumbers;
format mathematica;
*#write <testExpB.log>  "formRes=%E", expr1
* export for running fullTestEps
*format mathematica;
#Do i=1,889;
    #write <expansionTest.log>  "res'i'=%e", int'i'
#enddo;  

.end
