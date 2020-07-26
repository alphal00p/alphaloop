#-
Off statistics;

*--#[ setup :

#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define H "25"
#define GHO "82"
#define GHOBAR "-82"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"
#define L "11,12,13"
#define LBAR "-11,-12,-13"
#define HELSUM "false"

Auto S mass;
CTable masses(-30:30);
CTable charges(-30:30);

#ifndef `OPTIMITERATIONS'
    #define OPTIMITERATIONS "100"
#endif

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = 0;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = 0;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;
#else
Fill masses(1) = mass_d;
Fill masses(2) = mass_u;
Fill masses(3) = mass_c;
Fill masses(4) = mass_s;
Fill masses(5) = mass_b;
Fill masses(-1) = mass_d;
Fill masses(-2) = mass_u;
Fill masses(-3) = mass_c;
Fill masses(-4) = mass_s;
Fill masses(-5) = mass_b;
Fill masses(11) = mass_e;
Fill masses(12) = mass_mu;
Fill masses(13) = mass_tau;
Fill masses(-11) = mass_e;
Fill masses(-12) = mass_mu;
Fill masses(-13) = mass_tau;
#endif

Fill masses(6) = mass_t;
Fill masses(-6) = mass_t;
Fill masses(25) = mass_h;

Fill charges(1) = -1/3; * d
Fill charges(2) = 2/3; * u
Fill charges(3) = -1/3; * s
Fill charges(4) = 2/3; * c
Fill charges(5) = -1/3; * b
Fill charges(6) = 2/3; * t
Fill charges(11) = -1; * e
Fill charges(-1) = 1/3; * d
Fill charges(-2) = -2/3; * u
Fill charges(-3) = 1/3; * s
Fill charges(-4) = -2/3; * c
Fill charges(-5) = 1/3; * b
Fill charges(-6) = -2/3; * t
Fill charges(-11) = 1; * e

S D, ep;
V p1,...,p40,k1,...,k40,c1,...,c40; * force this internal ordering in FORM
Auto V p,k,c, eps,ceps, sV,sVbar,sU,sUbar;
Auto S lm,ext;
Auto I mu=D,s=D;
Symbol ge, gs, gy, ghhh, type, in, out, virtual;
Auto S y,x, idx, t, n;

Set dirac: s1,...,s80;
Set diracdummy: sd1,...,sd80;
Set lorentz: mu1,...,mu80;
Set lorentzdummy: mud1,...,mud80;

CF gamma, spinor ,vector,g(s),delta(s),T, counter,color, prop;
CF f, vx, vec, vec1;
CF subs, configurations, conf, cmb, der, energy, spatial(s), spatialComp;
CF subgraph, uvconf, uvconf1, uvprop, uv;
S integratedctflag, mUV, logmUV, mi1L1, alarmMi1L1;
CF integratedct, rat, num, den;
CF gammaCONF;

Set ts: t0,...,t20;
CT penergy;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli80;
Set colA: cOlj1,...,cOlj80;
Set colAdum: cOljj1,...,cOljj40;

CF hermconjugate, pol, cpol, uSpinor, ubarSpinor, vSpinor, vbarSpinor, sp(s), d, sunA, sunF, diracInd, lorentzInd,sunTF,sunTA,vec;
S ii,m,y;
Auto I i1, i2, j1,j2;


Polyratfun rat;

*--#] setup :

* Load the diagrams
#include- input_`SGID'.h
#$epsCount = 0;
#$cepsCount = 0;
#$vCount = 0;
#$vbarCount =0;
#$uCount = 0;
#$ubarCount = 0;
*--#[ feynman-rules :
.sort
*do pre-counting in case the amplitude was proccessed already (e.g. by colorbasis. frm)
#do i = 1,50
        if ( (occurs(eps`i')>0) ) ; 
            #$epsCount=$epsCount+1;
        endif;
        if ( (occurs(ceps`i')>0));
            #$cepsCount=$cepsCount+1;
        endif;
        if ( (occurs(sU`i')>0) );
          #$uCount=$uCount+1;
        endif;
        if ( (occurs(sV`i')>0) );
          #$vCount=$vCount+1;
        endif;
        if ( (occurs(sUbar`i')>0) );  
            #$ubarCount=$ubarCount+1;
        endif;
        if ( (occurs(sVbar`i')>0) );
          #$vbarCount=$vbarCount+1;
        endif;        
        .sort        
#enddo
.sort

*id hermconjugate(x?) = x;
************************************************
* Implement the inteferences of scalar integrals
************************************************

*id hermconjugate(x?) = x;
************************************************
* Implement the inteferences of scalar integrals
************************************************

* now on the left-diagram
* extract polarizations and spinors
id d(i1?,i2?) = d_(i1,i2);
* construct gamma string
    repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
    repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p,?b,s2);
    Multiply counter(1);
* vectors
    repeat id gamma(?a,pol(?b),?c)*counter(i?) = pol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id gamma(?a,cpol(?b),?c)*counter(i?) = cpol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id sp(pol(?a),p?)*counter(i?) = pol(?a,lorentzdummy[i])*vec(p,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),p?)*counter(i?) = cpol(?a,lorentzdummy[i])*vec(p,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),pol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*pol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),cpol(?b))*counter(i?) = cpol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),cpol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    id counter(n?) = 1;

* spinors
    Multiply counter(1);
    repeat id gamma(ubarSpinor(?a),?b)*counter(i?) = ubarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(vbarSpinor(?a),?b)*counter(i?) = vbarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,uSpinor(?a))*counter(i?) = uSpinor(?a,diracdummy[i])*gamma(?b,diracdummy[i])*counter(i+1);
    repeat id gamma(?b,vSpinor(?a))*counter(i?) = vSpinor(?a,diracdummy[i])*gamma(?b,diracdummy[i])*counter(i+1);
    id counter(n?) =1;
* replace by proper indices from the corresponding set
argument;
    id sunA(x?) = colA[x];
    id sunF(x?) = colF[x];
    id diracInd(x?) = dirac[x];
    id lorentzInd(x?) = lorentz[x];
endargument;
Print +ss;
.sort
id d(i1?,i2?) = d_(i1,i2);
* count polarizations and spinors on left diagram
#do i = 1,50
        if ( (count(pol,1)>0) ) ; 
            #$epsCount=$epsCount+1;
        endif;
        repeat id pol($epsCount,?a,i1?) = vec(eps`i',i1);
        if ( count(pol,1)==0 ) redefine i "50";
        .sort        
#enddo

#do i = 1,50
        if ( (count(cpol,1)>0));
            #$cepsCount=$cepsCount+1;
        endif;
        id cpol($cepsCount,?a,i1?) = vec(ceps`i',i1);
        if ( count(cpol,1)==0 ) redefine i "50";
        .sort        
#enddo
Print +ss;
.sort
#do i = 1,50
        if ( (count(uSpinor,1)>0) );
          #$uCount=$uCount+1;
        endif;
        id uSpinor($uCount,?a,i1?) = spinor(sU`i',i1);
        if ( count(uSpinor,1)==0 ) redefine i "50";
        .sort        
#enddo

#do i = 1,50
        if ( (count(vSpinor,1)>0) );
          #$vCount=$vCount+1;
        endif;
        id uSpinor($vCount,?a,i1?) = spinor(sV`i',i1);
        if ( count(vSpinor,1)==0 ) redefine i "50";
        .sort        
#enddo

#do i = 1,50
        if ( (count(ubarSpinor,1)>0) );  
            #$ubarCount=$ubarCount+1;
        endif;
        id ubarSpinor($ubarCount,?a,i1?) = spinor(sUbar`i',i1);
        if ( count(ubarSpinor,1)==0 ) redefine i "50";
        .sort        
#enddo

#do i = 1,50
        if ( (count(vbarSpinor,1)>0) );
          #$vbarCount=$vbarCount+1;
        endif;
        id vbarSpinor($vbarCount,?a,i1?) = spinor(sVbar`i',i1);
        if ( count(vbarSpinor,1)==0 ) redefine i "50";
        .sort        
#enddo
**Print +s; 
.sort


* play the same game with the right-diagram
argument hermconjugate;
* construct gamma string
    repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
    repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p,?b,s2);
    Multiply counter(10);
* vectors
    repeat id gamma(?a,pol(?b),?c)*counter(i?) = pol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id gamma(?a,cpol(?b),?c)*counter(i?) = cpol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id sp(pol(?a),p?vector_)*counter(i?) = pol(?a,lorentzdummy[i])*vec(p,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),p?vector_)*counter(i?) = cpol(?a,lorentzdummy[i])*vec(p,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),pol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*pol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),cpol(?b))*counter(i?) = cpol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),cpol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    id counter(n?) = 1;
* spinors
    Multiply counter(10);
    repeat id gamma(ubarSpinor(?a),?b)*counter(i?) = ubarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(vbarSpinor(?a),?b)*counter(i?) = vbarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,uSpinor(?a))*counter(i?) = uSpinor(?a,diracdummy[i])*gamma(?b,diracdummy[i])*counter(i+1);
    repeat id gamma(?b,vSpinor(?a))*counter(i?) = vSpinor(?a,diracdummy[i])*gamma(?b,diracdummy[i])*counter(i+1);
    id counter(n?) =1;
endargument;
.sort 

Splitarg,hermconjugate;
repeat id hermconjugate(x?,y?,?a)= hermconjugate(x)+hermconjugate(y,?a);
FactArg,hermconjugate;
repeat id hermconjugate(x?,y?,?a)= hermconjugate(x)*hermconjugate(y,?a);
.sort
* manipulate indices for the right diagrams
argument hermconjugate;
    argument;
        id sunA(x?) = sunA(x+40);
        id sunF(x?) = sunF(x+40);
        id diracInd(x?) = diracInd(x+40);
        id lorentzInd(x?) = lorentzInd(x+40);
        id sunA(x?) = colA[x];
        id sunF(x?) = colF[x];
        id diracInd(x?) = dirac[x];
        id lorentzInd(x?) = lorentz[x];
    endargument;
    id d(i1?,i2?) = d_(i1,i2);
endargument;
.sort
id  hermconjugate(ii) = -ii;
id hermconjugate(vSpinor(?a)) =vbarSpinor(?a);
id hermconjugate(vbarSpinor(?a)) =vSpinor(?a);
id hermconjugate(uSpinor(?a)) =ubarSpinor(?a);
id hermconjugate(ubarSpinor(?a)) =uSpinor(?a);
id hermconjugate(pol(?a)) = cpol(?a);
id hermconjugate(cpol(?a)) = pol(?a);

* assume no gamma5
argument hermconjugate;
    Transform, gamma, reverse(1,last);
endargument;
id hermconjugate(x?)=x;
*Print;
id sunTF(?a) =T(?a);
id sunTA(?a) = cOlf(?a);
.sort

* update the spinor/polarization counts for spinors from the right-hand diagram.
#$offsetVB = $vbarCount;
#$offsetV = $vCount;
#$offsetU = $uCount;
#$offsetUB = $ubarCount;
#$offsetEP = $epsCount;
#$offsetCEP = $cepsCount;
.sort
#do i=1,50
    if(match(vbarSpinor(`i'-$offsetVB,?a))>=1);
        id vbarSpinor(`i'-$offsetVB,?a,i1?) =  spinor(sVbar`i',i1);
        #$vbarCount = $vbarCount+1;
    endif;
    if(match(vSpinor(`i'-$offsetV,?a))>=1);
        id vSpinor(`i'-$offsetV,?a,i1?) =  spinor(sV`i',i1);
        #$vCount = $vCount+1;
    endif;
    if(match(uSpinor(`i'-$offsetU,?a))>=1);
        id uSpinor(`i'-$offsetU,?a,i1?) =  spinor(sU`i',i1);
        #$uCount = $uCount+1;
    endif;
    if(match(ubarSpinor(`i'-$offsetUB,?a))>=1);
        id ubarSpinor(`i'-$offsetUB,?a,i1?) =  spinor(sUbar`i',i1);
        #$ubarCount = $ubarCount+1;
    endif;
    if(match(pol(`i'-$offsetEP,?a))>=1);
        id pol(`i'-$offsetEP,?a,i1?) =  vec(eps`i',i1);
        #$epsCount = $epsCount+1;
    endif;    
    if(match(cpol(`i'-$offsetCEP,?a))>=1);
        id cpol(`i'-$offsetCEP,?a,i1?) =  vec(ceps`i',i1);
        #$cepsCount = $cepsCount+1;
    endif;
    if((occurs(vbarSpinor,vSpinor,ubarSpinor,uSpinor,pol,cpol)==0));
        redefine i "50";
    endif;
    .sort
#enddo
id ii=i_;
Print +ss;

if (count(pol,1,cpol,1,uSpinor,1,vSpinor,1,ubarSpinor,1,vbarSpinor,1));
    Print "Unsubstituted polarization: %t";
    exit "Critical error";
endif;

* construct gamma string
id vec(p?,mu?) = p(mu);
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p1,?b,s2);
repeat id gamma(?a,p1?+p2?,?b) = gamma(?a,p1,?b)+gamma(?a,p2,?b);
repeat id gamma(?a,p1?-p2?,?b) = gamma(?a,p1,?b)-gamma(?a,p2,?b);
repeat id sp(p1?,p2?) = p1.p2;
.sort

Multiply counter(1);
repeat id gamma(i1?,?a,p?,?b,i2?)*counter(i?) = vec(p,lorentzdummy[i])*gamma(i1,?a,lorentzdummy[i],?b,i2)*counter(i+1);
id counter(n?) = 1;
*Do gamma traces
#do i=1,10
    id once gamma(mu?,?a,mu?) = g_(`i',?a);
#enddo
.sort

#do i=1,10
    tracen `i';
    .sort:trace-`i';
#enddo

#call ChisholmIdentities

id vec(p?,mu?) = p(mu);

.sort

.sort
CF gammaAll;
S iter, intSym ;

id gamma(?a) = gammaAll(gamma(?a));
B+ gammaAll;
.sort
keep brackets;
argument gammaAll;
    #do i=1,1
        label retry2;
        id once ifmatch->retry2 gamma(?a,p?,p?,?b) = gamma(?a,?b)*p.p;
        id once ifmatch->retry2 gamma(?a,p?,?b,mu?,p?,?c) = 2*p(mu)*gamma(?a,p,?b,?c)-gamma(?a,p,?b,p,mu,?c);
        id once ifmatch->retry2 gamma(?a,p?,?b,k?,p?,?c) = 2*p.k*gamma(?a,p,?b,?c)-gamma(?a,p,?b,p,k,?c);
    #enddo
endargument;
id gammaAll(x?) = x;
.sort
Table explSpinor(1:4,p?);
Fill explSpinor(1) = penergy(p);
Fill explSpinor(2) = spatialComp(p,1);
Fill explSpinor(3) = spatialComp(p,2);
Fill explSpinor(4) = spatialComp(p,3);
Print +ss;
.sort
repeat id once gamma(?a,mud1?,?b)*gamma(?c,mud1?,?d) = sum_(y,0,3,gamma(?a,y,?b)*gamma(?c,y,?d));
repeat id once spinor(p?,sd?diracdummy)*gamma(sd?diracdummy,?a) = sum_(y,1,4,spinor(p,y)*gamma(y,?a));
repeat id once spinor(p?,sd?diracdummy)*gamma(?a,sd?diracdummy) = sum_(y,1,4,spinor(p,y)*gamma(?a,y));
repeat id once spinor(p?,s?dirac)*gamma(s?dirac,?a) = sum_(y,1,4,spinor(p,y)*gamma(y,?a));
repeat id once spinor(p?,s?dirac)*gamma(?a,s?dirac) = sum_(y,1,4,spinor(p,y)*gamma(?a,y));

id spinor(p?,intSym?int_) = explSpinor(intSym,p);
.sort

id gamma(?a) = gammaAll(gamma(?a));
B+ gammaAll;
.sort
nTable slash(1:4,1:4,p?);
nTable gam(1:4,1:4,0:3);

#include- definition_gamma_explicit.h
*********************** EXPAND GAMMA FUNCTIONS ************************************************************
keep brackets;

argument gammaAll;
    #do i=1,1
        label retry4;        
        id once ifmatch->retry4  gamma(x?int_, p?vector_,  y?int_) = slash(x,y,p);
        id once ifmatch->retry4  gamma(x?int_, intSym?int_,  y?int_) = gamma(x,y,intSym);    
        id once ifmatch->retry4  gamma(x?int_,?a, p?vector_,  y?int_) = gamma(x,?a,1)*slash(1,y,p)+gamma(x,?a,2)*slash(2,y,p)+gamma(x,?a,3)*slash(3,y,p)+gamma(x,?a,4)*slash(4,y,p);
*        sum_(iter,1,4,gamma(x,?a,iter)*slash(iter,y,p));
        id once ifmatch->retry4  gamma(x?int_,?a, intSym?int_, y?int_) = sum_(iter,1,4,gamma(x,?a,iter)*gam(iter,y,intSym));   
*       if ( count(gamma,1) ) redefine i "0";
    #enddo
endargument;
.sort
#do i=1,1
    id once gammaAll(x?) =x;
    if ( count(gammaAll,1) > 0 ) redefine i "0";
    .sort
#enddo

.sort

*************************************************
* Process different configurations (bubbles, etc)
*************************************************

* process all the configurations
id configurations(x?) = x;
* multiply the numerator contribution of derivatives
id conf(?a,p?) = conf(?a) * penergy(p);
id conf(?a,x?) = conf(?a) * x; * note the type difference
* Factor out the mass
.sort
Polyratfun;
id rat(x1?,x2?) = num(x1)*den(x2);
.sort
FactArg den, num;
ChainOut,num;
ChainOut,den;
Multiply replace_(D, 4 - 2 * ep);

id num(x1?)=x1;
id den(x1?number_)=1/x1;
id den(mUV) = 1/mUV;
id den(x1?) = rat(1,x1);
id ep^n1?  = rat(ep^n1,1);
.sort:mass-factorization;
PolyRatFun rat(expand,ep,10);
.sort:ep-expansion;
PolyRatFun;
id rat(x1?) = x1;
if (count(ep, 1)) Discard; * keep only the ep^0 piece

.sort

* convert the polynomial to the cut momentum basis
id conf(x?,cmb(?a),?b) = conf(x,?b)*replace_(?a);
.sort:cmb;

* now extract the energy components of the LTD loop variables
id k1?.k2? = g(k1, k2);
repeat id conf(x?,?a,k1?,?b)*penergy(k1?) = conf(x,?a,k1,?b)*energy(k1); * detect loop energies from the bubble derivative
repeat id conf(x?,?a,k1?,?b)*g(k1?,k1?) = conf(x,?a,k1,?b)*(energy(k1)*energy(k1)-spatial(k1,k1));
repeat id conf(x?,?a,k1?,?b,k2?,?c)*g(k1?,k2?) = conf(x,?a,k1,?b,k2,?c)*(energy(k1)*energy(k2)-spatial(k1,k2));
repeat id conf(x?,?a,k?,?b)*g(k?,p?) = conf(x,?a,k,?b)*(energy(k)*penergy(p)-spatial(k,p));
id g(p1?,p2?) = p1.p2;

repeat id energy(?a)*energy(?b) = energy(?a,?b);
symmetrize energy;
id energy(?a) = energy(f(?a));
if (count(energy,1) == 0) Multiply energy(f(c0)); * signal with c0 that we are dealing with the constant term
.sort:energy-splitoff;

*********************************************
* Construction of optimized numerator C code
*********************************************
* Convert the dot products and energies to a symbol
#$MAXK = `NFINALMOMENTA';
#$MAXP = `NINITIALMOMENTA';
* for polarized cross-sections
#$MAXEPS = $epsCount;
#$MAXCEPS = $cepsCount;
#$MAXV = $vCount;
#$MAXVBAR = $vbarCount;
#$MAXU = $uCount;
#$MAXUBAR = $ubarCount;

#$OFFSET = 0;
#do i=1,`$MAXP'
    id penergy(p`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(p`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo

#enddo

#do i=1,`$MAXK'
    id penergy(c`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(c`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=1,`$MAXP'
        id c`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        id spatial(p`j', c`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=`i',`$MAXK'
        id c`i'.c`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
        id spatial(c`i', c`j') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXEPS'
    id penergy(eps`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(eps`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo


    #do j=`i', `$MAXEPS'
        id eps`i'.eps`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
    #enddo

    #do j=1, `$MAXK'
        id eps`i'.c`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
        id spatial(c`j', eps`i') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
    
    #do j=1, `$MAXP'
        id eps`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
    #enddo
    
    #do j=1, `$MAXCEPS'
        id eps`i'.ceps`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
    #enddo

#enddo

#do i=1,`$MAXCEPS'
    id penergy(ceps`i') =   lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(ceps`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo


    #do j=`i', `$MAXCEPS'
        id ceps`i'.ceps`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
    #enddo

    #do j=1, `$MAXK'
        id ceps`i'.c`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
        id spatial(c`j', ceps`j') = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
    
    #do j=1, `$MAXP'
        id eps`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET+1;
    #enddo    
#enddo
#do i=1,`$MAXV'
    id penergy(sV`i') =   lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(sV`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo
#do i=1,`$MAXVBAR'
    id penergy(sVbar`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(sVbar`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo
#do i=1,`$MAXU'
    id penergy(sU`i')  = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(sU`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo
#do i=1,`$MAXUBAR'
    id penergy(sUbar`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j =1,3
        id spatialComp(sUbar`i',`j') =  lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo


* split off every energy configuration into a new expression
id conf(?a) = conf(conf(?a));
argtoextrasymbol tonumber,conf,1;
#redefine oldextrasymbols "`extrasymbols_'"
B+ conf;
.sort:conf-1;
Hide F;
#redefine energysymbolstart "`extrasymbols_'"
#do ext={`oldextrasymbols'+1}, `energysymbolstart'
    delete extrasymbols>`energysymbolstart'; * clear all extra symbols from the last configuration
    #$conf = extrasymbol_(`ext');
    #write<out_`SGID'.proto_c> "#CONF\n%$", $conf;
    L FF`ext' = F[conf(`ext')];
    .sort
    argtoextrasymbol energy,1;
    id energy(x?) = x;
    .sort
    #define energysymbolend "`extrasymbols_'"
    B+ <Z{`energysymbolstart' + 1}_>,...,<Z`energysymbolend'_>;
    .sort:conf-2;

* Optimize the output
    Format C;
    Format O1,stats=off,saIter=`OPTIMITERATIONS';
    #Optimize FF`ext'
    #write<out_`SGID'.proto_c> "%O"
    B+ <Z{`energysymbolstart' + 1}_>,...,<Z`energysymbolend'_>;
    .sort:optim-`ext'-1;

    #do symb={`energysymbolstart' + 1}, `energysymbolend'
        #$energyexpr = FF`ext'[Z`symb'_];
        #$energyconf = extrasymbol_(`symb');
        #write<out_`SGID'.proto_c> "#NEWMONOMIAL\n%$\nreturn %$;",$energyconf, $energyexpr
    #enddo

    #clearoptimize;
    .sort:optim-`ext'-2;
    Format O0;
    Format normal;
    Drop FF`ext';
#enddo
.end
