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
Auto V p,k,c;
Auto S lm,ext;
Auto I mu=D,s=D;
Symbol ge, gs, gy, ghhh, type, in, out, virtual;
Auto S x, idx, t, n;

Set dirac: s1,...,s80;
Set diracdummy: sd1,...,sd80;
Set lorentz: mu1,...,mu80;
Set lorentzdummy: mud1,...,mud80;

CF gamma, vector,g(s),delta(s),T, counter,color, prop;
CF f, vx, vec, vec1;
CF subs, configurations, conf, cmb, der, energy, spatial(s);
CF subgraph, uvconf, uvconf1, uvprop, uv;
S integratedctflag, mUV, logmUV, mi1L1, alarmMi1L1;
CF integratedct, rat, num, den;


Set ts: t0,...,t20;
CT penergy;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli80;
Set colA: cOlj1,...,cOlj80;
Set colAdum: cOljj1,...,cOljj40;

CF hermconjugate, pol, cpol, uSpinor, ubarSpinor, vSpinor, vbarSpinor, sp, d, sunA, sunF, diracInd, lorentzInd;
S ii,m,x,y;
Auto I i1, i2, j1,j2;

Polyratfun rat;

*--#] setup :

* Load the diagrams
#include- input_`SGID'.h

*--#[ feynman-rules :


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
    repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p1,?b,s2);
    Multiply counter(1);
* vectors
    repeat id gamma(?a,pol(?b),?c)*counter(i?) = pol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id gamma(?a,cpol(?b),?c)*counter(i?) = cpol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id sp(pol(?a),p?vector_)*counter(i?) = pol(?a,lorentzdummy[i])*p(lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),p?vector_)*counter(i?) = cpol(?a,lorentzdummy[i])*p(lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),pol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*pol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),cpol(?b))*counter(i?) = cpol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),cpol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    id counter(n?) = 1;
* spinors
    Multiply counter(1);
    repeat id gamma(ubarSpinor(?a),?b)*counter(i?) = ubarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(vbarSpinor(?a),?b)*counter(i?) = vbarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,uSpinor(?a))*counter(i?) = uSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,vSpinor(?a))*counter(i?) = vSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    id counter(n?) =1;
* replace by proper indices from the corresponding set
argument;
    id sunA(x?) = colA[x];
    id sunF(x?) = colF[x];
    id diracInd(x?) = dirac[x];
    id lorentzInd(x?) = lorentz[x];
endargument;
.sort
* play the same game with the right-diagram
argument hermconjugate;
* construct gamma string
    repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
    repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p1,?b,s2);
    Multiply counter(1);
* vectors
    repeat id gamma(?a,pol(?b),?c)*counter(i?) = pol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id gamma(?a,cpol(?b),?c)*counter(i?) = cpol(?b,lorentzdummy[i])*gamma(?a,lorentzdummy[i],?c)*counter(i+1);
    repeat id sp(pol(?a),p?vector_)*counter(i?) = pol(?a,lorentzdummy[i])*p(lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),p?vector_)*counter(i?) = cpol(?a,lorentzdummy[i])*p(lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),pol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*pol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(cpol(?a),cpol(?b))*counter(i?) = cpol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    repeat id sp(pol(?a),cpol(?b))*counter(i?) = pol(?a,lorentzdummy[i])*cpol(?b,lorentzdummy[i])*counter(i+1);
    id counter(n?) = 1;
* spinors
    Multiply counter(1);
    repeat id gamma(ubarSpinor(?a),?b)*counter(i?) = ubarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(vbarSpinor(?a),?b)*counter(i?) = vbarSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,uSpinor(?a))*counter(i?) = uSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    repeat id gamma(?b,vSpinor(?a))*counter(i?) = vSpinor(?a,diracdummy[i])*gamma(diracdummy[i],?b)*counter(i+1);
    id counter(n?) =1;
endargument;
.sort 

id d(i1?,i2?) = d_(i1,i2);
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
id ii=i_;
.sort
*polarization sums
#if HELSUM
    repeat id cpol(p1?,0,i1?)*pol(p1?,0,i2?) = -d_(i1,i2);
    repeat id uSpinor(p1?,m?,s1?)*ubarSpinor(p1?,m?,s2?) = gamma(s1,p,s2)+m*d_(s1,s2);
    repeat id vSpinor(p1?,m?,s1?)*vbarSpinor(p1?,m?,s2?) = gamma(s1,p,s2)-m*d_(s1,s2);repeat id cpol(p1?,0,j1?,i1?)*pol(p1?,0,j2?,i2?) = -d_(j1,j2)*d_(i1,i2);
#endif
* construct gamma string
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p1,?b,s2);

*if (count(pol,1,cpol,1,uSpinor,1,vSpinor,1,ubarSpinor,1,vbarSpinor,1));
*    Print "Unsubstituted polarization sum: %t";
*    exit "Critical error";
*endif;

repeat id gamma(?a,p1?+p2?,?b) = gamma(?a,p1,?b)+gamma(?a,p2,?b);
repeat id gamma(?a,p1?-p2?,?b) = gamma(?a,p1,?b)-gamma(?a,p2,?b);
repeat id sp(p1?,p2?) = p1.p2;
*Print +s;
.sort

id D^n? = rat(D^n, 1);

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, ?a, cOli1?) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color(1);
repeat id d_(cOlj1?,cOlj2?) * color(x?) = color(x*d_(cOlj1,cOlj2));
repeat id d_(cOli1?,cOli2?) * color(x?) = color(x*d_(cOli1,cOli2));
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
repeat id cOlf(cOlj1?,cOlj2?,cOlj3?)*color(x?) = color(x * cOlf(cOlj1,cOlj2,cOlj3));
.sort

B+ color;
.sort:color-prep;

Keep brackets;
* TODO: do we need to set the dimension?
* Only evaluate this part per unique color by bracketing
Argument color;
    #call color
    #call simpli
    id  cOlI2R = cOlcR*cOlNR/cOlNA;
    id  cOlNR/cOlNA*cOlcR = cOlI2R;
    id  cOld33(cOlpR1,cOlpR2) = [dabc^2/n]*cOlNR;
    id  cOlNR/cOlNA = nf/cf/2;
    id  cOlcR = cf;
    id  cOlcA = ca;
    id  cOlI2R = nf/2;
	id	cOld44(cOlpA1,cOlpA2) = [d4AA/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpR2) = [d4RR/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpA1) = [d4RA/n]*cOlNA;

* set the SU(3) values
    id [dabc^2/n] = 13/18;
    id cOlNR = 3;
    id cOlNA = 8;
    id cf = 4 / 3;
    id ca = 3;
    id nf = 1;
EndArgument;

.sort:color;
Splitarg,color;
repeat id color(x?,y?,?a)= color(x)+color(y,?a);
FactArg,color;
repeat id color(x?,y?,?a)= color(x)*color(y,?a);
id color(x?number_) = x;

.sort
* set the SU(3) values
id cOlNA =8;
id cOlNR = 3;

* construct gamma string
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);

Multiply counter(1);
repeat id gamma(mu?,?a,p?,?b,mu?)*counter(i?) = vec(p,lorentzdummy[i])*gamma(mu,?a,lorentzdummy[i],?b,mu)*counter(i+1);
id counter(n?) = 1;

#do i=1,10
    id once gamma(mu?,?a,mu?) = g_(`i',?a);
#enddo

id vec(p?,mu?) = p(mu);
.sort

#do i=1,10
    tracen `i';
    .sort:trace-`i';
#enddo

*if (count(gamma, 1));
*    Print "Unsubstituted gamma string: %t";
*    exit "Critical error";
*endif;

id D^n? = rat(D^n, 1);
.sort:gamma-traces;

* in case there is no color structure
Multiply color(1);
repeat id  color(x?)*color(y?) = color(x*y);
brackets+ color;
.sort
* split off every color configuration into a new expression
id color(?a) = color(color(?a));
argtoextrasymbol tonumber,color,1;
#redefine oldextrasymbols "`extrasymbols_'"
B+ color;
.sort:color-1;
Hide F;
#redefine colorsymbolstart "`extrasymbols_'"
#write<SG_`SGID'_color_decomp.txt> "{\n"
#do cc={`oldextrasymbols'+1}, `colorsymbolstart'
    delete extrasymbols>`colorsymbolstart'; * clear all extra symbols from the last configuration
    #$color = extrasymbol_(`cc');
    #write<SG_`SGID'_color_decomp.txt> "\"%$\":", $color;
    #$struc = F[color(`cc')];
    #write<SG_`SGID'_color_decomp.txt> "\"%$\",\n", $struc;
#enddo
    #write<SG_`SGID'_color_decomp.txt> "}"
.end