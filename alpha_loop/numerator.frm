#-
* PDGs:
*d,u,c,s,b,t -> 1 to 6
*d~,u~,c~,s~,b~,t~ -> -1 to -6
*g, photon -> 21, 22
*e+ e- > -11, 11
*mu+, mu-, ta+ ,ta- > -12, 12, -13, 13
*z w+ w- z -> 23, 24, -24
*H 25

#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"

Auto S mass;
CTable masses(-30:30);
Fill masses(1)= mass_d;
Fill masses(2)= mass_u;
Fill masses(3)= mass_c;
Fill masses(4)= mass_s;
Fill masses(5)= mass_b;
Fill masses(6)= mass_t;
Fill masses(11) = mass_e;
Fill masses(12) = mass_mu;
Fill masses(13) = mass_tau;
Fill masses(25) = mass_h;
Fill masses(-1) = mass_d;
Fill masses(-2) = mass_u;
Fill masses(-3) = mass_c;
Fill masses(-4) = mass_s;
Fill masses(-5) = mass_b;
Fill masses(-6) = mass_t;
Fill masses(-11) = mass_e;
Fill masses(-12) = mass_mu;
Fill masses(-13) = mass_tau;

Auto V p,k;
Auto I mu=4,s=4;
Symbol mT,ge,ee,gs,ii,yt,n,m;
Auto S x, idx;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;

CF gamma, vector,g(s),delta(s),T, counter,color, prop;
CF u,ubar,v,vbar,f,vx;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  type, in, out, virtual, L, A;


#include- diacolor.h

* Load the diaggrams
#include- input.h

************************************************
* Substitute the Feynman rules for the numerator
************************************************

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);

* incoming edges
id prop(`EM', in, p?, idx?) = u(p, masses(`EM'), dirac[idx]);
id prop(x?{`Q'}, in, p?, idx?) = u(p, masses(x), dirac[idx]);
id prop(`EP', in, p?, idx?) = vbar(p, masses(`EP'), dirac[idx]);
id prop(x?{`QBAR'}, in, p?, idx?) = vbar(p, masses(x), dirac[idx]);

* outgoing edges
id prop(`EM', out, p?, idx?) = ubar(p, masses(`EM'), dirac[idx]);
id prop(`EP', out, p?, idx?) = v(p, masses(`EP'), dirac[idx]);
id prop(x?{`Q'}, out, p?, idx?) = ubar(p, masses(x), dirac[idx]);
id prop(x?{`Q'}, out, p?, idx?) = v(p, masses(x), dirac[idx]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]) * d_(cOljjj[idx1], cOljjj[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`FERM'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2])) * d_(cOliii[idx1], cOliii[idx2]);

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gs * i_ * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * T(cOliii[idx1], cOljjj[idx2], cOliii[idx3]);
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = 2/3 * ee * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * d_(cOliii[idx1], cOliii[idx3]);
id vx(`EP', `PHO', `EM', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ee * i_ * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);

if (count(vx, 1, prop, 1));
    Print "Unsubstituted propagator or vertex: %t";
    exit "Critical error";
endif;

Print +s;
.end
.sort:feynman-rules-vertices;

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, cOli1?, ?a) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color;
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
Print +s;
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
EndArgument;
Print +s;
.sort:color;

* Do the spin sum
repeat id u(p?,m?,s1?)*ubar(p?,m?,s2?) = gamma(s1, p, s2) + m*gamma(s1, s2);
repeat id v(p?,m?,s1?)*vbar(p?,m?,s2?) = gamma(s1, p, s2) - m*gamma(s1, s2);

* construct gamma string
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);

#do i=1,10
    Print "%t";
    id once gamma(mu?,?a,mu?) = g_(`i',?a);
    trace4 `i';
#enddo

if (count(gamma, 1));
    Print "Unsubstituted gamma string: %t";
    exit "Critical error";
endif;

.sort:gamma-traces;

id color(x?) = x;

id pzero = 0; * Substitute the 0-momentum by 0

.sort:feynman-rules-final;

*********************************************
* Construction of optimized numerator C code
*********************************************

* TODO: split off the energy components of the LTD loop momenta
* and optimize a monomial in the energies
#do i =1,10
    #do j=1,10
        id k`i'.p`j' = k`i'p`j';
        id k`i'.k`j' = k`i'`j';
    #enddo
#enddo
.sort
Format O4,stats=on,saIter=1000;
Format C;
#Optimize F
#write<out.proto_c> "%O\nreturn %e",F

Print +s;
.end