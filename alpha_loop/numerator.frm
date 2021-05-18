#-
On statistics;
On nospacesinnumbers;

*--#[ setup :

#define PSI "1337"
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
#define QMASSLESS "1,2,3,4,5"
#define QBARMASSLESS "-1,-2,-3,-4,-5"
#define QMASSIVE "6,"
#define QBARMASSIVE "-6,"
#define L "11,12,13"
#define LBAR "-11,-12,-13"

S vev, pi;

Auto S mass, mUV;
Auto S yukawa;
CTable masses(-99:1338);
CTable gyq(-30:30);
CTable logmasses(-30:30);
CTable charges(-30:30);

#ifndef `OPTIMLVL'
    #define OPTIMLVL "4"
#endif 

#ifndef `OPTIMITERATIONS'
    #define OPTIMITERATIONS "100"
#endif

Fill gyq(1) = yukawad; * d
Fill gyq(2) = yukawau; * u
Fill gyq(3) = yukawas; * s
Fill gyq(4) = yukawac; * c
Fill gyq(5) = yukawab; * b
Fill gyq(6) = yukawat; * t
Fill gyq(11) = 0; * e-
Fill gyq(12) = 0; * mu-
Fill gyq(13) = 0; * ta-
Fill gyq(-1) = yukawad; * d
Fill gyq(-2) = yukawau; * u
Fill gyq(-3) = yukawas; * s
Fill gyq(-4) = yukawac; * c
Fill gyq(-5) = yukawab; * b
Fill gyq(-6) = yukawat; * t
Fill gyq(-11) = 0; * e+
Fill gyq(-12) = 0; * mu+
Fill gyq(-13) = 0; * ta+

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = massb;
Fill masses(6) = masst;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = massb;
Fill masses(-6) = masst;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;
#else
Fill masses(1) = massd;
Fill masses(2) = massu;
Fill masses(3) = massc;
Fill masses(4) = masss;
Fill masses(5) = massb;
Fill masses(6) = masst;
Fill masses(-1) = massd;
Fill masses(-2) = massu;
Fill masses(-3) = massc;
Fill masses(-4) = masss;
Fill masses(-5) = massb;
Fill masses(-6) = masst;
Fill masses(11) = masse;
Fill masses(12) = massmu;
Fill masses(13) = masstau;
Fill masses(-11) = masse;
Fill masses(-12) = massmu;
Fill masses(-13) = masstau;
#endif

Fill masses(-82) = 0;
Fill masses(82) = 0;
Fill masses(21) = 0;
Fill masses(22) = 0;
Fill masses(25) = massh;
Fill masses(1337) = 0;

* note: this set needs to be complete for the UV expansion
Set allmasses: massu, massd, massc, masss, masst, massb, masse, massmu, masstau, massh, massw, massz, mUV;

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

S D, ep(:3);
V energyselector,p1,...,p40,ps1,...,ps40,k1,...,k40,c1,...,c40,cs1,...,cs40,fmb1,...,fmb40,fmbs1,...,fmbs40; * force this internal ordering in FORM
Auto V p,k,c;
Auto S lm,ext;
#ifdef `FOURDIM'
    Auto I mu=4,s=4;
#else
    Auto I mu=D,s=D;
#endif
Symbol ge, gs, ghhh, type, in, out, virtual;
Auto S x, idx, t, n;

Set spatialparts:ps1,...,ps40,cs1,...,cs40,fmbs1,...,fmbs40;
Set fmbs: fmb1,...,fmb40;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, gammatrace(c), GGstring, NN, vector,g(s),delta(s),T, counter,color, prop, replace;
CF f, vx, vxs(s), uvx, vec, vec1;
CF subs, configurations, conf, cmb, cbtofmb, fmbtocb, diag, forestid, der, energy, spatial(s);
CF subgraph, uvconf, uvconf1, uvconf2, uvprop, uv, uvtopo, integrateduv;
CT gammatracetensor(c),opengammastring;

S UVRenormFINITE;
S ICT, mUV, logmu, logmUV, logmt, mi1L1, alarmMi1L1;
Fill logmasses(6) = logmt;
Fill logmasses(-6) = logmt;

CF integratedct, rat, num, den;
Set ts: t0,...,t20;
CT penergy,energyct;
NF energync;
Symbol ca,cf,nf,[dabc^2],[d4RR],[d4RA],[d4AA];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli40;
Set colA: cOlj1,...,cOlj40;
Set colAdum: cOljj1,...,cOljj40;

Polyratfun rat;

*--#] setup :

#include tensorreduce.frm
#include integrateduv.frm

* Load the diagrams
#include- input_`SGID'.h

.sort:load-input;
Hide CONF;

*--#[ feynman-rules :

************************************************
* Substitute the Feynman rules for the numerator
************************************************

#procedure FeynmanRulesGlobal()
* extract the global factors from the Feynman rules, including colour

* make a copy of the Feynman rules
id prop(?a) = prop(?a)*tmp(prop(?a));
id vx(?a) = vx(?a)*tmp(vx(?a));
repeat id tmp(x1?)*tmp(x2?) = tmp(x1*x2);

* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* do the spin sum external particles
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = 1;
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = 1;
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = 1;
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = 1;
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = 1;

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = i_;
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = i_ * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = - i_ * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(`PSI', virtual, p?, idx1?, idx2?) = -i_;
id prop(`PSI', in, p?, idx1?) = 1;
id prop(`PSI', out, p?, idx1?) = 1;

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * cOlf(colA[idx3], colA[idx2], colA[idx1]) * (1/2);
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_* d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_;
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_;
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ghhh * i_;

id vx(`H', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = - i_ * d_(colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );
id vx(`H', `GLU', `GLU', `GLU', p4?, p1?, p2?, p3?, idx4?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );

#do i=3,6
    id vx(<x1?{`PSI',}>,...,<x`i'?{`PSI',}>, p1?, ...,p`i'?, idx1?, ..., idx`i'?) = (-1*i_)^(`i'-2);
#enddo

* delta_Z vertex

* The first multiplicity factor is always the loop multiplicity factor! It must be adjusted w.r.t to n_f!

* dZ massless quark
id vx(x1?{`QBARMASSLESS'}, x2?{`QMASSLESS'}, p1?, p2?, idx1?, idx2?) = (1/1) * (-1) * i_ * ((4/3)*gs^2/16/pi^2) * (1/ep) * d_(colF[idx1], colF[idx2]);

* the finite part needs to be checked, also because the factor 4/3 on the pole of the mass correction is pure fudge for now.
* dZ massive quark
id vx(x1?{`QBARMASSIVE'}, x2?{`QMASSIVE'}, p1?, p2?, idx1?, idx2?) = (1/1) * (-1) * i_ * ((4/3)*gs^2/16/pi^2) * d_(colF[idx1], colF[idx2]);

* dZ gluon

* The version below is for contributions to the gluon wavefunction from g, gh and down quark only, so it is good for e+ e- > j j j / u c s b t
id vx(`GLU', `GLU', p1?, p2?, idx1?, idx2?) = (1/3) * (-1) * i_ * d_(colA[idx1], colA[idx2]) * (gs^2/16/pi^2);

id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]);

* For the quartic gluon vertex we need an extra dummy index
* We also split up the vertex in 3 components with different colour
Multiply counter(1);
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(i?) = - counter(i + 1) * gs^2 * i_ *(
    + cOlf(colAdum[i], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[i])
        * vx(`GLU', `GLU', `GLU', `GLU', 1, p1, p2, p3, p4, idx1, idx2, idx3, idx4)
    + cOlf(colAdum[i], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[i])
        * vx(`GLU', `GLU', `GLU', `GLU', 2, p1, p2, p3, p4, idx1, idx2, idx3, idx4)
    + cOlf(colAdum[i], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[i])
        * vx(`GLU', `GLU', `GLU', `GLU', 3, p1, p2, p3, p4, idx1, idx2, idx3, idx4)
);
repeat id vx(`H', `GLU', `GLU', `GLU', `GLU', p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?)*counter(i?) = - counter(i + 1) * gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) * (
    + cOlf(colAdum[i], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[i])
        * vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4)
    + cOlf(colAdum[i], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[i])
        * vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4)
    + cOlf(colAdum[i], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[i])
        * vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4)
);

id counter(x?) = 1;

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id tmp(x?) = x;
.sort:feynman-rules-global;

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, ?a, cOli1?) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color(1);
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
repeat id cOlf(cOlj1?,cOlj2?,cOlj3?)*color(x?) = color(x * cOlf(cOlj1,cOlj2,cOlj3));
repeat id cOlNA*color(x?) = color(x * cOlNA);
repeat id cOlNR*color(x?) = color(x * cOlNR);

B+ color;
.sort:color-prep;
Keep brackets;

* Only evaluate this part per unique color by bracketing
Argument color;
    #call color
    #call simpli
    id  cOlI2R = cOlcR*cOlNR/cOlNA;
    id  cOlNR/cOlNA*cOlcR = cOlI2R;
    id  cOld33(cOlpR1,cOlpR2) = [dabc^2];
    id  cOlNR/cOlNA = nf/cf/2;
    id  cOlcR = cf;
    id  cOlcA = ca;
    id  cOlI2R = nf/2;
	id	cOld44(cOlpA1,cOlpA2) = [d4AA];
	id	cOld44(cOlpR1,cOlpR2) = [d4RR];
	id	cOld44(cOlpR1,cOlpA1) = [d4RA];

* set the SU(3) values
    id [dabc^2] = 15/18;
    id [d4AA] = 135;
    id [d4RR] = 5/12;
    id [d4RA] = 15/2;
    id cOlNR = 3;
    id cOlNA = 8;
    id cf = 4 / 3;
    id ca = 3;
    id nf = 1;
EndArgument;
.sort:color;

* set the SU(3) values
id cOlNR = 3;

id color(x?) = x;

#endprocedure

#procedure FeynmanRulesMomentum()
* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pdgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);

* do the spin sum external particles
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = 1;
id prop(`PHO', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = 1;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(`PSI', virtual, p?, idx1?, idx2?) = 1;
id prop(`PSI', in, p?, idx1?) = 1;
id prop(`PSI', out, p?, idx1?) = 1;

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) ;
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = (p3(lorentz[idx2])-p1(lorentz[idx2]));
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = 1;

id vx(`H', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = p3(lorentz[idx2])*p2(lorentz[idx3]) - p2.p3 * d_(lorentz[idx2], lorentz[idx3]);
id vx(`H', `GLU', `GLU', `GLU', p4?, p1?, p2?, p3?, idx4?, idx1?, idx2?, idx3?) =
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2])
;

#do i=3,6
    id vx(<x1?{`PSI',}>,...,<x`i'?{`PSI',}>, p1?, ...,p`i'?, idx1?, ..., idx`i'?) = 1;
#enddo

* delta_Z vertex

* The first multiplicity factor is always the loop multiplicity factor! It must be adjusted w.r.t to n_f!

* dZ massless quark
id vx(x1?{`QBARMASSLESS'}, x2?{`QMASSLESS'}, p1?, p2?, idx1?, idx2?) = gamma(dirac[idx1], p2, dirac[idx2]);

* the finite part needs to be checked, also because the factor 4/3 on the pole of the mass correction is pure fudge for now.
* dZ massive quark
id vx(x1?{`QBARMASSIVE'}, x2?{`QMASSIVE'}, p1?, p2?, idx1?, idx2?) =
      (1/ep + UVRenormFINITE*(4 + 3*(logmu - logmasses(x1))) ) * ( -gamma(dirac[idx1], p1, dirac[idx2]) - masses(x1) * gamma(dirac[idx1], dirac[idx2]) )
    + (-3/ep + UVRenormFINITE*(-4 - 3*(logmu - logmasses(x1))) ) * masses(x1) * gamma(dirac[idx1], dirac[idx2]);

* dZ gluon

* The version below is for contributions to the gluon wavefunction from g, gh and down quark only, so it is good for e+ e- > j j j / u c s b t
id vx(`GLU', `GLU', p1?, p2?, idx1?, idx2?) = (
    p1(lorentz[idx1]) * p1(lorentz[idx2]) * (
* gluon contribution
        ( (-11)*(1/ep) )
* ghost contribution
      + ( (-1/2)*(1/ep) )
    )
    - (p1.p1) * d_(lorentz[idx1], lorentz[idx2]) * (
* gluon contribution
        ( (-19/2)*(1/ep) )
* ghost contribution
      + ( (-1/2)*(1/ep) )
    )
* one massless quark contribution
    +(p1(lorentz[idx1]) * p1(lorentz[idx2]) - (p1.p1) * d_(lorentz[idx1], lorentz[idx2])) * (
        (1)*( (+4/3)*(1/ep) )
    )
);

id D^n? = rat(D^n, 1);
.sort:feynman-rules-vertices-1;

* construct gamma string, drop odd-length gamma traces and symmetrize the trace
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
id gamma(s1?,?a,s1?) = gammatrace(?a)*delta_(mod_(nargs_(?a), 2));

.sort:gamma-filter;

* TODO: use momentum conservation to reduce the number of different terms
#do i=1,1
    id once ifnomatch->skip vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) =
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2]);

    redefine i "0";
    label skip;

    B+ vx;
    .sort:3g;
    Keep brackets;
#enddo

id vx(`GLU', `GLU', `GLU', `GLU', 1, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 2, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 3, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

id vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id D^n? = rat(D^n, 1);
.sort:feynman-rules-vertices-2;

Multiply counter(1);
repeat id gammatrace(?a,p?,?b)*counter(i?) = vec(p,lorentzdummy[i])*gammatrace(?a,lorentzdummy[i],?b)*counter(i+1);
id counter(n?) = 1;

id gammatrace(?a) = gammatracetensor(?a);
id vec(p?,mu?) = p(mu);
Multiply replace_(D, 4-2*ep);

#ifdef `FOURDIM'
    .sort:fourdim-poly;
    Polyratfun;
    id rat(x1?,x2?) = x1/x2;
#endif

B+ gammatracetensor;
.sort:gamma-to-tensor;
Keep brackets;

#ifdef `FOURDIM'
    #do i=1,10
        id once gammatracetensor(?a) = g_(`i',?a);
        trace4 `i';
        B+ gammatracetensor;
        .sort:trace-`i';
        Keep brackets;
    #enddo
    .sort:fourdim-trace;
    Polyratfun rat;
#else
* at this stage all indices should be inside the gammatracetensor only
* FIXME: not the case anymore!
    #call Gstring(gammatracetensor,1)
    Multiply replace_(D, 4-2*ep);
#endif

id pzero = 0; * Substitute the 0-momentum by 0
.sort:feynman-rules-final;
#endprocedure

#procedure FeynmanRules()
    #call FeynmanRulesGlobal()
    #call FeynmanRulesMomentum()
#endprocedure

#procedure ExtractMomenta()
* strip momentum tags
    repeat id f?{vx,prop}(?a,p?) = f(?a);

* extract momenta from the Feynman rules
    id prop(x1?,x2?,?a,p1?,x3?,?b) = prop(x1,x2,?a,p1,x3,?b,conf(?a,p1));
    id vx(?a,x1?,p1?,?b,p2?,x3?,?c) = vx(?a,x1,p1,?b,p2,x3,?c,conf(p1,?b,p2));
    argument prop,vx;
        splitarg conf;
        repeat id conf(?a,-p?vector_,?b) = conf(?a,p,?b);
        repeat id conf(?a,p?,?b,p?,?c) = conf(?a,p,?b,?c);
    endargument;
    id f?{prop,vx}(?a,conf(?b)) = f(?a,?b);
#endprocedure
*--#] feynman-rules :

* Fill in the Feynman rules for the parts that do not depend on kinematics
#call FeynmanRulesGlobal()

* If the expression is empty due to color, we still write a file
#if ( termsin(F) == 0 )
    #write<out_`SGID'.proto_c> "#0 due to color\n"
    #write<out_integrand_`SGID'.proto_c> "#0 due to color\n"
    setexitflag;
#endif

.sort:empty-expression-check;

* process all the configurations
L F = F * CONF;
.sort:load-conf;
Drop CONF;

* TODO: needs to be moved
#if 0
* multiply the numerator contribution of derivatives
* FIXME: this is not global
id conf(?a,p?) = conf(?a) * penergy(p);
id conf(?a,x?) = conf(?a) * x; * note the type difference

AB+ forestid,subs,conf,cmb;
.sort:bubble-treatment-prep;
Keep brackets;

* Take a single derivative in pbubble_i^0
* FIXME: if the number of terms is 0 after differentiating, this configuration will not be created
repeat;
    id ifnomatch->bubbleend once der(p?) = replace_(p, p + pzero * x);
    if (count(x, 1) != 1) Discard;
    id x = 1;
    id pzero.p? = penergy(p);
    id penergy(pzero) = 1; * d/dpbubble_i^0 pbubble_i^0 = 1
endrepeat;
label bubbleend;
.sort:bubble-treatment-derivative;

id subs(p1?,p2?) = replace_(p1, p2);
.sort:bubble-treatment-replacement;
#else
* FIXME: this messes with the power counting though...
id conf(?a,p?) = conf(?a);
id conf(?a,x?) = conf(?a);
id subs(p1?,p2?) = replace_(p1, p2);
#endif

* transform the graph to the cut basis
* store a copy for the integrand routine in cmb
id conf(x?,x1?,cmb(?a),?b) = conf(x,x1,?b)*replace(?a)*cmb(?a);

* gradually transform the basis to prevent term blow-up
#do i=1,1
    id replace = 1;
    if (count(replace,1)) redefine i "0";

    AB+ cmb;
    .sort:cmb-1;
    Keep brackets; * make sure cmb is not replaced

    id replace(p1?,p2?,?a) = replace_(p1,p2)*replace(?a);
#enddo

.sort:cmb-2;

* collect the edges and vertices belonging to the forests
#call ExtractMomenta()
id forestid(?a) = forestid(?a,1);
repeat id forestid(?a,k1?,?b,x?)*f?{prop,vx}(?c,k1?,?d) = forestid(?a,k1,?b,x*f(?c,k1,?d));
repeat id cmb(?a)*forestid(?b,x?) = f(forestid(?b,x*cmb(?a)))*cmb(?a);
id cmb(?a) = 1;
id f(?a) = forestid(?a);
argtoextrasymbol tonumber,forestid,1;
#redefine tensorforeststart "`extrasymbols_'"

.sort:tensorforest-splitoff-1;
#redefine tensorforestend "`extrasymbols_'"
#do ext={`tensorforeststart'+1},`tensorforestend'
    L tensorforest{`ext'-`tensorforeststart'} = extrasymbol_(`ext');
#enddo
#define tensorforestcount "{`tensorforestend'-`tensorforeststart'}"

id forestid(x1?,?a,x3?) = forestid(?a)*forest(x1)*x3*conf(-1,x1);

.sort:tensorforest-splitoff-2;
Hide F;

id forestid(?a) = 1;
id conf(?a) = 1;

* Apply cmb to the uv subgraph
id cmb(?a) = cmb(?a)*replace(?a);
AB+ cmb;
.sort:cmb-replace-1;
Keep brackets;
id replace(?a) = replace_(?a);
.sort:cmb-replace-2;

AB+ cmb;
.sort:fmb-replace-1;
Keep brackets;
* convert to the forest mb
id cbtofmb(?a) = replace_(?a);
.sort:fmb-replace-2;

* get the momentum dependence in the fmb
#call ExtractMomenta()

* move the vertices and propagators into a subgraph, starting from the most nested ones first
id subgraph(?a,uvconf(?b,x?,x1?)) = subgraph(?a,uvconf(?b,x,x1),?b); * prevent pattern matching issues by lowering the ?b by one level
#do i=1,5
    repeat id subgraph(n1?,...,n`i'?,uvconf(?a,x1?),?b,fmb1?,?c)*f?{prop,vx}(?d,fmb1?,?e) = subgraph(n1,...,n`i',uvconf(?a,x1*f(?d,fmb1,?e)),?b,fmb1,?c);
#enddo
id subgraph(?a,uvconf(?b),?c) = subgraph(?a,uvconf(?b));

.sort:graphs;

* Create the local UV counterterm
#define uvdiagtotalstart "`extrasymbols_'"
#do i = 1,1
* Split off all UV subgraphs without dependencies to a separate expression
    id subgraph(x1?, x2?) = -uvconf1(x1, x2);
    argtoextrasymbol tonumber,uvconf1,2;
    #redefine uvdiagstart "`extrasymbols_'"
    .sort:uvdiag-1;
    #redefine uvdiagend "`extrasymbols_'"
    #do ext={`uvdiagstart'+1},`uvdiagend'
        L uvdiag`ext' = extrasymbol_(`ext');
    #enddo
    .sort:uvdiag-2;

    id uvconf(?a,x2?,x3?) = tmax^x2*x3;

* fill in results from subdiagrams
    #do ext={`uvdiagtotalstart'+1},`uvdiagstart'
        id uvconf1(`ext') = uvdiag`ext';
    #enddo
    .sort:uv-subgrap-fillin;
    Hide tensorforest1,...,tensorforest`tensorforestcount';

* Apply Feynman rules to the UV subgraph
    #call FeynmanRulesMomentum()

* linearize the gamma matrices and convert them to tensors
* this makes them suitable for differentiation
    repeat id once gamma(s1?,?a,p?!vector_,?b,s2?) = p(mudummy)*gamma(s1,?a,mudummy,?b,s2);
    id gamma(?a) = opengammastring(?a);

    id uvprop(k?,t1?,p?,m?) = uvprop(k,t1,p,m)*uvconf2(p);

    splitarg uvconf2;
    chainout uvconf2;
    id uvconf2(-p?vector_) = uvconf2(p);
    repeat id uvconf2(p?)*uvconf2(p?) = uvconf2(p);

    id uvconf2(p?) = replace_(p, t * p);

    argument uvprop,1,vxs,uvtopo,diag;
        id t = 1;
    endargument;

* rescale all masses in the numerator
    repeat id m?allmasses = tmp(m);
    id tmp(m?) = t*m;

* Taylor expand the propagators to the right depth
* p carries a t dependence that determines the order
* t1 determines the powers of the UV propagator

* expand the propagators without loop momentum dependence
    id uvprop(k?,t1?,0,m?) = uvprop(k,t1,1,m) * (1 - (mUV^2*t^2-m^2*t^2) * t1 + (mUV^2*t^2-m^2*t^2)^2 * t1^2 + ALARM * t^5);
    id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
    repeat;
        id once ifnomatch->skiptruncation uvprop(k?,t1?,p?,m?)*t^x1?*tmax^x2? = uvprop(k,t1,1,m) * t^x1*tmax^x2 * theta_(x2-x1) *
            (1 +
                (-2*p.k-(p.p+mUV^2*t^2-m^2*t^2)) * t1 +
                (+4*p.k^2+4*p.k*(p.p+mUV^2*t^2-m^2*t^2)+(p.p+mUV^2*t^2-m^2*t^2)^2) * t1^2 +
                (-8*p.k^3-12*p.k^2*(p.p+mUV^2*t^2-m^2*t^2)) * t1^3 +
                (16*p.k^4) * t1^4 +
                ALARM * t^5);
        id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
        label skiptruncation;
    endrepeat;

    id t = 1;
    id tmax = 1;

    if (count(ALARM, 1));
        Print "UV Taylor expansion depth exceeded.";
        exit "";
    endif;

* match the denominator structure to a diagram
    id uvprop(?a,m?) = uvprop(?a);
    repeat id t1?ts^n2? = tmp(t1,n2);
    repeat id uvprop(k?,t1?,n1?)*uvprop(k?,t1?,n2?) = uvprop(k,t1,n1+n2); 
    repeat id uvprop(k?,t1?,n1?)*uvtopo(x?,x1?,?a)*tmp(t1?,n2?) = uvprop(k,n1 + n2)*uvtopo(x,x1*t1^n2,?a);
    id uvprop(k?,t1?ts,n?) = uvprop(k, n);

    #call uvmap()
    if (count(uvtopo, 1, tmp, 1));
        Print "Unsubstituted UV topology: %t";
        exit "Critical error";
    endif;
    
* collect all uv propagators of the subgraph
    Multiply replace_(vxs, uvx);
    if (count(ICT, 1) == 0);
        id uvprop(?a) = 1;
        id uvx(?a) = 1;
    else;
        id 1/ICT = ICT;
    endif;

    chainin uvprop;
    id uvprop(?a) = uvprop(?a, 1);
    repeat id uvx(?a)*uvprop(?b,x?) = uvprop(?b, x*uvx(?a));
    Multiply replace_(uvx, vxs);
    id uvprop(?a) = integrateduv(?a);


* TODO: do integrated counterterm here
#if 0
* compute the integrated UV counterterm
    Multiply counter(1);
    repeat id k1?.k2?*counter(n?) = vec(k1,n)*vec(k2,n)*counter(n + 1);
* convert every k^0 into k.p0select, where p0select is effectively (1,0,0,0)
    repeat id penergy(k1?)*counter(n?) = vec(k1,n)*vec(p0select,n)*counter(n + 1);
    id counter(x?) = 1;

    #do j=1,1
        id ICT = -1; * we add back the counterterm

        id once integrateduv(?a,x?) = uvprop(?a)*x;
        if (count(uvprop,1)) redefine j "0";

        Multiply replace_(vec, vec1); * consider all vectors as external
        repeat id uvprop(k1?,n?)*vec1(k1?,n1?) = uvprop(k1,n)*vec(k1,n1);

* TODO: also drag in gamma matrices in the projection

        #call TensorReduce()

* contract all metrics
        repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);
        id g(n1?,n1?) = rat(4-2*ep,1);
        repeat id vec1(p1?,n1?)*g(n1?,n2?) = vec1(p1,n2);

        .sort:tensor-projection-loop;

        if (count(uvprop,1));
* divide by the normalizing factor of the denominator that is added to the topology
* this is always 1/(k^2 - m_UV^2)^3 = -i / (4 pi)^2 * 1/2 * 1/mUV^2
            Multiply i_ * (4 * pi)^2 * 2 * mUV^2;

            #call IntegrateUV()
        endif;

        Multiply replace_(D, 4 - 2 * ep);
        id ep^n1? = rat(ep^n1,1);

        .sort:ibp-reduction;
    #enddo

    id vec1(k1?,n?)*vec1(k2?,n?) = k1.k2;
    id k1?.p0select = penergy(k1);

* Substitute the masters and expand in ep
    #call SubstituteMasters()
    .sort:integrated-ct-1;

    PolyRatFun;
    id rat(x1?) = x1;
    if (count(ep, 1) != `SELECTEDEPSILONORDER') Discard; * keep only the ep^0 piece
    id ep^n? = 1;
    #if `UVRENORMFINITEPOWERTODISCARD' > 0
        if (count(UVRenormFINITE, 1) >= `UVRENORMFINITEPOWERTODISCARD') Discard; * Discard UVRenormFinite pieces up to some order
    #elseif  `UVRENORMFINITEPOWERTODISCARD' < 0
        if (count(UVRenormFINITE, 1) < -`UVRENORMFINITEPOWERTODISCARD') Discard; * Keep UVRenormFinitiePieces up to some order
    #endif
    id UVRenormFINITE^n? = 1;

    .sort:integrated-ct-2;
#endif

    .sort:uv-subgraph-done;
    UnHide tensorforest1,...,tensorforest`tensorforestcount';
    Hide uvdiag{`uvdiagstart'+1},...,uvdiag`uvdiagend';

* in the next loop the UV subgraph will be added to the tensorforest or to a UV graph that embeds it
    if (count(uvconf1, 1)) redefine i "0";

* fill in the unsubstituted subgraph into the supergraph
    repeat id subgraph(x1?,?a,n?,?b,uvconf(?c,x2?))*uvconf1(n?,x3?) = subgraph(x1,?a,?b,uvconf(?c,x2*uvconf1(x3)));
    id uvconf1(n?,x?) = uvconf1(x);
    .sort:uv3;
#enddo

.sort:local-uv-done;
Drop uvdiag`uvdiagtotalstart',...,uvdiag`uvdiagend';
delete extrasymbols>`uvdiagtotalstart';

if (count(subgraph, 1));
    Print "Unsubstituted UV subgraph: %t";
    exit "Critical error";
endif;

* Apply Feynman rules to remaining graph
#call FeynmanRulesMomentum()

* Simplify all open gamma strings
#call Gstring(opengammastring,0)

* split off the energy part: energyselector=(1,0,0,0,..), fmbs = spatial part
* fmbs1.fmbs2 will later get a minus sign to honour the Minkowksi metric
AB+ energy,diag,fmbtocb;
.sort:energy-splitoff-1;
Keep brackets;
#do i=1,40
    Multiply replace_(fmb`i', energy(fmb`i')*energyselector - fmbs`i');
#enddo
* linearize the gamma string to extract the energy
repeat id once gamma(s1?,?a,p?!vector_,?b,s2?) = p(mudummy)*gamma(s1,?a,mudummy,?b,s2);
id energyselector.energyselector = 1;
id energyselector.p?spatialparts = 0;
.sort:energy-splitoff-2;

* TODO: study grouping and effects on term size
Multiply f(1);
repeat id f(x?)*f1?{cmb,diag,energy}(?a) = f(x*f1(?a));
.sort:tensorforest-0;

argtoextrasymbol tonumber,f,1;
#redefine foreststart "`extrasymbols_'"

.sort:tensorforest-1;
#redefine forestend "`extrasymbols_'"
#do ext={`foreststart'+1},`forestend'
    L forest{`ext'-`foreststart'} = extrasymbol_(`ext');
#enddo
#define forestcount "{`forestend'-`foreststart'}"

* convert fmb spatial part back to cmb in the forest
id fmbtocb(?a) = replace_(?a);

.sort:tensorforest-2;
Drop tensorforest1,...,tensorforest`tensorforestcount';
Hide forest1,...,forest`forestcount';
UnHide F;

* fill in the tensor forest into F and fill in the Feynman rules
#do i=1,`forestcount'
    id forestid(`i') = tensorforest`i';
#enddo
id f(x?) = forestid(x-`foreststart');
.sort:tensorforest-3;

id opengammastring(?a) = gamma(?a);

#call FeynmanRulesMomentum()

if (count(gamma,1));
    Print "Unsubstituted gamma: %t";
    exit "Critical error";
endif;

* TODO: keep the spatial vectors since it's faster
id energyselector.energyselector = 1;
id energyselector.p?spatialparts = 0;
id p?.energyselector = penergy(p);
id p1?spatialparts.p2?spatialparts = -p1.p2; * add a -1 to fix the metric
id p1?spatialparts.p? = spatial(p1,p);
argument spatial;
    Multiply replace_(<ps1,p1>,...,<ps40,p40>,<cs1,c1>,...,<cs40,c40>);
endargument;

.sort:feynman-complete-graph;

* Expand
Multiply replace_(D, 4 - 2 * ep);
id ep^n1? = rat(ep^n1,1);
.sort:rat-expand;
PolyRatFun rat(expand,ep,{`MAXPOLE'+`SELECTEDEPSILONORDER'});
Keep brackets;
.sort:integrated-ct-1;

PolyRatFun;
id rat(x1?) = x1;
if (count(ep, 1) != `SELECTEDEPSILONORDER') Discard; * keep only the ep^0 piece
id ep^n? = 1;
#if `UVRENORMFINITEPOWERTODISCARD' > 0
    if (count(UVRenormFINITE, 1) >= `UVRENORMFINITEPOWERTODISCARD') Discard; * Discard UVRenormFinite pieces up to some order
#elseif  `UVRENORMFINITEPOWERTODISCARD' < 0
    if (count(UVRenormFINITE, 1) < -`UVRENORMFINITEPOWERTODISCARD') Discard; * Keep UVRenormFinitiePieces up to some order
#endif
id UVRenormFINITE^n? = 1;
.sort:rat-truncate;

* If the expression is empty (due to epsilon pole selection), we still write a file
#if ( termsin(F) == 0 )
    #write<out_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
    #write<out_integrand_PF_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
    #write<out_integrand_LTD_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
#endif

id cmb(?a) = cmb(?a)*replace(?a);
AB+ cmb;
.sort:cmb-replace;
Keep brackets;
id replace(?a) = replace_(?a);
.sort:pf-splitoff;
UnHide forest1,...,forest`forestcount';
Hide F;

id rat(x1?,x2?) = x1/x2;

* collect all the energies in the diagram
repeat id diag(?a,p1?,p2?,?b) = diag(?a,p1,0,p2,?b);
id diag(?a,p1?) = diag(?a,p1,0);

repeat id energy(k?)*diag(?a,k?,n?,?b) = diag(?a,k,n+1,?b);

if (count(energy,1));
    Print "Energy left: %t";
    exit "Critical error";
endif;

if (count(cmb,1) == 0);
    Print "CMB missing: %t";
    exit "Critical error";
endif;
.sort:energy-collect;

*********************************************
* Construction of the integrand
*********************************************

#if (`INTEGRAND' == "LTD") || (`INTEGRAND' == "both")
    Hide forest1,...,forest`forestcount';

* copy the forest since LTD and PF may map the diag calls differently
* TODO: unify call
    #do i=1,`forestcount'
        L forestltd`i' = forest`i';
    #enddo

    .sort:ltd-forest-copy;
    #include- ltdtable_`SGID'.h

* map all diagrams to their unique representative
    id diag(x1?,x2?,?a) = diag(x1,ltdmap(x1,x2),?a);
    id diag(x1?,diag(?a),?b) = diag(x1,?a,?b);

    repeat id cmb(?a)*diag(?b) = f(diag(?b,cmb(?a)))*cmb(?a);
    id f(?a) = diag(?a);
    argtoextrasymbol tonumber,diag,1;
    #redefine diagstart "`extrasymbols_'"

    .sort:diag-1;
    #redefine diagend "`extrasymbols_'"
    #do ext={`diagstart'+1},`diagend'
        L diag{`ext'-`diagstart'} = extrasymbol_(`ext');
    #enddo
    #define diagcount "{`diagend'-`diagstart'}"

    id diag(x?) = diag(x-`diagstart'+`forestcount');
    id diag(xcut?,x1?,x2?,?a,x3?) = diag(?a)*ltdtopo(x1,x2)*x3*conf(-1,xcut,-1);
    id cmb(?a) = replace_(?a);
    .sort:ltd-splitoff;
    Hide forestltd1,...,forestltd`forestcount';

* transform the LTD energies back to normal energies
    id energy(f(?a)) = energy(?a);
    chainout energy;
    id energy(c0) = 1;

    repeat id ltdcbtolmb(?a,fmb1?,x?,?b)*diag(fmb1?,n?,?c) = x^n*ltdcbtolmb(?a,fmb1,x,?b)*diag(?c);
    id diag = 1;
    id ltdcbtolmb(?a) = 1;

    if (count(diag,1));
        Print "Unsubstituted energies: %t";
        exit "Critical error";
    endif;

    .sort:energy-replacement;

    id prop(?a) = prop(-1,?a);
    repeat id prop(x1?,x?,?a)*prop(x2?,x?,?a) = prop(x1+x2,x,?a);

    #do i = 1,`NFINALMOMENTA'
        id once der(ltd1?,n?,?a) = f(ltd1,n)*der(?a)/fac_(n);
* write the numerator as a propagator
        id ltd1?^n?*f(ltd1?,n1?) = prop(n,-1,1,ltd1,0,0)*f(ltd1,n1);
        repeat id prop(x1?,x?,?a)*prop(x2?,x?,?a) = prop(x1+x2,x,?a);
* perform the derivative
        repeat id all prop(n2?,x1?,?a,n1?,ltd1?,?b)*f(ltd1?,n3?{>0}) = n2*n1*prop(n2-1,x1,?a,n1,ltd1,?b)*f(ltd1,n3 - 1);

        id prop(n?{>=0},-1,1,E?,0,0) = E^n;

        id f(ltd1?,0) = 1;

        if (count(f,1));
            Print "Derivative left in result: %t";
            exit "Critical error";
        endif;
        .sort:der-`i';
    #enddo
    id der = 1;

* rewrite the propagators
    id energies(0) = 0;
    id prop(n?,x?,?a) = x^(-n);

* set the ltd energies (including cut sign)
    id ltdenergy(?a) = replace_(?a);

    argument ellipsoids;
        id energies(0) = 0;
    endargument;

    .sort:ltd-num-0;
    UnHide forestltd1,...,forestltd`forestcount';
    .sort:ltd-num-1;
    Drop diag1,...,diag`diagcount',forestltd1,...,forestltd`forestcount';

* now add all LTD structures as a special conf
    #if `diagcount' > 0
        L FINTEGRANDLTD = F + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})> + <forestltd1*conf(-1,-1)>+...+<forestltd`forestcount'*conf(-`forestcount',-1)>;
    #endif

    id conf(x?)*conf(x1?{<0},?a) = conf(x,?a);

* TODO: deprecated
    #if (`SUMDIAGRAMSETS' == "onlysum")
        id conf(x?{>=0},x1?,?a) = conf(1000 + x1);
    #elseif (`SUMDIAGRAMSETS' == "both")
        id conf(x?{>=0},x1?,?a) = conf(x,?a) + conf(1000 + x1);
    #endif

    .sort:integrand-ltd;
    #if (`INTEGRAND' == "both")
        Hide FINTEGRANDLTD;
        delete extrasymbols>`diagstart';
    #endif
#endif

#if (`INTEGRAND' == "PF") || (`INTEGRAND' == "both")
    .sort:pf-start-1;
    UnHide forest1,...,forest`forestcount';
    .sort:pf-start-2;

    #include- pftable_`SGID'.h
    .sort:load-pf;

* map all diagrams to their unique representative
    id diag(x1?,x2?,?a) = diag(x1,pfmap(x1,x2),?a);
    id diag(x1?,diag(?a),?b) = diag(x1,?a,?b);

* TODO: here we create a unique diagram per cut as we include the cmb (and conf) into the key
* it is possible that diagrams are the same after the cmb is applied, for example for UV topologies
    repeat id cmb(?a)*diag(?b) = f(diag(?b,cmb(?a)))*cmb(?a);
    id f(?a) = diag(?a);
    argtoextrasymbol tonumber,diag,1;
    #redefine diagstart "`extrasymbols_'"

    .sort:diag-1;
    #redefine diagend "`extrasymbols_'"
    #do ext={`diagstart'+1},`diagend'
        L diag{`ext'-`diagstart'} = extrasymbol_(`ext');
    #enddo
    #define diagcount "{`diagend'-`diagstart'}"

    id diag(x?) = diag(x-`diagstart'+`forestcount');
    id diag(xcut?,x1?,x2?,?a,x3?) = diag(?a)*pftopo(x1,x2)*x3*conf(-1,xcut,-1);
    id cmb(?a) = replace_(?a);
    .sort:pf-splitoff;
    Hide forest1,...,forest`forestcount';

* apply the numerator procedure
    #do i = 0,{`NFINALMOMENTA'-1}
        id diag(p?,n?,?a) = diag(?a)*x1^n*f(p);
        repeat id energy(p?)*f(p?) = x1*f(p); * energies could come from previous iteration
        id f(p?) = 1;
        id once num(ncmd(?z),?x) =ncmd(?z)*num(?x);
        B+ ncmd, x1;
        .sort:collect-ncmd;
        keep brackets;

        id x1^r?*ncmd(?z,x?) = sum_(s,0, r - nargs_(?z), a(r,s,?z)*x^s);
        B+ a;
        .sort:collect-a;
        keep brackets;

        repeat id a(r?,s?,?z,x?) = -sum_(s1, s+1,r-nargs_(?z), a(r,s1,?z)*x^(s1-s-1));
        id a(r?,s?) = delta_(r,s);
        .sort:energy-`i';
    #enddo
    id diag = 1;
    id num = 1;

* check that the substitution is complete
    if (count(ncmd, 1,num,1,a, 1, energy, 1, diag, 1));
        Print "Unsubstituted ncmd: %t";
        exit "Critical error";
    endif;

    .sort:pf-num-1;

    UnHide forest1,...,forest`forestcount';
    .sort:pf-num-2;
    Drop diag1,...,diag`diagcount',forest1,...,forest`forestcount';

* now add all PF structures as a special conf
    #if `diagcount' > 0
        L FINTEGRANDPF = F + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})> + <forest1*conf(-1,-1)>+...+<forest`forestcount'*conf(-`forestcount',-1)>;
    #endif

    id conf(x?)*conf(x1?{<0},?a) = conf(x,?a);

* TODO: deprecated
    #if (`SUMDIAGRAMSETS' == "onlysum")
        id conf(x?{>=0},x1?,?a) = conf(1000 + x1);
    #elseif (`SUMDIAGRAMSETS' == "both")
        id conf(x?{>=0},x1?,?a) = conf(x,?a) + conf(1000 + x1);
    #endif

    .sort:pf-integrand-collect;
    #if (`INTEGRAND' == "both")
        UnHide FINTEGRANDLTD;
    #endif
#endif

* fill in the shifts
id replace(?a) = replace_(?a);
id energy(p?) = penergy(p);
id energies(p?) = penergy(p);

B+ penergy,spatial,energies,allenergies, ellipsoids, constants;
.sort:func-prep;
Keep brackets;

argument ellipsoids;
    id energies(p?) = penergy(p);
endargument;

repeat id allenergies(x?,p?,m?,?a) = energync(x,p.p+m*m)*allenergies(?a);
argument energync, 2;
    id p1?.p2? = spatial(p1, p2);
endargument;
chainin energync;
id energync(?a)*allenergies = allenergies(?a);
id allenergies(?a) = energies(?a);

repeat id constants(p?,m?,?a,x?) = energync(p.p-m*m)*constants(?a,x);
chainin energync;
id energync(?a)*constants(x?) = constants(?a,x);

argument ellipsoids, constants;
    id energy(p?) = penergy(p);
endargument;

* Convert the dot products and energies to a symbol
#$MAXK = `NFINALMOMENTA';
#$MAXP = `NINITIALMOMENTA';
#$OFFSET = 0;
#do i=1,`$MAXP'
    id penergy(p`i') = lm`$OFFSET';
    argument energies, ellipsoids, constants;
        id penergy(p`i') = lm`$OFFSET';
    endargument;
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        argument energies, ellipsoids, constants;
            id p`i'.p`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`i', p`j') = lm`$OFFSET';
        argument energies;
            id spatial(p`i', p`j') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXK'
    id penergy(c`i') = lm`$OFFSET';
    argument energies, ellipsoids, constants;
        id penergy(c`i') = lm`$OFFSET';
    endargument;
    #$OFFSET = $OFFSET + 1;
    #do j=1,`$MAXP'
        argument energies, ellipsoids, constants;
            id c`i'.p`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`j', c`i') = lm`$OFFSET';
        argument energies;
            id spatial(p`j', c`i') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=`i',`$MAXK'
        argument energies, ellipsoids, constants;
            id c`i'.c`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(c`i', c`j') = lm`$OFFSET';
        argument energies;
            id spatial(c`i', c`j') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

.sort:conv-func;

#$OFFSET = 0;
#do i=1,`$MAXP'
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 2;
    #enddo
#enddo

#do i=1,`$MAXK'
    #$OFFSET = $OFFSET + 1;
    #do j=1,`$MAXP'
        id c`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 2;
    #enddo

    #do j=`i',`$MAXK'
        id c`i'.c`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 2;
    #enddo
#enddo
.sort:conv-dots;

* split off every energy configuration into a new expression
B+ conf;
.sort:conf-collect;
Keep brackets;

#do INTEGRANDTYPE = {PF,LTD}
    #if (`INTEGRAND' == `INTEGRANDTYPE') || (`INTEGRAND' == "both")
        Hide;
        UNHide FINTEGRAND`INTEGRANDTYPE';
        NHide FINTEGRAND`INTEGRANDTYPE';

        B+ conf;
        .sort:conf-0;
        Keep brackets;

        id conf(?a) = conf(conf(?a));
        argtoextrasymbol tonumber,conf,1;
        #redefine oldextrasymbols "`extrasymbols_'"
        B+ conf;
        .sort:conf-1;
        Hide FINTEGRAND`INTEGRANDTYPE';
        #redefine energysymbolstart "`extrasymbols_'"
        #do ext={`oldextrasymbols'+1}, `energysymbolstart'
            #write "START integrand"
            #$conf = extrasymbol_(`ext');
            #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#CONF\n%$", $conf;
            L FF`ext' = FINTEGRAND`INTEGRANDTYPE'[conf(`ext')];

            $isdenominator = 0;
            inside $conf;
                if (match(conf(x?{<0},x1?,x2?))) $isdenominator = 1;
            endinside;
            $ellipsoids = 0;
            $energies = 0;
            $constants = 0;
            id ellipsoids(?a$ellipsoids) = 1;
            id energies(?a$energies) = 1;
            id constants(?a$constants) = 1;
            .sort:conf-`ext'-0;
            #if (`$isdenominator' == 1)
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#CONSTANTS\n%$\n#CONSTANTS", $constants
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#ENERGIES\n%$\n#ENERGIES", $energies
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#ELLIPSOIDS\n%$\n#ELLIPSOIDS",$ellipsoids

* FIXME: for now, we don't optimize the PF denominators
                Format C;
                Format O1,stats=on;
            #else
                Format C;
                #if `OPTIMLVL' > 1
                    Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
                #else
                    Format O1,stats=on;
                #endif
            #endif 

* Optimize the output
            #Optimize FF`ext'
            #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "%O"
            #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "\n\treturn %E;",FF`ext'
            .sort:optim-`ext'-1;

            #clearoptimize;
            .sort:optim-`ext'-2;
            Format O0;
            Format normal;
            Drop FF`ext';
            #write "END integrand"
        #enddo
        .sort:conf-2;
        Drop FINTEGRAND`INTEGRANDTYPE';
        delete extrasymbols>`oldextrasymbols';
    #endif
#enddo

.end
