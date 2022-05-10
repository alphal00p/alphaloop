#-
#:ContinuationLines 0
On statistics;
On nospacesinnumbers;

*--#[ setup :

#define PSI "3370,3371,3372,3373,3374,3375,3376,3377,3378,3379"
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
#define Z "23"

**************************************************
* START SE PDGs
**************************************************
#define PHOPRIME "1022"
#define QMASSIVEPRIME "1006,"
#define QBARMASSIVEPRIME "-1006,"
#define SDUMMY "1122"
**************************************************
* END SE PDGs
**************************************************

S vev, pi, cw, sw ,sw2 , gw;

Auto S mass;
Auto S yukawa;
CTable masses(-10000:10000);
CTable gyq(-10000:10000);
CTable logmasses(-10000:10000);
CTable charges(-10000:10000);
CTable zVcoupling(-10000:10000);
CTable zAcoupling(-10000:10000);

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
Fill masses(3370) = 0;
Fill masses(3371) = massdummya;
Fill masses(3372) = massdummyb;
Fill masses(3373) = massdummyc;
Fill masses(3374) = massdummyd;
Fill masses(3375) = massdummye;
Fill masses(3376) = massdummyf;
Fill masses(3377) = massdummyg;
Fill masses(3378) = massdummyh;
Fill masses(3379) = massdummyi;

**************************************************
* START EW parameters
**************************************************

Fill masses(23) = massz;

Fill zVcoupling(1) = -1/2+2/3*sw2; * d
Fill zVcoupling(2) = 1/2-4/3*sw2; * u
Fill zVcoupling(3) = -1/2+2/3*sw2; * s
Fill zVcoupling(4) = 1/2-4/3*sw2; * c
Fill zVcoupling(5) = -1/2+2/3*sw2; * b
Fill zVcoupling(6) = 1/2-4/3*sw2; * t
Fill zVcoupling(11) = -1/2+2*sw2; * e-
Fill zVcoupling(12) = -1/2+2*sw2; * mu-
Fill zVcoupling(13) = -1/2+2*sw2; * ta-

Fill zAcoupling(1) = -1/2; * d
Fill zAcoupling(2) = 1/2; * u
Fill zAcoupling(3) = -1/2; * s
Fill zAcoupling(4) = 1/2; * c
Fill zAcoupling(5) = -1/2; * b
Fill zAcoupling(6) = 1/2; * t
Fill zAcoupling(11) = -1/2; * e-
Fill zAcoupling(12) = -1/2; * mu-
Fill zAcoupling(13) = -1/2; * ta-

**************************************************
* END EW parameters
**************************************************

**************************************************
* START SE parameters
**************************************************
Fill masses(-1006) = masst;
Fill masses(1006) = masst;
Fill charges(-1006) = -2/3;
Fill charges(1006) = 2/3;
Fill gyq(-1006) = yukawat;
Fill gyq(1006) = yukawat;
**************************************************
* END SE parameters
**************************************************

* note: this set needs to be complete for the UV expansion
Set allmasses: massu, massd, massc, masss, masst, massb, masse, massmu, masstau, massh, massw, massz;

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

S D, ep(:{`MAXPOLE'+`SELECTEDEPSILONORDER'});
V energyselector,p1,...,p40,ps1,...,ps40,k1,...,k40,c1,...,c40,cs1,...,cs40,fmb1,...,fmb40,fmbs1,...,fmbs40; * force this internal ordering in FORM
Auto V p,k,c;
Auto S lm,ext,E;
#ifdef `FOURDIM'
    Auto I mu=4,s=4;
#else
* note: a change of spinor dimension here needs to be reflected in Gstring and in the Feynman rules
    Auto I mu=D,s=4;
#endif
Symbol ge, gs, ghhh, type, in, out, virtual;
Auto S x, idx, t, n;

Set spatialparts:ps1,...,ps40,cs1,...,cs40,fmbs1,...,fmbs40;
Set fmbs: fmb1,...,fmb40;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, gamma5, gammatrace(c), levicivita, GGstring, NN, vector,g(s),delta(s),tmps(s), counter,color, prop, replace;
CF dot, f, vx, vxs(s), uvx, vec, vec1;
CF subs, configurations, conf, tder, cmb, cbtofmb, fmbtocb, diag, forestid, nloops, der, energy, spatial(s), onshell;
CF subgraph, uvconf, uvconf1, uvconf2, uvprop, uv, uvtopo, irtopo, intuv, integrateduv;
CT gammatracetensor(c),opengammastring;

CF logmUVmu; * make it a function so that the optimizer makes sure it is only computed once
S UVRenormFINITE, ICT, mUV, logmu, logmUV, logmt, z3, mi1L1, alarmMi1L1, gamma0;
Fill logmasses(6) = logmt;
Fill logmasses(-6) = logmt;

**************************************************
* START SE parameters
**************************************************
Fill logmasses(1006) = logmt;
Fill logmasses(-1006) = logmt;
**************************************************
* END SE parameters
**************************************************

CF integratedct, rat, num, den, tmp;
Set ts: t0,...,t20;
CT penergy,energyct;
NF energync;
Symbol ca,cf,nf,[dabc^2],[d4RR],[d4RA],[d4AA];

S  i, m, n, ALARM;

CT colf,T,colTr;
Auto I coli=3,colj=8;
Set colF: coli1,...,coli40;
Set colA: colj1,...,colj40;

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

* construct a counter that yields new indices
id vx(?a,p?,idx1?,?b) = tmp(idx1,?b)*vx(?a,p,idx1,?b);
chainin tmp;
repeat id tmp(idx1?,idx2?,?a) = tmp(max_(idx1,idx2),?a);
id tmp(x?) = counter(x + 1);

* split up the quartic gluon vertex into distinct colour factors
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(idx5?) = counter(idx5 + 1) *(
    +vx(`GLU', `GLU', `GLU', `GLU', 1, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
    +vx(`GLU', `GLU', `GLU', `GLU', 2, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
    +vx(`GLU', `GLU', `GLU', `GLU', 3, p1, p2, p3, p4, idx1, idx2, idx3, idx4, idx5)
);
repeat id vx(`H', `GLU', `GLU', `GLU', `GLU', p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?)*counter(idx6?) = counter(idx6 + 1) * (
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
    +vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5, p1, p2, p3, p4, idx5, idx1, idx2, idx3, idx4, idx6)
);

repeat id vx(x1?, `Z', x2?, p1?, p2?, p3?, idx1?, idx2?, idx3?)*counter(idx4?) = counter(idx4 + 5) *
    vx(x1, `Z', x2, p1, p2, p3, idx1, idx2, idx3, idx4, idx4 + 1, idx4 + 2, idx4 + 3, idx4 + 4);

id counter(x?) = 1;

* make a copy of the Feynman rules
id prop(?a) = prop(?a)*tmp(prop(?a));
id vx(?a) = vx(?a)*tmp(vx(?a));
repeat id tmp(x1?)*tmp(x2?) = tmp(x1*x2);

* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* do the spin sum external particles
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = 1;
repeat id prop(`GLU', in, p?, idx1?)*prop(`GLU', out, p?, idx2?) = d_(colA[idx1], colA[idx2]);
repeat id prop(`Z', in, p?, idx1?)*prop(`Z', out, p?, idx2?) = 1;
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = 1;
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = 1;
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = d_(colF[idx1], colF[idx2]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_;
id prop(`Z', virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = i_;
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = - i_;
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = i_ * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = - i_ * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(x?{`PSI'}, virtual, p?, idx1?, idx2?) = -i_;
id prop(x?{`PSI'}, in, p?, idx1?) = 1;
id prop(x?{`PSI'}, out, p?, idx1?) = 1;

**************************************************
* START SE prop couplings Feynman rules
**************************************************
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = 1;
repeat id prop(x1?{`QBARMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVE'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2]);
repeat id prop(x1?{`QBARMASSIVE'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVEPRIME'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2]);
repeat id prop(x1?{`QMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QMASSIVE'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
repeat id prop(x1?{`QMASSIVE'}, in, p?, idx1?)*prop(x2?{`QMASSIVEPRIME'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
id prop(`SDUMMY', virtual, p?, idx1?, idx2?) = 1;
id prop(`SDUMMY', in, p?, idx1?) = 1;
id prop(`SDUMMY', out, p?, idx1?) = 1;
id prop(`PHOPRIME', virtual, p?, idx1?, idx2?) = 1;
id prop(x?{`QMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = i_ * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBARMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = - i_ * d_(colF[idx1], colF[idx2]);
**************************************************
* END SE prop couplings Feynman rules
**************************************************

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

**************************************************
* START SE vx couplings Feynman rules
**************************************************
id vx(`PHOPRIME', `PHOPRIME', `SDUMMY', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_;
id vx(x1?{`QBARMASSIVEPRIME'}, `SDUMMY', x2?{`QMASSIVEPRIME'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_* d_(colF[idx1], colF[idx3]);
id vx(x1?{`QBAR'}, `PHOPRIME', x2?{`Q'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = charges(x2) * ge * i_* d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHOPRIME', x2?{`L'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = charges(x2) * ge * i_;
id vx(x1?{`QBARMASSIVE'}, `H', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = -gyq(x1) * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`QBARMASSIVE'}, `GLU', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(x1?{`QBARMASSIVE'}, `PHO', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = charges(x2) * ge * i_* d_(colF[idx1], colF[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `H', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = -gyq(x1) * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `GLU', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `PHO', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = charges(x2) * ge * i_* d_(colF[idx1], colF[idx3]);
**************************************************
* END SE vx couplings Feynman rules
**************************************************

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * colf(colA[idx3], colA[idx2], colA[idx1]) * (1/2);
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_;
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_ * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gyq(x1) * i_;
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ghhh * i_;
id vx(x1?{`QBAR'}, `Z', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?, idx4?, idx5?, idx6?, idx7?, idx8?) = -i_ * gw / cw / 2 * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `Z', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?, idx4?, idx5?, idx6?, idx7?, idx8?) = -i_ * gw / cw / 2;

id vx(`H', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = - i_ * d_(colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );
id vx(`H', `GLU', `GLU', `GLU', p4?, p1?, p2?, p3?, idx4?, idx1?, idx2?, idx3?) = i_ * gs * colf(colA[idx1], colA[idx2], colA[idx3]) * ( -gs^2/12/vev/pi^2 );

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

id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * colf(colA[idx1], colA[idx2], colA[idx3]);

id vx(`GLU', `GLU', `GLU', `GLU', 1, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    colf(colA[idx5], colA[idx1], colA[idx2]) * colf(colA[idx3], colA[idx4], colA[idx5]);
id vx(`GLU', `GLU', `GLU', `GLU', 2, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    colf(colA[idx5], colA[idx1], colA[idx3]) * colf(colA[idx2], colA[idx4], colA[idx5]);
id vx(`GLU', `GLU', `GLU', `GLU', 3, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?, idx5?) = -gs^2 * i_ *
    colf(colA[idx5], colA[idx1], colA[idx4]) * colf(colA[idx2], colA[idx3], colA[idx5]);

id vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    colf(colA[idx6], colA[idx1], colA[idx2]) * colf(colA[idx3], colA[idx4], colA[idx6]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    colf(colA[idx6], colA[idx1], colA[idx3]) * colf(colA[idx2], colA[idx4], colA[idx6]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?, idx6?) = -gs^2 * i_ * ( -gs^2/12/vev/pi^2 ) *
    colf(colA[idx6], colA[idx1], colA[idx4]) * colf(colA[idx2], colA[idx3], colA[idx6]);

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id tmp(x?) = x;

B+ T, colf;
.sort:feynman-rules-global;
Keep brackets;

******************
* Color evaluation
******************
repeat id T(coli1?,?a,coli2?)*T(coli2?,?b,coli3?) = T(coli1,?a,?b,coli3); * collect the colour string
id T(coli1?, ?a, coli1?) = colTr(?a);
repeat;
    id colTr = 3;
    id colTr(coli1?) = 0;
    id colTr(?a,coli1?,?b,coli1?,?c) = 1/2*colTr(?a,?c)*colTr(?b) - 1/6*colTr(?a,?b,?c);
    id colf(colj1?,colj2?,colj3?) = -2*i_*colTr(colj1,colj2,colj3) + 2*i_*colTr(colj3,colj2,colj1);
    id colTr(?a,colj1?,?b)*colTr(?c,colj1?,?d) = 1/2*colTr(?a,?d,?c,?b) - 1/6*colTr(?a,?b)*colTr(?c,?d);
endrepeat;
.sort:color;

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
repeat id prop(`GLU', in, p?, idx1?)*prop(`GLU', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
repeat id prop(`Z', in, p?, idx1?)*prop(`Z', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = 1;
id prop(`PHO', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(`Z', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = 1;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;
id prop(x?{`PSI'}, virtual, p?, idx1?, idx2?) = 1;
id prop(x?{`PSI'}, in, p?, idx1?) = 1;
id prop(x?{`PSI'}, out, p?, idx1?) = 1;

**************************************************
* START SE prop Lorentz Feynman rules
**************************************************
* The original spin-sum can be used:
*repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
* or one can include a projector, like below
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = -d_(lorentz[idx1], lorentz[idx2]) + energyselector(lorentz[idx1]) * energyselector(lorentz[idx2]);
* same for all repeat ID below:
repeat id prop(x1?{`QBARMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVE'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x1)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x1?{`QBARMASSIVE'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVEPRIME'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x1)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x1?{`QMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QMASSIVE'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x1)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x1?{`QMASSIVE'}, in, p?, idx1?)*prop(x2?{`QMASSIVEPRIME'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x1)*gamma(dirac[idx1], dirac[idx2]);
id prop(`SDUMMY', virtual, p?, idx1?, idx2?) = 1;
id prop(`SDUMMY', in, p?, idx1?) = 1;
id prop(`SDUMMY', out, p?, idx1?) = 1;
id prop(`PHOPRIME', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`QMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]);
id prop(x?{`QBARMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]);
**************************************************
* END SE prop Lorentz Feynman rules
**************************************************

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

**************************************************
* START SE vx Lorentz Feynman rules
**************************************************
id vx(`PHOPRIME', `PHOPRIME', `SDUMMY', p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(lorentz[idx1], lorentz[idx2]);
id vx(x1?{`QBARMASSIVEPRIME'}, `SDUMMY', x2?{`QMASSIVEPRIME'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`QBAR'}, `PHOPRIME', x2?{`Q'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`LBAR'}, `PHOPRIME', x2?{`L'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBARMASSIVE'}, `H', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`QBARMASSIVE'}, `GLU', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBARMASSIVE'}, `PHO', x2?{`QMASSIVEPRIME'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `H', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `GLU', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBARMASSIVEPRIME'}, `PHO', x2?{`QMASSIVE'}, `SDUMMY', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
**************************************************
* END SE vx Lorentz Feynman rules
**************************************************

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) ;
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = (p3(lorentz[idx2])-p1(lorentz[idx2]));
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(dirac[idx1], dirac[idx3]);
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = 1;
id vx(x1?{`QBAR'}, `Z', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?, idx4?, idx5?, idx6?, idx7?, idx8?) = zVcoupling(x2) * gamma(dirac[idx1], lorentz[idx2], dirac[idx3])
                            - zAcoupling(x2) * gamma(dirac[idx1], lorentz[idx2], dirac[idx4])*gamma5(dirac[idx4], lorentz[idx5], lorentz[idx6], lorentz[idx7], lorentz[idx8], dirac[idx3]);
id vx(x1?{`LBAR'}, `Z', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?, idx4?, idx5?, idx6?, idx7?, idx8?) = zVcoupling(x2) * gamma(dirac[idx1], lorentz[idx2], dirac[idx3])
                            - zAcoupling(x2) * gamma(dirac[idx1], lorentz[idx2], dirac[idx4])*gamma5(dirac[idx4], lorentz[idx5], lorentz[idx6], lorentz[idx7], lorentz[idx8], dirac[idx3]);

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

id D = rat(4-2*ep, 1);
.sort:feynman-rules-vertices-1;

* construct gamma string, drop odd-length gamma traces, pair gamma5s and symmetrize the trace
repeat;
    id once gamma(s1?, ?a, s2?)*gamma5(s2?, ?b, s3?) = gamma5(s1, ?b, s2)*gamma(s2, ?a, s3)*(-1)^nargs_(?a);
    id once gamma5(s1?, ?a, s2?)*gamma5(s2?, ?b, s3?) = gamma(s1, s3);
    id once gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
    id once gamma(s1?,?a,s1?) = gammatrace(?a)*delta_(mod_(nargs_(?a), 2));
    id once gamma(s1?, ?a, s2?)*gamma5(s2?,?b, s1?) = e_(?b)*gammatrace(?a,?b)*delta_(mod_(nargs_(?a),2))/fac_(4);
endrepeat;
contract 0;
id gammatrace = 4;

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

id vx(`GLU', `GLU', `GLU', `GLU', 1, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 2, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`GLU', `GLU', `GLU', `GLU', 3, p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?,idx5?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

id vx(`H', `GLU', `GLU', `GLU', `GLU', 1, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 2, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);
id vx(`H', `GLU', `GLU', `GLU', `GLU', 3, p5?, p1?, p2?, p3?, p4?, idx5?, idx1?, idx2?, idx3?, idx4?,idx6?) =
    d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]);

if (count(vx, 1));
    Print "Unsubstituted vertex: %t";
    exit "Critical error";
endif;

id D = rat(4-2*ep, 1);
.sort:feynman-rules-vertices-2;

Multiply counter(1);
repeat id gammatrace(?a,p?,?b)*counter(i?) = vec(p,lorentzdummy[i])*gammatrace(?a,lorentzdummy[i],?b)*counter(i+1);
id counter(n?) = 1;

id gammatrace(?a) = gammatracetensor(?a);
id vec(p?,mu?) = p(mu);
id D = rat(4-2*ep, 1);

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
    id D = rat(4-2*ep, 1);
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
repeat id cmb(?a)*tder(?c)*forestid(?b,x?) = f(forestid(?b,x*cmb(?a)*tder(?c)))*cmb(?a)*tder(?c);
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

id forestid(x1?,?a,x3?) = nloops(nargs_(?a))*forest(x1)*x3;

.sort:tensorforest-splitoff-2;
Hide F;

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
#do i=0,5
    repeat id subgraph(`i',?e,x1?,fmb1?,?a,fmb2?,?b)*f?{prop,vx}(?c,fmb2?,?d) = subgraph(`i',?e,x1*f(?c,fmb2,?d),fmb1,?a,fmb2,?b);
    repeat id subgraph(`i',?e,x1?,fmb1?,?a)*f?{prop,vx}(?b,fmb1?,?c) = subgraph(`i',?e,x1*f(?b,fmb1,?c),fmb1,?a);
#enddo
repeat id subgraph(?a,p?) = subgraph(?a);
.sort:graphs;

* Create the local UV counterterm
#define uvdiagtotalstart "`extrasymbols_'"
#do i = 1,1
* Split off all subgraphs without dependencies to a separate expression
    id subgraph(x1?, x2?) = uvconf1(x1, x2);
    argtoextrasymbol tonumber,uvconf1,2;
    #redefine uvdiagstart "`extrasymbols_'"
    .sort:uvdiag-1;
    #redefine uvdiagend "`extrasymbols_'"
    #do ext={`uvdiagstart'+1},`uvdiagend'
        L uvdiag`ext' = extrasymbol_(`ext');
    #enddo
    .sort:uvdiag-2;

* fill in results from subdiagrams
    #do ext={`uvdiagtotalstart'+1},`uvdiagstart'
        id uvconf1(`ext') = uvdiag`ext';
    #enddo
    .sort:uv-subgrap-fillin;
    Hide tensorforest1,...,tensorforest`tensorforestcount';

* Apply Feynman rules to the UV subgraph
    id opengammastring(?a) = gamma(?a);
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

* evaluate massive self-energies around p_os=(m,0,0,0) for the 0th order term and
* take derivatives in p^mu but use (p - gamma^0*p_os)^mu as a prefactor for the higher-order terms
    AB+ cmb,diag,fmbtocb;
    .sort:onshell-treatment-1;
    Keep brackets;

    if (count(onshell,1));
        id uvconf2(p?) = 1;

* we use that the cmb momenta that make up the bubble external momentum only appear in that combination
        splitfirstarg onshell;
        id onshell(p1?,-p?vector_,E?) = onshell(-p1,p,-E);
        id onshell(-p?vector_,E?) = onshell(p,-E);

        Multiply xnoos + xos + xneg; * create the OS terms
        if (count(xos,1) || count(xneg,1));
            if (count(xneg,1));
                id onshell(k?,p?,E?) = replace_(p, -E*energyselector - k);
                id onshell(p?,E?) = replace_(p, -E*energyselector);
            else;
                id onshell(k?,p?,E?) = replace_(p, E*energyselector - k);
                id onshell(p?,E?) = replace_(p, E*energyselector);
            endif;
            id uvprop(k?,t1?,p?,m?) = uvprop(k,t1,0,m); * prevent expansion of the propagators
            id xos = 1;
        else;
* wrap the OS mass in a function so that it does not get set to 0 in the mass expansion later on
            id onshell(k?,p?,E?) = replace_(p, t * (p + k - gamma0*tmp(E)*energyselector) - k);
            id onshell(p?,E?) = replace_(p, t * (p - gamma0*tmp(E)*energyselector));
        endif;
    else;
        id uvconf2(p?) = replace_(p, t * p);
    endif;
    .sort:onshell-treatment-2;

    argument uvprop,1,vxs,uvtopo,irtopo,diag,onshell,intuv;
        Multiply replace_(t, 1);
    endargument;

* Taylor expand the numerator and propagators to order tmax
* The orders up to tmax - 1 will only be expanded in the external momenta
* and the tmax order will be expanded in the masses too and every propagator will get an mUV mass
* this way, spurious IR divergences in the external momenta of self-energies will be canceled as well
* note: p carries a t dependence that determines the order, t1 determines the powers of the UV propagator
    id uvprop(k?,t1?,0,m?) = uvprop(k,t1,1,m);
    repeat;
        id once ifnomatch->skiptruncation uvprop(k?,t1?,p?,m?)*t^x1?*tmax^x2? = uvprop(k,t1,1,m) * t^x1*tmax^x2 * theta_(x2-x1) *
            (1 + (-2*p.k - p.p) * t1 + 4*p.k^2 * t1^2) + ALARM * t^3;
        id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
        label skiptruncation;
    endrepeat;

    if (count(ALARM, 1));
        Print "UV Taylor expansion depth exceeded.";
        exit "";
    endif;

* drop the 0th order of the non-OS CT
    if ((count(t,1) == 0) && count(xnoos,1)) Discard;
    id xnoos = 1;

    if (count(t,1) == count(tmax,1));
* drop IR topology and all masses and keep UV topology
        id xirexpand = 0;
        id xuvexpand = 1;
        id m?allmasses = 0;
    else;
        #ifdef `UVTEST'
            Discard;
        #endif
* keep a UV version of the IR topology
        id uvtopo(?a) = irtopo(?a)*xirexpand + uvtopo(?a)*xuvexpand;
        id xirexpand*xuvexpand = 0;
        id xneg*xuvexpand = 0; * only include the UV version once as xos + xneg gives one times the UV CT

* drop the integrated counterterm as they always evaluate to 0
        id intuv(x?) = x;
        if ((count(irtopo,1)) && (count(vxs,1))) Discard;
    endif;

    id tmp(m?) = m;

    id t = 1;
    id tmax = 1;

* match the denominator structure to a diagram
    id uvprop(?a,m?) = uvprop(?a);
    id uvtopo?{uvtopo,irtopo}(x?,?a) = uvtopo(x,1,?a);
    repeat id t1?ts^n2? = tmp(t1,n2);
    repeat id uvprop(k?,t1?,n1?)*uvprop(k?,t1?,n2?) = uvprop(k,t1,n1+n2);
    repeat id uvprop(k?,t1?,n1?)*uvtopo?{uvtopo,irtopo}(x?,x1?,?a)*tmp(t1?,n2?) = uvprop(k,n1 + n2)*uvtopo(x,x1*t1^n2,?a);
    id uvprop(k?,t1?ts,n?) = uvprop(k, n);

*   Unwrap the integrated UV expression
    id intuv(x?) = x;

* select what to keep for the local UV vs integrated UV CT
    if (count(vxs, 1) == 0);
        id xnomsbar = 0;
        id xlct = 1;
        id uvprop(?a) = xlct;
    else;
* if a subdiagram of the graph has a local CT, drop it and replace it by its untruncated integrated version
        id xlct = 0;
        id xnomsbar = 1;
        id uvtopo(?a) = 1;
        id irtopo(?a) = 1;
    endif;

* select the OS topology with the negative mass
    if (count(xneg,1));
        id irtopo(x1?,x2?,?a) = irtopo(x1,x2,xneg,?a);
        id xneg = 1;
    endif;

    #call uvmap()

* insert the gamma^0_ij from the OS CT
    id opengammastring(s1?,?a,s2?)*gamma0 = opengammastring(s1,energyselector,?a,s2);

    if (count(uvtopo, 1, irtopo, 1, tmp, 1));
        Print "Unsubstituted UV topology: %t";
        exit "Critical error";
    endif;

* compute the integrated UV counterterm
    Multiply counter(1);
    repeat id k1?.k2?*counter(n?) = vec(k1,n)*vec(k2,n)*counter(n + 1);
    id opengammastring(?a) = gamma(?a);
    repeat id gamma(s1?,?a,k1?,?b,s2?)*counter(n?) = gamma(s1,?a,n,?b,s2)*vec(k1,n)*counter(n + 1);

    id e_(mu1?, mu2?, mu3?, mu4?) = levicivita(mu1, mu2, mu3, mu4); * the reverse substitution is delayed until after the call to IntegrateUV()
    repeat id levicivita(?a, k1?, ?b)*counter(n?) = levicivita(?a,n,?b)*vec(k1,n)*counter(n+1);

    repeat id k1?(mu?)*counter(n?) = vec(k1,n)*vec1(mu,n)*counter(n + 1);

* convert every k^0 into k.p0select, where p0select is effectively (1,0,0,0)
    repeat id penergy(k1?)*counter(n?) = vec(k1,n)*vec(p0select,n)*counter(n + 1);
    id counter(x?) = 1;

    Multiply replace_(vec, vec1); * consider all vectors as external
    repeat id uvprop(k1?,n?)*vec1(k1?,n1?) = uvprop(k1,n)*vec(k1,n1);

    #call TensorReduce()

* contract all metrics
    repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);
    id g(n1?,n1?) = rat(4-2*ep,1);
    repeat id vec1(p1?,n1?)*g(n1?,n2?) = vec1(p1,n2);

    .sort:tensor-projection-loop;

    if (count(uvprop,1));
        id k1?.k2? = g(k1,k2);

* divide by the normalizing factor of the denominator that is added to the topology
* this is always 1/(k^2 - m_UV^2)^3 = -i / (4 pi)^2 * 1/2 * 1/mUV^2
* note that mUV is ignored in the Taylor expansion of masses
* determine the number of loops with E - V + 1
        Multiply replace_(vxs,vx);
        id vx(?a) = x1*x2^nargs_(?a)*vxs(?a);
        id x1^n1?*x2^n2? = (i_ * (4 * pi)^2 * 2 * mUV^2)^(n2/2 - n1 + 1);
    endif;

    #call IntegrateUV()

    id D = rat(4-2*ep,1);
    id ep^n1? = rat(ep^n1,1);
    .sort:ibp-reduction;

    id g(k1?,k2?) = k1.k2; * k1,k2 should be external only
    id vec1(p?,n1?)*vec1(mu?,n1?) = p(mu);
    id vec1(k1?,n?)*vec1(k2?,n?) = k1.k2;
    id k1?.p0select = penergy(k1);

    repeat id gamma(s1?,?a,n?,?b,s2?) = gamma(s1,?a,lorentzdummy[n],?b,s2);
    repeat id levicivita(?a,n?,?b) = levicivita(?a,lorentzdummy[n],?b);
    id g(n1?,n2?) = d_(lorentzdummy[n1], lorentzdummy[n2]);
    id vec1(k1?,n?) = k1(lorentzdummy[n]);
    id vec1(mu?,n1?) = d_(mu, lorentzdummy[n1]);
    id gamma(?a) = opengammastring(?a);
    id levicivita(?a) = e_(?a);

* Simplify all open gamma strings so that the Ds are contracted out
    #call Gstring(opengammastring,0)
    id energyselector.energyselector = 1;

    if (match(opengammastring(s1?,?a,mu?lorentzdummy,?b,s2?)) || count(vec1,1));
        Print "Warning: dummy index left in gamma string after projecting: %t";
* this happens when there are two gamma strings in a subgraph
* TODO: is there a risk for dummy index collisions?
*        exit "Critical error";
    endif;

* Substitute the masters and expand in ep
    #call SubstituteMasters()

* subtract MS-bar contributions (the poles)
    #ifndef `NOMSBARSUBTRACTION'
* store an unsubtracted version that can be used as the integrated version of the local counterterm
        if (count(xlct, 1) == 0) Multiply xmsbar - xnomsbar;
        if (count(xmsbar, 1));
            argument rat;
                id ep^n? = ep^n*theta_(n);
            endargument;
        endif;
        id xmsbar = 1;
    #endif

    .sort:uv-subgraph-done;
    UnHide tensorforest1,...,tensorforest`tensorforestcount';
    Hide uvdiag{`uvdiagstart'+1},...,uvdiag`uvdiagend';

* in the next loop the UV subgraph will be added to the tensorforest or to a UV graph that embeds it
    if (count(uvconf1, 1)) redefine i "0";

* fill in the unsubstituted subgraph into the supergraph
    repeat id subgraph(x1?,?a,n?,?b,x2?)*uvconf1(n?,x3?) = subgraph(x1,?a,?b,x2*uvconf1(x3));
    id uvconf1(n?,x?) = uvconf1(x);
    .sort:uv3;
#enddo

id xnomsbar = 0;
id xlct = 1;
id xuvexpand = 0;
id xirexpand = 1;
.sort:local-uv-done;
Drop uvdiag`uvdiagtotalstart',...,uvdiag`uvdiagend';
delete extrasymbols>`uvdiagtotalstart';

if (count(subgraph, 1));
    Print "Unsubstituted UV subgraph: %t";
    exit "Critical error";
endif;

* Apply Feynman rules to remaining graph
id opengammastring(?a) = gamma(?a);
#call FeynmanRulesMomentum()

* Linearize the gamma string
repeat id once gamma(s1?,?a,p?!vector_,?b,s2?) = p(mudummy)*gamma(s1,?a,mudummy,?b,s2);
id gamma(?a) = opengammastring(?a);

* Simplify all open gamma strings
#call Gstring(opengammastring,0)

* split off the energy part: energyselector=(1,0,0,0,..), fmbs = spatial part
* fmbs1.fmbs2 will later get a minus sign to honour the Minkowksi metric
AB+ cmb,energy,diag,fmbtocb;
.sort:energy-splitoff-1;
Keep brackets;
#do i=1,`NFINALMOMENTA'
    Multiply replace_(fmb`i', fmb`i'*xorig + xsplit*energy(fmb`i')*energyselector - xsplit*fmbs`i');
    id xorig * xsplit = 0;
#enddo
id energyselector.energyselector = 1;
id energyselector.p?spatialparts = 0;

id xorig^n?{>0} = xorig;
id xsplit^n?{>0} = xsplit;
.sort:energy-splitoff-2;

* extract all scalar forests into new expressions, except for when the forest
* does not have any loops or has all the loops
if ((match(nloops(0)) == 0) && (match(nloops({`NFINALMOMENTA'-1})) == 0));
    Multiply f(1);
    repeat id f(x?)*f1?{cmb,diag,energy,tder}(?a) = f(x*f1(?a));
endif;
id tder(?a) = 1;
id nloops(x?) = 1;

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
id xorig * xsplit = 0;
if (count(xorig,1));
    Multiply replace_(<fmbs1,fmb1>,...,<fmbs40,fmb40>);
    argument fmbtocb;
        Multiply replace_(<ps1,p1>,...,<ps40,p40>,<cs1,c1>,...,<cs40,c40>);
    endargument;
endif;

id fmbtocb(?a) = replace_(?a);
.sort:tensorforest-2;
Hide forest1,...,forest`forestcount';
UnHide F;

id xorig * xsplit = 0;
id xorig^n?{>0} = xorig;
id xsplit^n?{>0} = xsplit;

B+ forestid;
.sort:tensorforest-3;
Drop tensorforest1,...,tensorforest`tensorforestcount';
Keep brackets;

* fill in the tensor forest into F and fill in the Feynman rules
#do i=1,`tensorforestcount'
    id forestid(`i') = tensorforest`i';
#enddo
id f(x?) = forestid(x-`foreststart');
repeat id cmb(?a)*cmb(?a) = cmb(?a);

* here we can already truncate, as there are no more poles
argument rat;
    id ep^{`SELECTEDEPSILONORDER' + 1} = 0;
endargument;
.sort:tensorforest-4;

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
PolyRatFun;

argument rat;
    if (count(ep, 1) != `SELECTEDEPSILONORDER') Discard;
endargument;
id rat(x1?) = x1;
id ep^n? = 1;
#if `UVRENORMFINITEPOWERTODISCARD' > 0
    if (count(UVRenormFINITE, 1) >= `UVRENORMFINITEPOWERTODISCARD') Discard; * Discard UVRenormFinite pieces up to some order
#elseif  `UVRENORMFINITEPOWERTODISCARD' < 0
    if (count(UVRenormFINITE, 1) < -`UVRENORMFINITEPOWERTODISCARD') Discard; * Keep UVRenormFinitiePieces up to some order
#endif
id UVRenormFINITE^n? = 1;

B+ forestid;
.sort:rat-truncate;
Keep brackets;

#$activeforestcount = `forestcount';
* Drop unused forests
#do i=1,`forestcount'
    #define HASFOREST`i' "0"
    if (match(forestid(`i'))) redefine HASFOREST`i' "1";
#enddo
.sort:forest-flatten-3;
UnHide forest1,...,forest`forestcount';
.sort:forest-flatten-4;

#do i=1,`forestcount'
    #if `HASFOREST`i'' == 0
        #$activeforestcount = $activeforestcount - 1;
        L forest`i' = 0;
    #endif
#enddo

.sort:forest-flatten-5;
Hide forest1,...,forest`forestcount';

.sort

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
.sort:cmb-transform;
UnHide forest1,...,forest`forestcount';

id rat(x1?,x2?) = x1/x2;
id rat(x?) = x;

* collect all the energies in the diagram
repeat id diag(?a,p1?,p2?,?b) = diag(?a,p1,0,p2,?b);
id diag(?a,p1?) = diag(?a,p1,0);

repeat id energy(k?)*diag(?a,k?,n?,?b) = diag(?a,k,n+1,?b);

if (count(energy,1));
    Print "Energy left: %t";
    exit "Critical error";
endif;

if ((count(cmb,1) == 0) && count(diag,1));
    Print "CMB missing: %t";
    exit "Critical error";
endif;
B xsplit,xorig;
.sort:energy-collect;

*********************************************
* Construction of the integrand
*********************************************

#if (`INTEGRAND' == "LTD") || (`INTEGRAND' == "both")
    #include- ltdtable_`SGID'.h
    Hide forest1,...,forest`forestcount',F;

* TODO: unify call and make decision at runtime
    #do i=1,`forestcount'
        L forestltd`i' = forest`i';
    #enddo
    L Fltd = F;

* only apply compressed numerator to case without nested UV CT
    if (count(forestid, 1) || (match(diag(x1?,x2?,c1?,?a)*diag(x3?,x4?,c2?,?b))));
        id xsplit = 1;
        id xorig = 0;
    else;
        id xsplit = 0;
    endif;

    id xorig^n?{>0} = xorig;
    B xorig;

    .sort
    Hide Fltd,forestltd1,...,forestltd`forestcount';
    L rest = Fltd[1];
    L diagorig = Fltd[xorig];

* split compressed numerator per diag as well, as they may have different energy dependence
    B conf,tder,cmb,diag;
    .sort;
    Hide rest;
    collect f; * let's hope this works for large expressions!

    repeat id f(x?)*diag(?a) = f(x*diag(?a));

    argument f;
        id diag(?a) = diag(?a,1);
        id p1?.p2? = dot(p1,p2,p1.p2);
        argument dot,3;
            Multiply replace_(<p1,ps1>,...,<p40,ps40>,<c1,cs1>,...,<c40,cs40>);
        endargument;
    endargument;
    id f(x?) = diag(x,-1); * tag compressed numerator
    .sort
    UnHide Fltd,forestltd1,...,forestltd`forestcount';
    .sort
    Drop rest, diagorig;
    L Fltd = rest + diagorig;
    .sort:ltd-forest-copy;

    repeat id cmb(?a)*tder(?c)*diag(?b) = f(diag(?b,cmb(?a)*tder(?c)))*cmb(?a)*tder(?c);
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

    id diag(x?,-1,y?) = x*y*xorig;
    id diag(x?) = diag(x-`diagstart'+`forestcount');
    id diag(xcut?,x2?,?a,x3?) = diag(?a)*ltdtopo(x2)*x3*conf(-1,xcut,-1);
    id cmb(?a) = replace_(?a);

    repeat id conf(?a)*conf(?a) = conf(?a);
    repeat id constants(x?)*constants(x1?) = constants(x*x1);
    repeat id allenergies(?a)*allenergies(?b) = allenergies(?a,?b);
    .sort:ltd-splitoff;
    Hide forestltd1,...,forestltd`forestcount',Fltd;

    Multiply norm(1)*num();
    id prop(?a) = prop(?a,1);

    B+ dot,onshell,diag,xorig;
    .sort:fmb-conv;
    Keep brackets;
    #do i=1,`NFINALMOMENTA'
* NOTE: taking derivatives is only possible if `fmb1` is a vector_ and not a sum of them
        if (count(xorig,1));
            id diag(fmb1?,n?,?a) = replace_(fmb1,ltd`i'*energyselector)*diag(?a);
        else;
            id diag(fmb1?,n?,?a) = ltd`i'^n*diag(?a); * order should be consistent
        endif;
        argument dot 1,2;
            id energyselector = 1;
        endargument;
        id diag = 1;
    #enddo

    argument dot 1,2;
        id p? = penergy(p);
    endargument;

    argument dot,3;
        Multiply replace_(<ps1,p1>,...,<ps40,p40>,<cs1,c1>,...,<cs40,c40>);
    endargument;
    id dot(?a,p1?.p2?) = dot(?a,spatial(p1,p2));

    #do i=1,`NFINALMOMENTA'
        symmetrize dot 1,2;
        id dot(x1?,x2?,x3?) = dot(x1,xsplit,x2,x3);
        splitarg (ltd`i') dot;
        id dot(x?,xsplit,?a) = dot(0,x,xsplit,?a);
        id dot(x1?,x2?,xsplit,x3?,x4?) = dot(x1,x2,xsplit,0,x3,x4);
        repeat id dot(?a,x?,-ltd`i',?b,x4?) = -dot(?a,-x,ltd`i',?b,-x4);

        B+ prop, ltd`i', norm, num, diag, onshell, xorig, dot;

        .sort:loop-{`i'};
        Keep brackets;

        repeat id prop(?a,n1?)*prop(?a,n2?) = prop(?a,n1+n2);

* Identify poles in the loop momentum energy
        splitarg (ltd`i') prop;
        id prop(ltd`i',E?,n?) = prop0(ltd`i',0,E,n);
        id prop(p0?,ltd?,E?,n?) = prop0(ltd,p0,E,n);

        factarg (-1) prop0 1;
        id prop0(ltd`i',y?{-1,1},E?,n?) = prop0(ltd`i',y,0,E,n);
        id all prop0(ltd?,y?{-1,1},p0?,E?,n?)*norm(n0?) =
            sum_(m,0,n-1,
                (-1)^(n-m-1)*fac_(2*n-m-2)/fac_(n-m-1)/fac_(m)*
                der(ltd`i',m)*norm(n0*(2*E)^(2*n-m-2))
                )/fac_(n-1)*
            a(y*p0,0,E)*a(y*p0,E,0);

* Apply derivative
        repeat;
            repeat id der(ltd`i',n1?{>0})*prop0(ltd`i',y?{-1,1},p0?,E?,n2?) =
                        +der(ltd`i',n1)*prop(ltd`i',y,p0,E,n2)
                        -der(ltd`i',n1-1,0)*(2*n2*(ltd`i'+y*p0))*prop0(ltd`i',y,p0,E,n2+1);
            id der(ltd`i',n1?{>0})*ltd`i'^n3? = der(ltd`i',n1)*ltd`i'^n3 + n3*der(ltd`i',n1-1,0)*ltd`i'^(n3-1);
            id all der(ltd`i',n1?{>0})*dot(?a,x1?,ltd`i',?b,x2?) = der(ltd`i',n1-1,0)*dot(?a,0,1,?b,0);
            id der(ltd`i',n?{>0}) = 0;
            id der(ltd`i',n?,0) = der(ltd`i',n);
            id prop(ltd`i',y?{-1,1},p0?,E?,n2?) = prop0(ltd`i',y,p0,E,n2);
        endrepeat;
        id der(ltd`i',0) = 1;
        id prop0(ltd?,y?{-1,1},?a) = prop(y*ltd,?a);

* check for errors
        if (count(der,1));
            print "Remaining diff operators: %t";
            exit "Critical ERROR";
        endif;

* Split the poles by their position in the complex plane
        splitarg a;
        repeat id a(?y1,+E?energyset,?y2,E1?,E2?) = a(?y1,?y2,E1+E,E2);
        repeat id a(?y1,-E?energyset,?y2,E1?,E2?) = a(?y1,?y2,E1,E2+E);
        id a(E1?,E2?) = a(0,E1,E2);

* check for errors
        id a(?y,0,0) = tmp(a(?y,0,0));
        if (count(tmp,1));
            print "[ltd`i'] Cannot find energy signature %t";
            exit "Critical ERROR";
        endif;
        transform, a, prop, addargs(1,last-2);

* take residue
        id a(p0?,0,E?)*a(p0?,E1?,E2?)*norm(n?)*num(?y) = replace_(ltd`i',E-p0)*norm(n*(E+E1-E2))*num(?y,E-p0);
        if ( count(a,1) ) discard;

* fuse dots
        repeat id dot(x1?,x2?,?a,xsplit,?b,x3?) = dot(x1+x2,?a,xsplit,?b,x3);
        repeat id dot(x1?,xsplit,x2?,x3?,?b,x4?) = dot(x1,xsplit,x2+x3,?b,x4);
        id dot(x1?,xsplit,x2?,x3?) = dot(x1,x2,x3);
    #enddo
    id prop(?a,n?) = prop(?a)^n;
    id prop(p0?,E?) = den(p0^2-E^2);
    id num(?a) = 1;
    id norm(n?) = 1/n;
    id diag = 1;
    id xorig = 1;

    .sort:ltd;
    argtoextrasymbol tonumber, den;
    #define denstart "`extrasymbols_'"
    .sort:invd-creation;
    #do i = {`denstart'+1},`extrasymbols_'
        id den(`i') = invd{`i'-`denstart'-1};
    #enddo
    #if `denstart' < `extrasymbols_'
        Multiply ellipsoids(<invd0,extrasymbol_({`denstart'+1})>,...,<invd{`extrasymbols_'-`denstart'-1},extrasymbol_(`extrasymbols_')>);
    #else
        Multiply ellipsoids;
    #endif

    B allenergies, ellipsoids,conf,tder,constants,pi;
    .sort:ltd-num-0;
    UnHide forestltd1,...,forestltd`forestcount',Fltd;
    delete extrasymbols>`denstart';

* Some expressions may have evaluated to 0, so remove the call
    #do i = 1,`diagcount'
        #if (termsin(diag`i') == 0)
            id diag({`i'+`forestcount'}) = 0;
        #endif
    #enddo

    .sort:ltd-num-1;
    Drop diag1,...,diag`diagcount',forestltd1,...,forestltd`forestcount',Fltd;

* now add all LTD structures as a special conf
    #if `diagcount' > 0
        #if `forestcount' > 0
            L FINTEGRANDLTD = Fltd + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})> + <forestltd1*conf(-1,-1)>+...+<forestltd`forestcount'*conf(-`forestcount',-1)>;
        #else
            L FINTEGRANDLTD = Fltd + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})>;
        #endif
    #endif

* Some forests may be 0 due to an LTD expression being 0 inside it, so remove the call
    #do i = 1,`forestcount'
        #if (termsin(forestltd`i') == 0)
            id forestid({`i'}) = 0;
        #endif
    #enddo

    id conf(x?)*conf(x1?{<0},?a) = conf(x,?a);

    .sort:integrand-ltd;
    #if (`INTEGRAND' == "both")
        Hide FINTEGRANDLTD;
        delete extrasymbols>`diagstart';
    #endif
#endif

#if (`INTEGRAND' == "PF") || (`INTEGRAND' == "both")
    .sort:pf-start-1;
    UnHide forest1,...,forest`forestcount',F;
    .sort:pf-start-2;

    #include- pftable_`SGID'.h
    .sort:load-pf;
    id xorig = 0;
    id xsplit = 1;

* TODO: here we create a unique diagram per cut as we include the cmb (and conf) into the key
* it is possible that diagrams are the same after the cmb is applied, for example for UV topologies
    repeat id cmb(?a)*tder(?c)*diag(?b) = f(diag(?b,cmb(?a)*tder(?c)))*cmb(?a)*tder(?c);
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
    .sort:pf-splitoff-1;

    id diag(xcut?,x2?,?a,x3?) = diag(?a)*pftopo(x2)*x3*conf(-1,xcut,-1);
    id cmb(?a) = replace_(?a);

    .sort:pf-splitoff-2;
    Hide forest1,...,forest`forestcount',F;

* apply the numerator procedure
    #do i = 0,{`NFINALMOMENTA'-1}
        id diag(p?,n?,?a) = diag(?a)*x1^n*f(p);
        repeat id energy(p?)*f(p?) = x1*f(p); * energies could come from previous iteration
        id f(p?) = 1;
        id once num(ncmd(?z),?x) =ncmd(?z)*num(?x);
        B+ ncmd, x1;
        .sort:collect-ncmd-`i';
        keep brackets;

        $aa = 1;
        id x1^np?*ncmd(?z,x?) = tmp(x1^np*ncmd(?z,x));
        id tmp(x?$aa) = 1;

        inside $aa;
            id x1^np?*ncmd(?z,x?) = sum_(n,0, np - nargs_(?z), a(np,n,?z)*tmp(x,n));
            repeat id a(np?,n?,?z,x?) = -sum_(n1, n+1,np-nargs_(?z), a(np,n1,?z)*tmp(x,(n1-n-1)));
            id a(np?,n?) = delta_(np,n);

            id tmp(x?,n?) = x^n;
        endinside;

        Multiply $aa;
        ModuleOption local $aa;
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
    UnHide forest1,...,forest`forestcount',F;

* Some PF expressions may have evaluated to 0, so remove the call
    #do i = 1,`diagcount'
        #if (termsin(diag`i') == 0)
            id diag({`i'+`forestcount'}) = 0;
        #endif
    #enddo

    .sort:pf-num-2;
    Drop diag1,...,diag`diagcount',forest1,...,forest`forestcount',F;

* now add all PF structures as a special conf
    #if `diagcount' > 0
        #if `forestcount' > 0
            L FINTEGRANDPF = F + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})> + <forest1*conf(-1,-1)>+...+<forest`forestcount'*conf(-`forestcount',-1)>;
        #else
            L FINTEGRANDPF = F + <diag1*conf(-{1+`forestcount'})>+...+<diag`diagcount'*conf(-{`diagcount'+`forestcount'})>;
        #endif
    #endif

    id conf(x?)*conf(x1?{<0},?a) = conf(x,?a);

* Some forests may be 0 due to a PF expression being 0 inside it, so remove the call
    #do i = 1,`forestcount'
        #if (termsin(forest`i') == 0)
            id forestid({`i'}) = 0;
        #endif
    #enddo

    .sort:pf-integrand-collect;
    #if (`INTEGRAND' == "both")
        UnHide FINTEGRANDLTD;
    #endif
#endif

id replace(?a) = replace_(?a);
id energy(p?) = penergy(p);
id energies(p?) = penergy(p);

B+ penergy,spatial,energies,allenergies,ellipsoids,constants,dot;
.sort:func-prep;
Keep brackets;

argument ellipsoids, dot;
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
    argument energies, ellipsoids, constants, dot;
        id penergy(p`i') = lm`$OFFSET';
    endargument;
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        argument energies, ellipsoids, constants, dot;
            id p`i'.p`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`i', p`j') = lm`$OFFSET';
        argument energies, dot;
            id spatial(p`i', p`j') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXK'
    id penergy(c`i') = lm`$OFFSET';
    argument energies, ellipsoids, constants, dot;
        id penergy(c`i') = lm`$OFFSET';
    endargument;
    #$OFFSET = $OFFSET + 1;
    #do j=1,`$MAXP'
        argument energies, ellipsoids, constants, dot;
            id c`i'.p`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`j', c`i') = lm`$OFFSET';
        argument energies, dot;
            id spatial(p`j', c`i') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=`i',`$MAXK'
        argument energies, ellipsoids, constants, dot;
            id c`i'.c`j' = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(c`i', c`j') = lm`$OFFSET';
        argument energies, dot;
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

        B+ conf,tder;
        .sort:conf-0;
        Keep brackets;

        id conf(?a) = conf(conf(?a));
        id conf(x?)*tder(?a) = conf(x*tder(?a));
        argtoextrasymbol tonumber,conf,1;
        #redefine oldextrasymbols "`extrasymbols_'"
        B+ conf;
        ModuleOption noparallel; * make sure the graph ordering stays intact
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
            ModuleOption noparallel;
            .sort:conf-`ext'-0;


            #if (`$isdenominator' == 1)
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#CONSTANTS\n%$\n#CONSTANTS", $constants
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#ENERGIES\n%$\n#ENERGIES", $energies
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "#ELLIPSOIDS\n%$\n#ELLIPSOIDS",$ellipsoids
            #endif

            Format C;
            #if (`$activeforestcount' < `MAXVARSFOROPTIM')
                #if `OPTIMLVL' > 1
                    Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
                #else
                    Format O1,stats=on;
                #endif

                #Optimize FF`ext'
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "%O"
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "\n\treturn %+E;",FF`ext'
            #else
                .sort:strip-content-0;
                L content = content_(FF`ext');
                ToPolynomial;
                .sort:strip-content-1;
                L FF`ext' = FF`ext' / content;
                FromPolynomial;
                .sort:strip-content-2;
                #write<out_integrand_`INTEGRANDTYPE'_`SGID'.proto_c> "\n//CMODE\n\treturn %E*(%+E);\n",content,FF`ext'
                Drop content;
            #endif

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
