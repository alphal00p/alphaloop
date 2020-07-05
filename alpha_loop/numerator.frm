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

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, vector,g(s),delta(s),T, counter,color, prop;
CF f, vx, vec, vec1;
CF subs, configurations, conf, cmb, der, energy, spatial(s);
CF subgraph, uvconf, uvconf1, uvprop, uv, conjugate;
S integratedctflag, mUV, logmUV, mi1L1, alarmMi1L1;
CF integratedct, rat, num, den;
Set ts: t0,...,t20;
CT penergy;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  i, m, n, ALARM;

#include- diacolor.h
Set colF: cOli1,...,cOli40;
Set colA: cOlj1,...,cOlj40;
Set colAdum: cOljj1,...,cOljj40;

Polyratfun rat;

*--#] setup :

* Load the diagrams
#include- input_`SGID'.h

*--#[ feynman-rules :

************************************************
* Substitute the Feynman rules for the numerator
************************************************

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pdgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);
id conjugate(x?) = x;
* do the spin sum external particles
repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2])*(gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]));
repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = d_(colF[idx1], colF[idx2])*(gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]));

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]) * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`L'}, virtual, p?, idx1?, idx2?) = i_ * (gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1]));
id prop(x?{`LBAR'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2]));
id prop(x?{`Q'}, virtual, p?, idx1?, idx2?) = i_ * (gamma(dirac[idx2], p, dirac[idx1]) + masses(x) * gamma(dirac[idx2], dirac[idx1])) * d_(colF[idx2], colF[idx1]);
id prop(x?{`QBAR'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2])) * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;
id prop(`H', in, p?, idx1?) = 1;
id prop(`H', out, p?, idx1?) = 1;

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * cOlf(colA[idx3], colA[idx2], colA[idx1]) * (1/2) * (p3(lorentz[idx2])-p1(lorentz[idx2]));
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `PHO', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gy * i_ * d_(dirac[idx1], dirac[idx3]) * d_(colF[idx1], colF[idx3]);
id vx(x1?{`LBAR'}, `H', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gy * i_ * d_(dirac[idx1], dirac[idx3]);
id vx(`H', `H', `H', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ghhh * i_;

id D^n? = rat(D^n, 1);
.sort:feynman-rules-vertices-1;

* TODO: use momentum conservation to reduce the number of different terms
id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * gs * cOlf(colA[idx1], colA[idx2], colA[idx3]) *(
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2]));

* For the quartic gluon vertex we need an extra dummy index
Multiply counter(1);
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(i?) = - counter(i + 1) * gs^2 * i_ *(
    + cOlf(colAdum[i], colA[idx1], colA[idx2]) * cOlf(colA[idx3], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]))
    + cOlf(colAdum[i], colA[idx1], colA[idx3]) * cOlf(colA[idx2], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]))
    + cOlf(colAdum[i], colA[idx1], colA[idx4]) * cOlf(colA[idx2], colA[idx3], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]))
);
id counter(x?) = 1;

if (count(vx, 1, prop, 1));
    Print "Unsubstituted propagator or vertex: %t";
    exit "Critical error";
endif;

id D^n? = rat(D^n, 1);
.sort:feynman-rules-vertices-2;

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

* set the SU(3) values
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

if (count(gamma, 1));
    Print "Unsubstituted gamma string: %t";
    exit "Critical error";
endif;

id D^n? = rat(D^n, 1);
.sort:gamma-traces;

id color(x?) = x;

* Set all external momenta on-shell
#do i=1,10
    id p`i'.p`i' = 0;
#enddo

id pzero = 0; * Substitute the 0-momentum by 0
.sort:feynman-rules-final;
*--#] feynman-rules :

* If the expression is empty (due to color), we still write a file
#if ( termsin(F) == 0 )
    #write<out_`SGID'.proto_c> "#0 due to color\n"
#endif

*************************************************
* Process different configurations (bubbles, etc)
*************************************************

* process all the configurations
id configurations(x?) = x;

* multiply the numerator contribution of derivatives
id conf(?a,p?) = conf(?a) * penergy(p);
id conf(?a,x?) = conf(?a) * x; * note the type difference

* Taylor expand in pbubble_i^0
* FIXME: if the number of terms is 0 after differentiating, this configuration will not be created
repeat;
    id once der(p?) = der*replace_(p, p + pzero * x * penergy(p));
    if ((count(der, 1)) && (count(x, 1) != 1)) Discard;
    id x = 1;
    id der = 1;
    id pzero.p? = penergy(p);

    argument subs;
        id x = 0; * undo the replacement in subs
    endargument;
endrepeat;

id subs(p1?,p2?) = replace_(p1, p2);

B+ uv;
.sort:bubble-treatment;
Keep brackets;

argument uv;
   repeat;
* all subgraphs without dependencies can be treated at the same time
* multiply each graph with -1 to correctly subtract it
        id subgraph(x1?, x2?) = -uvconf1(x1, x2);
        argument uvconf1,2;
            id uvconf(x1?,x2?,?a,x3?) = uvconf(x1) * tmax^x2 * uvconf1(?a) * x3;
            chainout uvconf1;
            repeat id uvconf1(p?)*uvconf1(p?) = uvconf1(p);
            id uvconf1(p?) = replace_(p, t * p);

            argument uvprop,1;
                id t = 1; * it could be that the LTD momentum also makes an appearance as an external momentum
            endargument;

* Taylor expand to the right depth
            repeat;
                id once ifnomatch->skiptruncation uvprop(k?,t1?,p?)*t^x1?*tmax^x2? = uvprop(k,t1,1) * t^x1*tmax^x2 * theta_(x2-x1) *
                    (1 - 2 * k.p * t1 - p.p * t1 + 4*p.k^2 * t1^2 + 4*p^2*p.k * t1^2 - 8 * p.k^3 * t1^3 + ALARM * t^4);
                id t^x1?*tmax^x2? = t^x1*tmax^x2 * theta_(x2-x1);
                label skiptruncation;
            endrepeat;

* select the right denominator structure
            repeat id uvprop(k?,t1?,n1?)*uvprop(k?,t1?,n2?) = uvprop(k,t1,n1+n2);
            id uvprop(k?,t1?,n1?)*t1?^n2? = uvprop(k,n1 + n2)*t1^n2;
            id uvprop(k?,t1?ts,n?) = uvprop(k, n);
            if (count(integratedctflag, 1) == 0) id uvprop(?a) = 1;
            id uvconf(x?) = 1/x;
            id t?ts^n? = 0;
            id t = 1;
            id tmax = 1;
        endargument;

* now fill in the subgraph evaluation into the supergraph
        repeat id subgraph(x1?,?a,n?,?b,uvconf(?c,x2?))*uvconf1(n?,x3?) = subgraph(x1,?a,?b,uvconf(?c,x2*x3));
    endrepeat;
   
    id once uvconf1(x?,x1?) = x1;
    if (count(subgraph, 1));
        Print "Unsubstituted UV subgraph: %t";
        exit "Critical error";
    endif;
endargument;
id uv(x?) = x;
.sort:uv-treatment;

* compute the integrated UV counterterm
* FIXME: we are not given the complete denominator: denominators that have no external momentum dependence are not in the subgraph
id uvprop(k?,n?) = uvprop(k, n + 1);
if (count(integratedctflag, 1) > 0);
    id integratedctflag = -1; * we add back the counterterm
    Multiply counter(mu1,...,mu20);
    repeat id k1?.k2?*counter(mu?,?a) = vec(k1,mu)*vec(k2,mu)*counter(?a);
    id counter(?a) = 1;
    repeat id uvprop(k1?,n?)*vec(k1?,mu?) = uvprop(k1,n)*vec1(k1,mu);
    chainin vec1;

* tensor reduce the vacuum bubble
    id vec1(k1?,mu1?) = 0;
    id vec1(k1?,mu1?,k1?,mu2?) = k1.k1 * rat(1, D) * d_(mu1,mu2);
    id vec1(k1?,mu1?,k1?,mu2?,k1?,mu3?) = 0;
    id vec1(k1?,mu1?,k1?,mu2?,k1?,mu3?,k1?,mu4?) = k1.k1^2 * rat(1, D * (2+D)) * (d_(mu1,mu2) * d_(mu3,mu4)
        + d_(mu1,mu3) * d_(mu2,mu4) + d_(mu1,mu4) * d_(mu2,mu3));
    id vec1(k1?,mu1?,k1?,mu2?,k1?,mu3?,k1?,mu4?,k1?,mu5?) = 0;
    id vec1(k1?,mu1?,k1?,mu2?,k1?,mu3?,k1?,mu4?,k1?,mu5?,k1?,mu6?) = k1.k1^3 * rat(1, D * (2+D) * (4+D)) * (
        + d_(mu1,mu2)*d_(mu3,mu4)*d_(mu5,mu6)
        + d_(mu1,mu2)*d_(mu3,mu5)*d_(mu4,mu6)
        + d_(mu1,mu2)*d_(mu3,mu6)*d_(mu4,mu5)
        + d_(mu1,mu3)*d_(mu2,mu4)*d_(mu5,mu6)
        + d_(mu1,mu3)*d_(mu2,mu5)*d_(mu4,mu6)
        + d_(mu1,mu3)*d_(mu2,mu6)*d_(mu4,mu5)
        + d_(mu1,mu4)*d_(mu2,mu3)*d_(mu5,mu6)
        + d_(mu1,mu4)*d_(mu2,mu5)*d_(mu3,mu6)
        + d_(mu1,mu4)*d_(mu2,mu6)*d_(mu3,mu5)
        + d_(mu1,mu5)*d_(mu2,mu3)*d_(mu4,mu6)
        + d_(mu1,mu5)*d_(mu2,mu4)*d_(mu3,mu6)
        + d_(mu1,mu5)*d_(mu2,mu6)*d_(mu3,mu4)
        + d_(mu1,mu6)*d_(mu2,mu3)*d_(mu4,mu5)
        + d_(mu1,mu6)*d_(mu2,mu4)*d_(mu3,mu5)
        + d_(mu1,mu6)*d_(mu2,mu5)*d_(mu3,mu4));

    if (count(vec1, 1) > 0);
        Print "Tensor reduction table insufficient: %t";
        exit "Critical error";
    endif;

* Set all external momenta on-shell
    #do i=1,10
        id p`i'.p`i' = 0;
    #enddo

    id vec(k1?,mu?)*vec(k2?,mu?) = k1.k2;

* divide by the normalizing factor of the denominator that is added to the topology
* this is always 1/(k^2 - m_UV^2)^3 = -i / (4 pi)^2 * 1/2 * 1/mUV^2
    Multiply i_ * rat(4 * 30246273033735921/9627687726852338,1)^2 * 2 * mUV^2;

* reduce the numerator
    repeat id k1?.k1?*uvprop(k1?,n1?) = uvprop(k1, n1-1) + mUV^2 * uvprop(k1, n1);
    id uvprop(k1?,n?) = uvprop(n);

* 1-loop IBP
    id uvprop(n1?{<1}) = 0;
    repeat id uvprop(n1?{>1}) = uvprop(-1 + n1)*rat((2 + D - 2*n1), 2 * mUV^2 * (-1 + n1));
    id uvprop(1) = mi1L1 * rat(2*mUV^2, D - 2);

* TODO: take to the power of loops
* normalize with 1/(4 pi e^-gamma)
    Multiply 1 - rat(53646286447601093/27457288774331243 * ep  + 7812755848557151/4093268398007683 * ep^2 - 4523275530886483/3638800576560925 * ep^3, 1);

    Multiply replace_(D, 4 - 2 * ep);

	id mi1L1 = (
        cMi1L1EpsM1logmUV0*rat(logmUV^0,1)*rat(1,ep^1)
        +cMi1L1Eps0logmUV0*rat(logmUV^0,1)*rat(ep^0,1)
        +cMi1L1Eps0logmUV1*rat(logmUV^1,1)*rat(ep^0,1)
        +cMi1L1Eps1logmUV0*rat(logmUV^0,1)*rat(ep^1,1)
        +cMi1L1Eps1logmUV1*rat(logmUV^1,1)*rat(ep^1,1)
        +cMi1L1Eps1logmUV2*rat(logmUV^2,1)*rat(ep^1,1)
        +cMi1L1Eps2logmUV0*rat(logmUV^0,1)*rat(ep^2,1)
        +cMi1L1Eps2logmUV1*rat(logmUV^1,1)*rat(ep^2,1)
        +cMi1L1Eps2logmUV2*rat(logmUV^2,1)*rat(ep^2,1)
        +cMi1L1Eps2logmUV3*rat(logmUV^3,1)*rat(ep^2,1)
        +cMi1L1Eps3logmUV0*rat(logmUV^0,1)*rat(ep^3,1)
        +cMi1L1Eps3logmUV1*rat(logmUV^1,1)*rat(ep^3,1)
        +cMi1L1Eps3logmUV2*rat(logmUV^2,1)*rat(ep^3,1)
        +cMi1L1Eps3logmUV3*rat(logmUV^3,1)*rat(ep^3,1)
        +cMi1L1Eps3logmUV4*rat(logmUV^4,1)*rat(ep^3,1)
        +alarmMi1L1*rat(ep^4,1)
	);

    id cMi1L1EpsM1logmUV0 = rat(132049792606502*i_,20852467428353107);
    id cMi1L1Eps0logmUV0 = rat(930733834984243*i_,75225176822851668);
    id cMi1L1Eps0logmUV1 = rat(-132049792606502*i_,20852467428353107);
    id cMi1L1Eps1logmUV0 = rat(238365146153033*i_,13782143439795685);
    id cMi1L1Eps1logmUV1 = rat(-1057101110939383*i_,85438623805262522);
    id cMi1L1Eps1logmUV2 = rat(66024896303251*i_,20852467428353107);
    id cMi1L1Eps2logmUV0 = rat(267137567007244*i_,17222977639926303);
    id cMi1L1Eps2logmUV1 = rat(-238365146153033*i_,13782143439795685);
    id cMi1L1Eps2logmUV2 = rat(31591818988785*i_,5106723491205427);
    id cMi1L1Eps2logmUV3 = rat(-26474513963629*i_,25084126035084908);
    id cMi1L1Eps3logmUV0 = rat(99723921272160*i_,7862278376418437);
    id cMi1L1Eps3logmUV1 = rat(-267137567007244*i_,17222977639926303);
    id cMi1L1Eps3logmUV2 = rat(238365146153033*i_,27564286879591370);
    id cMi1L1Eps3logmUV3 = rat(-10530606329595*i_,5106723491205427);
    id cMi1L1Eps3logmUV4 = rat(3349661396909*i_,12694975820195403);
endif;
.sort:integrated-ct-1;

* Factor out the mass
Polyratfun;
id rat(x1?,x2?) = num(x1)*den(x2);
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

* for now, just set the log piece to 0
id logmUV = 0;
.sort:integrated-ct-2;

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
#$OFFSET = 0;
#do i=1,`$MAXP'
    id penergy(p`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXK'
    id penergy(c`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
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
    Format O4,stats=off,saIter=`OPTIMITERATIONS';
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
