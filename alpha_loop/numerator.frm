#-
#:ContinuationLines 0
On statistics;
On nospacesinnumbers;

*--#[ setup :

*import the model parameters, model_parameters.frm should be located in the model directory 
#include model_parameters.frm

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
CF gaugevector, invdot, invdotgauge, ffguard;

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

#include feynman_rules.frm

#procedure FeynmanRulesGlobal()
* extract the global factors from the Feynman rules, including colour

* construct a counter that yields new indices
id vx(?a,p?,idx1?,?b) = tmp(idx1,?b)*vx(?a,p,idx1,?b);
chainin tmp;
repeat id tmp(idx1?,idx2?,?a) = tmp(max_(idx1,idx2),?a);
id tmp(x?) = counter(x + 1);

* split up the quartic gluon vertex into distinct colour factors
#call SplitQuarticVertex()

id counter(x?) = 1;

* make a copy of the Feynman rules
id prop(?a) = prop(?a)*tmp(prop(?a));
id vx(?a) = vx(?a)*tmp(vx(?a));
repeat id tmp(x1?)*tmp(x2?) = tmp(x1*x2);

* strip momentum tags
repeat id f?{vx,prop}(?a,p?) = f(?a);

* do the spin sum external particles
#call SpinSum()

* virtual edges
#call VirtualEdges()

**************************************************
* START SE prop couplings Feynman rules
**************************************************
#call SEPropCouplings()
**************************************************
* END SE prop couplings Feynman rules
**************************************************

**************************************************
* START amp prop couplings Feynman rules
**************************************************
#call AmpPropCouplings()
**************************************************
* START amp prop couplings Feynman rules
**************************************************

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

**************************************************
* START SE vx couplings Feynman rules
**************************************************
#call SEVxCouplings()
**************************************************
* END SE vx couplings Feynman rules
**************************************************

**************************************************
* START amp vx couplings Feynman rules
**************************************************
#call AmpVxCouplings()
**************************************************
* END amp vx couplings Feynman rules
**************************************************

* vertices
#call Vertices()

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
#call SpinSumMomentum()

* virtual edges
#call VirtualEdgesMomentum()

**************************************************
* START SE prop Lorentz Feynman rules
**************************************************
#call SEPropLorentzFeynmanRules()
**************************************************
* END SE prop Lorentz Feynman rules
**************************************************

**************************************************
* START amp prop Lorentz Feynman rules
**************************************************
#call AmpPropLorentzFeynmanRules()

**************************************************
* END amp prop Lorentz Feynman rules
**************************************************

if (count(prop, 1));
    Print "Unsubstituted propagator: %t";
    exit "Critical error";
endif;

**************************************************
* START SE vx Lorentz Feynman rules
**************************************************
#call SEVxLorentzFeynmanRules()
**************************************************
* END SE vx Lorentz Feynman rules
**************************************************

**************************************************
* START amp vx Lorentz Feynman rules
**************************************************

#call AmpVxLorentzFeynmanRules()

**************************************************
* END amp vx Lorentz Feynman rules
**************************************************

* substitute the gauge vector
id gaugevector(p?, idx?) = energyselector(lorentz[idx]);
id invdotgauge(p?) = invdot(energyselector.p);

*this can be used to simplify the FF x FF
*id p1.p1 = 0;
*id p2.p2 = 0;
*id k1.k1 = 0;


.sort:feynman-rules-edges;

**************************************************
* START amp FF substitution
**************************************************
#call AmpFFSubstitution()
**************************************************
* END amp FF substitution
**************************************************

* vertices
#call VerticesMomentum()

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

* remaining gluon vertices
#call GluonVerticesMomentum()

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

*    id cten(mu1?, mu2?, mu3?, mu4?) = ctenfunc(mu1, mu2, mu3, mu4); * do I need to the other line as well?
    id e_(mu1?, mu2?, mu3?, mu4?) = levicivita(mu1, mu2, mu3, mu4); * the reverse substitution is delayed until after the call to IntegrateUV()
    repeat id levicivita(?a, k1?, ?b)*counter(n?) = levicivita(?a,n,?b)*vec(k1,n)*counter(n+1);
*    repeat id ctenfunc(?a, k1?, ?b)*counter(n?) = ctenfunc(?a,n,?b)*vec(k1,n)*counter(n+1);
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
    #if (`INTEGRAND' == "LTD") || (`INTEGRAND' == "both")
        Multiply replace_(fmb`i', fmb`i'*xorig + xsplit*energy(fmb`i')*energyselector - xsplit*fmbs`i');
    #else
        Multiply replace_(fmb`i', energy(fmb`i')*energyselector - fmbs`i');
    #endif
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

#ifdef `FORMFACTORS'
argument APHOAMPFFSTU, APHOAMPFFTSU, APHOAMPFFUST, BPHOAMPFFSTU, BPHOAMPFFTSU, BPHOAMPFFUST, CPHOAMPFFSTU;
    id energyselector.energyselector = 1;
    id energyselector.p?spatialparts = 0;
    id p?.energyselector = penergy(p);
    id p1?spatialparts.p2?spatialparts = -p1.p2; * add a -1 to fix the metric
    id p1?spatialparts.p? = spatial(p1,p);
endargument;
#endif

*id invdot(x?) * x? = 1;

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

* split compressed numerator per diag, as they may have different energy dependence
* create new expressions per bracket as the bracket content may be too large to fit in a term
    B conf,tder,cmb,diag;
    .sort:diagorig-splitoff-1;
    Hide rest;
    Keep brackets;

    putinside f;
    argtoextrasymbol tonumber f;
    #define bracketstart "`extrasymbols_'"
    B+ f;
    .sort:diagorig-splitoff-2;
    #define bracketend "`extrasymbols_'"
    #do i ={`bracketstart'+1},`bracketend'
        L diagorig`i' = diagorig[f(`i')]*extrasymbol_(`i');
    #enddo

    id f?{conf,cmb,tder}(?a) = 1;
    id diag(?a) = diag(?a,1);
    id p1?.p2? = dot(p1,p2,p1.p2);
    argument dot,3;
        Multiply replace_(<p1,ps1>,...,<p40,ps40>,<c1,cs1>,...,<c40,cs40>);
    endargument;
    .sort:diagorig-splitoff-3;
    Hide diagorig{`bracketstart'+1},...,diagorig`bracketend';

    #if {`bracketstart'+1} <= `bracketend'
        L diagorig = <extrasymbol_({`bracketstart'+1})*f(diag({`bracketstart'+1}),-1)>+...+<extrasymbol_(`bracketend')*f(diag(`bracketend'),-1)>;
    #else
        L diagorig =  0;
    #endif

    id diag(?a) = 1; * remove the diagram from the bracket
    id f(?a) = diag(?a);
    .sort:diagorig-splitoff-4;
    UnHide Fltd,forestltd1,...,forestltd`forestcount';
    .sort:diagorig-splitoff-5;
    Drop rest, diagorig;
    L Fltd = rest + diagorig;
    .sort:ltd-forest-copy;

    repeat id cmb(?a)*tder(?c)*diag(?b) = f(diag(?b,cmb(?a)*tder(?c)))*cmb(?a)*tder(?c);
    repeat id cmb(?a)*diag(?b) = f(diag(?b,cmb(?a)))*cmb(?a);
    id f(?a) = diag(?a);
    argtoextrasymbol tonumber,diag,1;
    #redefine diagstart "`extrasymbols_'"
    .sort:diag-1;
    Drop diagorig{`bracketstart'+1},...,diagorig`bracketend';
    #redefine diagend "`extrasymbols_'"
    #do ext={`diagstart'+1},`diagend'
        L diag{`ext'-`diagstart'} = extrasymbol_(`ext');
    #enddo
    #define diagcount "{`diagend'-`diagstart'}"

    #do i ={`bracketstart'+1},`bracketend'
        id diag(diag(`i'),-1,y?) = diagorig`i'*y*xorig;
    #enddo

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
* TODO: do per cut so that the invd array will not have zeroed entries
    argtoextrasymbol tonumber, den;
    #define denstart "`extrasymbols_'"
    .sort:invd-creation;
    #do i = {`denstart'+1},`extrasymbols_'
        id den(`i') = invd{`i'-`denstart'-1};

        #do j = 1,`diagcount'
            #define HASINVD`i'IN`j' "0"
            if (expression(diag`j') && count(invd{`i'-`denstart'-1}, 1)) redefine HASINVD`i'IN`j' "1";
        #enddo
    #enddo
    .sort:invd-collect;

    #do j = 1,`diagcount'
        if (expression(diag`j'));
            Multiply ellipsoids;
            #do i = {`denstart'+1},`extrasymbols_'
                #if `HASINVD`i'IN`j''
                    Multiply ellipsoids(invd{`i'-`denstart'-1}, extrasymbol_(`i'));
                #endif
            #enddo
        endif;
        chainin ellipsoids;
    #enddo

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

* Append the final momenta to the initial ones
#$MAXK = `NFINALMOMENTA';
#$MAXP = `NINITIALMOMENTA';
#$OFFSET = `$MAXP';
#do i= 1,`$MAXK'
    #$OFFSET = $OFFSET + 1;
    repeat id e_(p1?,p2?,p3?,cs`i') = penergy(p`$OFFSET')*e_(p1,p2,p3,energyselector) - e_(p1,p2,p3,p`$OFFSET');
****maybe I should do this:
*    repeat id cten(p1?, p2?, p3?, cs`i') = penergy(p`$OFFSET')*cten(p1,p2,p3,energyselector) - cten(p1,p2,p3,p`$OFFSET');
    multiply replace_(c`i', p`$OFFSET');
#enddo
#$MAXP = $MAXP + $MAXK;

.sort:func-prep;

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

* Convert the dot products, energies and levi civita tensors to a symbol, add your custom form factors to the list of argumetns

*this should eventually not have to be hardcoded anymore
#ifdef `FORMFACTORS'
    #$OFFSET = 0;
    #do i=1,`$MAXP'
        id penergy(p`i') = lm`$OFFSET';
        argument invdot, energies, ellipsoids, constants, dot, APHOAMPFFSTU, APHOAMPFFTSU, APHOAMPFFUST, BPHOAMPFFSTU, BPHOAMPFFTSU, BPHOAMPFFUST, CPHOAMPFFSTU;
            id penergy(p`i') = lm`$OFFSET';
        endargument;
        #$OFFSET = $OFFSET + 1;
        #do j=`i',`$MAXP'
            id p`i'.p`j'^-1 = lm`$OFFSET'^-1;
            argument invdot, energies, ellipsoids, constants, APHOAMPFFSTU, APHOAMPFFTSU, APHOAMPFFUST, BPHOAMPFFSTU, BPHOAMPFFTSU, BPHOAMPFFUST, CPHOAMPFFSTU;
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
#else
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
#endif
id invdot(x?) = x^-1;
.sort:conv-func;

#do i1=1,`$MAXP'
    #do i2={`i1'+1},`$MAXP'
        #do i3={`i2'+1},`$MAXP'
            id e_(p`i1',p`i2',p`i3',energyselector) = lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
            #do i4={`i3'+1},`$MAXP'
                id e_(p`i1',p`i2',p`i3',p`i4') = lm`$OFFSET';
                #$OFFSET = $OFFSET + 1;
            #enddo
        #enddo
    #enddo
#enddo

*#do i1=1,`$MAXP'
*    #do i2=1,`$MAXP'
*        #do i3=1,`$MAXP'
*            #do i4=1,`$MAXP'
*                id cten(p`i1',p`i2',p`i3',p`i4') = lm`$OFFSET';
*                #$OFFSET = $OFFSET + 1;
*            #enddo
*        #enddo
*    #enddo
*#enddo

#$OFFSET = 0;
#do i=1,`$MAXP'
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
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
            Format 255;
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
