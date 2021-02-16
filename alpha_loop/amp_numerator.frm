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

Auto S mass;
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

Fill gyq(1) = yukawa_d; * d
Fill gyq(2) = yukawa_u; * u
Fill gyq(3) = yukawa_s; * s
Fill gyq(4) = yukawa_c; * c
Fill gyq(5) = yukawa_b; * b
Fill gyq(6) = yukawa_t; * t
Fill gyq(11) = 0; * e-
Fill gyq(12) = 0; * mu-
Fill gyq(13) = 0; * ta-
Fill gyq(-1) = yukawa_d; * d
Fill gyq(-2) = yukawa_u; * u
Fill gyq(-3) = yukawa_s; * s
Fill gyq(-4) = yukawa_c; * c
Fill gyq(-5) = yukawa_b; * b
Fill gyq(-6) = yukawa_t; * t
Fill gyq(-11) = 0; * e+
Fill gyq(-12) = 0; * mu+
Fill gyq(-13) = 0; * ta+

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = mass_b;
Fill masses(6) = mass_t;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = mass_b;
Fill masses(-6) = mass_t;
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
Fill masses(6) = mass_t;
Fill masses(-1) = mass_d;
Fill masses(-2) = mass_u;
Fill masses(-3) = mass_c;
Fill masses(-4) = mass_s;
Fill masses(-5) = mass_b;
Fill masses(-6) = mass_t;
Fill masses(11) = mass_e;
Fill masses(12) = mass_mu;
Fill masses(13) = mass_tau;
Fill masses(-11) = mass_e;
Fill masses(-12) = mass_mu;
Fill masses(-13) = mass_tau;
#endif

Fill masses(-82) = 0;
Fill masses(82) = 0;
Fill masses(21) = 0;
Fill masses(22) = 0;
Fill masses(25) = mass_h;
Fill masses(1337) = 0;

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
V p1,...,p40,k1,...,k40,c1,...,c40; * force this internal ordering in FORM
Auto V p,k,c, eps,ceps, sV,sVbar,sU,sUbar;
Auto S lm,ext;
#ifdef `FOURDIM'
    Auto I mu=4,s=4;
#else
    Auto I mu=D,s=D;
#endif
Symbol ge, gs, gy, ghhh, type, in, out, virtual;
Auto S y,x, idx, t, n;

Set dirac: s1,...,s80;
Set diracdummy: sd1,...,sd80;
Set lorentz: mu1,...,mu80;
Set lorentzdummy: mud1,...,mud80;

CF gamma, spinor ,NN, vector,g(s),delta(s),T, counter,color, prop, replace;
CF subs, configurations, conf, cmb, diag, der, energy, spatial(s);


CF f, vx, vec, vec1;
CF f, vx, vxs(s), vec, vec1;
CF subs, configurations, conf, cmb, diag, der, energy, spatial(s);
CF subgraph, uvconf, uvconf1, uvconf2, uvprop, uv, integrateduv;
CT gammatracetensor(c);


S UVRenormFINITE;
S ICT, mUV, logmu, logmUV, logmt, mi1L1, alarmMi1L1;
Fill logmasses(6) = logmt;
Fill logmasses(-6) = logmt;

CF integratedct, rat, num, den;
CF gammaCONF;

Set ts: t0,...,t20;
CT penergy;
NF energync;

S  i, m, n, ALARM;

CF hermconjugate, pol, cpol, uSpinor, ubarSpinor, vSpinor, vbarSpinor, sp(s), d, diracInd, lorentzInd,sunTF,sunTA,vec;
S ii,m,x,y,yy,xx;
CF spatialComp, spatialCompTmp;
Auto I i1, i2, j1,j2;


Polyratfun rat;

*--#] setup :

#include tensorreduce.frm
#include integrateduv.frm


* Load the diagrams
#include- input_`SGID'.h
.sort:load-input;
Hide CONF;




if (count(pol,1,cpol,1,uSpinor,1,vSpinor,1,ubarSpinor,1,vbarSpinor,1));
    Print "Unsubstituted polarization: %t";
    exit "Critical error";
endif;

* construct gamma string
id vec(p?,mu?) = p(mu);
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);
* treat gamma(p1+p2)
repeat id p(mu1?)*gamma(s1?,?a,mu1?,?b,s2?)= gamma(s1,?a,p1,?b,s2);
repeat; 
    id once sp(p?!vector_,x?) = p(mu)*sp(mu,x);
    id p?(mu?)*sp(mu?,x?) = sp(p,x);
endrepeat;
id sp(p1?,p2?) = p1.p2;
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

*************************************************
* Process different configurations (bubbles, etc)
*************************************************

* process all the configurations
L F = F * CONF;
.sort
Drop CONF;

* multiply the numerator contribution of derivatives
id conf(?a,p?) = conf(?a) * penergy(p);
id conf(?a,x?) = conf(?a) * x; * note the type difference

* Take a single derivative in pbubble_i^0
* FIXME: if the number of terms is 0 after differentiating, this configuration will not be created
repeat;
    id once der(p?) = der*replace_(p, p + pzero * x);
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

* TODO: do the local UV construction in the cmb
argument uv;
    id subgraph(?a,uvconf(?b,x?)) = subgraph(?a,uvconf(?b,uvdiag(x)));
* uvdiag is filled in from the table and yields a uvconf
    id subgraph(?a,uvconf(?b,uvconf(?c,x?))) = subgraph(?a,uvconf(?b,x));

   repeat;
* all subgraphs without dependencies can be treated at the same time
* multiply each graph with -1 to correctly subtract it
        id subgraph(x1?, x2?) = -uvconf1(x1, x2);
        argument uvconf1,2;
            id uvconf(x1?,x2?,?a,x3?) = uvconf2(x1) * tmax^x2 * uvconf1(?a) * x3;
            chainout uvconf1;
            repeat id uvconf1(p?)*uvconf1(p?) = uvconf1(p);
            id uvconf1(p?) = replace_(p, t * p);

            argument uvprop,1,vxs;
                id t = 1; * it could be that the LTD momentum also makes an appearance as an external momentum
            endargument;

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

            if (count(ALARM, 1));
                Print "UV Taylor expansion depth exceeded.";
                exit "";
            endif;

* select the right denominator structure
            id uvprop(?a,m?) = uvprop(?a);
            repeat id uvprop(k?,t1?,n1?)*uvprop(k?,t1?,n2?) = uvprop(k,t1,n1+n2);
            id uvprop(k?,t1?,n1?)*t1?^n2? = uvprop(k,n1 + n2)*t1^n2;
            id uvprop(k?,t1?ts,n?) = uvprop(k, n);
            id uvconf2(x?) = 1/x;
            id t?ts^n? = 0;
            id t = 1;
            id tmax = 1;

* collect all uv propagators of the subgraph
            Multiply replace_(vxs, vx);
            if (count(ICT, 1) == 0);
                id uvprop(?a) = 1;
                id vx(?a) = 1;
            else;
                id 1/ICT = ICT;
            endif;

            chainin uvprop;
            id uvprop(?a) = uvprop(?a, 1);
            repeat id vx(?a)*uvprop(?b,x?) = uvprop(?b, x*vx(?a));
            Multiply replace_(vx, vxs);
            id uvprop(?a) = integrateduv(?a);
        endargument;
        id uvconf1(?a) = uvconf(?a);

* now fill in the subgraph evaluation into the supergraph
        repeat id subgraph(x1?,?a,n?,?b,uvconf(?c,x2?))*uvconf(n?,x3?) = subgraph(x1,?a,?b,uvconf(?c,x2*x3));
    endrepeat;
   
    id uvconf(x?,x1?) = x1;
    if (count(subgraph, 1));
        Print "Unsubstituted UV subgraph: %t";
        exit "Critical error";
    endif;
endargument;
id uv(x?) = x;
.sort:local-uv-treatment;

* convert the polynomial to the cut momentum basis
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

* compute the integrated UV counterterm
Multiply counter(1);
repeat id k1?.k2?*counter(n?) = vec(k1,n)*vec(k2,n)*counter(n + 1);
* convert every k^0 into k.p0select, where p0select is effectively (1,0,0,0)
* TODO: is this safe in D dimensions?  
repeat id penergy(k1?)*counter(n?) = vec(k1,n)*vec(p0select,n)*counter(n + 1);
id counter(x?) = 1;

#do i=1,1
    id ICT = -1; * we add back the counterterm

    id once integrateduv(?a,x?) = uvprop(?a)*x;
    if (count(uvprop,1)) redefine i "0";

    Multiply replace_(vec, vec1); * consider all vectors as external
    repeat id uvprop(k1?,n?)*vec1(k1?,n1?) = uvprop(k1,n)*vec(k1,n1);

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

* If the expression is empty (due to epsilon pole selection), we still write a file
#if ( termsin(F) == 0 )
    #write<out_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
    #write<out_integrand_PF_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
    #write<out_integrand_LTD_`SGID'.proto_c> "#0 due to epsilon pole selection\n"
#endif
*Print +ss;
#ifdef `INTEGRAND'
    #if (`INTEGRAND' == "LTD") || (`INTEGRAND' == "both")
        .sort
        L FINTEGRANDLTD = F;
        .sort
        #include- ltdtable_`SGID'.h
        .sort:load-ltd;
        Hide F;

        id conf(n?,n1?,?a) = conf(n,n1,?a)*ltdtopo(n);
        id ltdtopo(n?) = 0; * unrecognized topology

        id cmb(?a) = cmb(?a)*replace(?a);
* gradually transform the basis to prevent term blow-up
        #do i=1,1
            id replace = 1;
            if (count(replace,1)) redefine i "0";

            AB+ cmb;
            .sort:cmb-1;
            Keep brackets; * make sure cmb is not replaced

            id replace(p1?,p2?,?a) = replace_(p1,p2)*replace(?a);
        #enddo
        .sort
        UnHide F;
    #endif
#endif


* now extract the energy components of the LTD loop variables
id k1?.k2? = g(k1, k2);
repeat id conf(x?,?a,k1?,?b)*g(k1?,k2?) = conf(x,?a,k1,?b)*f(k1,k2);
id g(p1?,p2?) = p1.p2;
symmetrize f;
id f(?a) = f(?a,1);
id f(?a,n1?)*f(?a,n2?) = f(?a,n1+n2);

B+ f;
.sort:energy-splitoff-1;
Keep brackets;
id f(k1?,k2?,n?) = (penergy(k1)*penergy(k2)-spatial(k1,k2))^n;

B+ conf,penergy,energy;
.sort:energy-splitoff-2;
Keep brackets;

repeat id conf(x?,?a,k1?,?b)*penergy(k1?) = conf(x,?a,k1,?b)*energy(k1);

repeat id energy(?a)*energy(?b) = energy(?a,?b);
symmetrize energy;
id energy(?a) = energy(f(?a));
if (count(energy,1) == 0) Multiply energy(f(c0)); * signal with c0 that we are dealing with the constant term
;
.sort:energy-splitoff-3;


*********** EXPAND GAMMA CHAINS INTO EXPLICIT COMPONENTS **********************
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

.sort
repeat id once gamma(?a,mud1?,?b)*gamma(?c,mud1?,?d) = sum_(y,0,3,gamma(?a,y,?b)*gamma(?c,y,?d));
repeat id once spinor(p?,sd?diracdummy)*gamma(sd?diracdummy,?a) = sum_(y,1,4,spinor(p,y)*gamma(y,?a));
repeat id once spinor(p?,sd?diracdummy)*gamma(?a,sd?diracdummy) = sum_(y,1,4,spinor(p,y)*gamma(?a,y));
repeat id once spinor(p?,s1?dirac)*gamma(s1?dirac,?a) = sum_(y,1,4,spinor(p,y)*gamma(y,?a));
repeat id once spinor(p?,s1?dirac)*gamma(?a,s1?dirac) = sum_(y,1,4,spinor(p,y)*gamma(?a,y));

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
*       gamma(1,p1+q1,2) expansion              
        id p?(mu?)*gamma(xx?number_,mu?,yy?number_) = gamma(xx,p,yy);
        id once ifmatch->retry4 gamma(xx?number_,?a,p?!vector_,yy?number_) = p(mu)*gamma(xx,?a,mu,yy);
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

.sort:gamma-chain-explicit;



*********************************************
* Construction of the integrand
*********************************************
#ifdef `INTEGRAND'
    Hide F;
    .sort

#if (`INTEGRAND' == "LTD") || (`INTEGRAND' == "both")
    id diag(?a) = 1;
    id cmb(?a) = 1;

* transform the LTD energies back to normal energies
    id energy(f(?a)) = energy(?a);
    chainout energy;
    id energy(c0) = 1;

    repeat id ltdcbtolmb(?a,c1?,x?,?b)*energy(c1?) = x*ltdcbtolmb(?a,c1,x,?b);
    id ltdcbtolmb(?a) = 1;

    if (count(energy,1));
        Print "Unsubstituted energies: %t";
        exit "Critical error";
    endif;

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
id prop(n?,x?,?a) = invdset[x]^(-n);

* set the ltd energies (including cut sign)
id ltdenergy(?a) = replace_(?a);


argument ellipsoids;
    id energies(0) = 0;
endargument;

.sort:integrand-ltd;

#if (`INTEGRAND' == "both")
    Hide FINTEGRANDLTD;
#endif
#endif

#if (`INTEGRAND' == "PF") || (`INTEGRAND' == "both")
    .sort
    #include- pftable_`SGID'.h
    L FINTEGRANDPF = F;

    B energy,diag,cmb;
    .sort:load-pf;
    Keep brackets;

    id energy(f(?a)) = energy(?a);
    chainout energy;
    id energy(c0) = 1;

* map all diagrams to their unique representative
    id diag(x1?,x2?,?a) = diag(pfmap(x1,x2),?a);
    id diag(diag(?a),?b) = diag(?a,?b);

* collect all the energies in the diagram
    repeat id diag(?a,p1?,p2?,?b) = diag(?a,p1,0,p2,?b);
    id diag(?a,p1?) = diag(?a,p1,0);

    repeat id energy(k?)*diag(?a,k?,n?,?b) = diag(?a,k,n+1,?b);

    if (count(energy,1));
        Print "Energy left: %t";
    endif;

    if (count(cmb,1) == 0);
        Print "FAIL: %t";
    endif;

    repeat id cmb(?a)*diag(?b) = f(diag(?b,cmb(?a)))*cmb(?a);
    id f(?a) = diag(?a);
    argtoextrasymbol tonumber,diag,1;
    #redefine oldextrasymbols "`extrasymbols_'"

    .sort:diag-1;
    #redefine diagend "`extrasymbols_'"
    #do ext={`oldextrasymbols'+1},`diagend'
        L diag{`ext'-`oldextrasymbols'} = extrasymbol_(`ext');
    #enddo
    #define diagcount "{`diagend'-`oldextrasymbols'}"

    id diag(x1?,x2?,?a,x3?) = diag(?a)*pftopo(x1,x2)*x3*conf(-1,x1,x2);
    id cmb(?a) = replace_(?a);
    .sort:pf-splitoff;
    Hide FINTEGRANDPF;


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

    .sort

    UnHide FINTEGRANDPF;
    
    .sort
    Drop diag1,...,diag`diagcount';

* now add all PF structures as a special conf
    #if `diagcount' > 0
        L FINTEGRANDPF = FINTEGRANDPF + <diag1*conf(-1)>+...+<diag`diagcount'*conf(-`diagcount')>;
    #endif
    id conf(x?)*conf(-1,?a) = conf(x,?a);

    #if (`SUMDIAGRAMSETS' == "onlysum")
        id conf(x?{>=0},x1?,?a) = conf(1000 + x1);
    #elseif (`SUMDIAGRAMSETS' == "both")
        id conf(x?{>=0},x1?,?a) = conf(x,?a) + conf(1000 + x1);
    #endif

    .sort
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

    repeat id allenergies(p?,?a) = energync(p.p)*allenergies(?a);
    argument energync;
        id p1?.p2? = spatial(p1, p2);
    endargument;
    chainin energync;
    id energync(?a)*allenergies = allenergies(?a);
    id allenergies(?a) = energies(?a);

    repeat id constants(p?,?a) = energync(p.p)*constants(?a);
    chainin energync;
    id energync(?a)*constants = constants(?a);

    argument ellipsoids, constants;
        id energy(p?) = penergy(p);
    endargument;
    B+ penergy,spatial,energies,allenergies, ellipsoids, constants;
*    print +ss;
    .sort
    Keep brackets;
***************** TRANSLATION OF FUNCTIONS INVOLVING ONLY MOMENTA (as in standard LTD)
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
********************** TRANSLATION OF ALL ADDITONAL OBJECTS ONLY EXISTING IN AMPLITUDES (spinors, spatial components etc)
* for polarized cross-sections
    #$MAXEPS =  `NPOL';
    #$MAXCEPS = `NCPOL';
    #$MAXV =`NSPINV';
    #$MAXVBAR=  `NSPINVBAR';
    #$MAXU  = `NSPINU';
    #$MAXUBAR = `NSPINUBAR';
* spatial components of momenta    
    #do i=1,`$MAXP'
        #do j =1,3
            id spatialComp(p`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
    #enddo
    #do i=1,`$MAXK'
        #do j =1,3
            id spatialComp(c`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
* polarization vectors
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
            id ceps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo    
    #enddo
* spinors    
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
    .sort:conv-amp-specific;
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
        .sort
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

            #if (`INTEGRANDTYPE' == "PF")
                $isdenominator = 0;
            #else
                $isdenominator = 1;
            #endif
            inside $conf;
                if (match(conf(x?{<0},?a))) $isdenominator = 1;
            endinside;
            $ellipsoids = 0;
            $energies = 0;
            $constants = 0;
            id ellipsoids(?a$ellipsoids) = 1;
            id energies(?a$energies) = 1;
            id constants(?a$constants) = 1;
            .sort
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
        .sort
        Drop FINTEGRAND`INTEGRANDTYPE';
        delete extrasymbols>`oldextrasymbols';
    #endif
    #enddo
    Unhide F;
#endif

*********************************************
* Construction of optimized numerator C code
*********************************************
* Convert the dot products and energies to a symbol
  
.sort

*********************************************
* Construction of optimized numerator C code
*********************************************

*#if ((isdefined(NUMERATOR)) && (`NUMERATOR' == 1))

id diag(?a) = 1;
id cmb(?a) = 1;

***************** TRANSLATION OF FUNCTIONS INVOLVING ONLY MOMENTA (as in standard LTD)
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
        id spatial(p`i', p`j') = lm`$OFFSET';
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
    
********************** TRANSLATION OF ALL ADDITONAL OBJECTS ONLY EXISTING IN AMPLITUDES (spinors, spatial components etc)
* for polarized cross-sections
    #$MAXEPS =  `NPOL';
    #$MAXCEPS = `NCPOL';
    #$MAXV =`NSPINV';
    #$MAXVBAR=  `NSPINVBAR';
    #$MAXU  = `NSPINU';
    #$MAXUBAR = `NSPINUBAR';
* spatial components of momenta    
    #do i=1,`$MAXP'
        #do j =1,3
            id spatialComp(p`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo       
    #enddo
    #do i=1,`$MAXK'
        #do j =1,3
            id spatialComp(c`i',`j') =  lm`$OFFSET';
            #$OFFSET = $OFFSET + 1;
        #enddo
    #enddo
* polarization vectors
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
            id ceps`i'.p`j' = lm`$OFFSET';
            #$OFFSET = $OFFSET+1;
        #enddo    
    #enddo
* spinors    
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
B+ conf;
.sort:conf-collect;
Keep brackets;

id conf(?a) = conf(conf(?a));
argtoextrasymbol tonumber,conf,1;
#redefine oldextrasymbols "`extrasymbols_'"
B+ conf;
.sort:conf-1;
Hide F;

#redefine energysymbolstart "`extrasymbols_'"
#do ext={`oldextrasymbols'+1}, `energysymbolstart'
    #write "START numerator"
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
    Format O1;*`OPTIMLVL',stats=on,method=`OPTIMISATIONSTRATEGY',saIter=`OPTIMITERATIONS';
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
    #write "END numerator"
#enddo
*#else
*    #write<out_`SGID'.proto_c> "#not generated\n"
*#endif
.end
