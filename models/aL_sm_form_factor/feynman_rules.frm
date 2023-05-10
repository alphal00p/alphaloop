******************************************
* Subprocedures for FeynmanRulesGlobal() *
******************************************

* split quartic vertex
#procedure SplitQuarticVertex()

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

#endprocedure

* do the spin sum external particles
#procedure SpinSum()
    
    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = 1;
    repeat id prop(`GLU', in, p?, idx1?)*prop(`GLU', out, p?, idx2?) = d_(colA[idx1], colA[idx2]);
    repeat id prop(`Z', in, p?, idx1?)*prop(`Z', out, p?, idx2?) = 1;
    repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = 1;
    repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
    repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = 1;
    repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = d_(colF[idx1], colF[idx2]);

#endprocedure

* virtual edges
#procedure VirtualEdges()

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

#endprocedure

**************************************************
* START SE prop couplings Feynman rules
**************************************************
#procedure SEPropCouplings()

    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = 1;
    repeat id prop(x1?{`QBARMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVE'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2]);
    repeat id prop(x1?{`QBARMASSIVE'}, in, p?, idx1?)*prop(x2?{`QBARMASSIVEPRIME'}, out, p?, idx2?) = d_(colF[idx1], colF[idx2]);
    repeat id prop(x1?{`QMASSIVEPRIME'}, in, p?, idx1?)*prop(x2?{`QMASSIVE'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
    repeat id prop(x1?{`QMASSIVE'}, in, p?, idx1?)*prop(x2?{`QMASSIVEPRIME'}, out, p?, idx2?) = d_(colF[idx2], colF[idx1]);
    repeat id prop(`GLUPRIME', in, p?, idx1?)*prop(`GLUPRIME', out, p?, idx2?) = d_(colA[idx1], colA[idx2]);
    repeat id prop(`GHOPRIME', in, p?, idx1?)*prop(`GHOPRIME', out, p?, idx2?) = d_(colA[idx1], colA[idx2]);
    repeat id prop(`GHOPRIMEBAR', in, p?, idx1?)*prop(`GHOPRIMEBAR', out, p?, idx2?) = d_(colA[idx1], colA[idx2]);
    id prop(`SDUMMY', virtual, p?, idx1?, idx2?) = 1;
    id prop(`SDUMMY', in, p?, idx1?) = 1;
    id prop(`SDUMMY', out, p?, idx1?) = 1;
    id prop(`PHOPRIME', virtual, p?, idx1?, idx2?) = 1;
    id prop(x?{`QMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = i_ * d_(colF[idx2], colF[idx1]);
    id prop(x?{`QBARMASSIVEPRIME'}, virtual, p?, idx1?, idx2?) = - i_ * d_(colF[idx1], colF[idx2]);
    id prop(x?{`GHOPRIME',`GHOPRIMEBAR'}, virtual, p?, idx1?, idx2?) = - i_ *d_(colA[idx1], colA[idx2]);
    id prop(`GLUPRIME', virtual, p?, idx1?, idx2?) = - i_ * d_(colA[idx1], colA[idx2]);

#endprocedure
**************************************************
* END SE prop couplings Feynman rules
**************************************************

**************************************************
* START amp prop couplings Feynman rules
**************************************************
#procedure AmpPropCouplings()

    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOAMPPRIME', out, p?, idx2?) = 1;
    id prop(`PHOAMPPRIME', virtual, p?, idx1?, idx2?) = - i_;

#endprocedure
**************************************************
* START amp prop couplings Feynman rules
*************************************************


**************************************************
* START SE vx couplings Feynman rules
**************************************************
#procedure SEVxCouplings()

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
    id vx(`GLU', `GLUPRIME', `SDUMMY', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * d_(colA[idx1], colA[idx2]);
    id vx(`GHOPRIMEBAR', `SDUMMY', `GHOPRIMEBAR', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * d_(colA[idx1], colA[idx3]);
    id vx(`GHOPRIME', `SDUMMY', `GHOPRIME', p1?, p2?, p3?, idx1?, idx2?, idx3?) = i_ * d_(colA[idx1], colA[idx3]);

#endprocedure
**************************************************
* END SE vx couplings Feynman rules
**************************************************

**************************************************
* START amp vx couplings Feynman rules
**************************************************
#procedure AmpVxCouplings()

    id vx(x1?{`QBAR'}, `PHOAMPPRIME', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_ * d_(colF[idx1], colF[idx3]);
    id vx(x1?{`LBAR'}, `PHOAMPPRIME', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = charges(x2) * ge * i_;
    id vx(`PHOAMPPRIME', `PHOAMPPRIME', `PHO', `PHO', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?) = ge^4 * i_;

#endprocedure
**************************************************
* END amp vx couplings Feynman rules
**************************************************


* vertices
#procedure Vertices()

    id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * T(colF[idx1], colA[idx2], colF[idx3]);
*   * The version below is the text-book ghost Feynman rules, but we prefer the symmetrized version
*   * id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * i_ * colf(colA[idx3], colA[idx2], colA[idx1]);
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

*  * delta_Z vertex

*  * The first multiplicity factor is always the loop multiplicity factor! It must be adjusted w.r.t to n_f!

*  * dZ massless quark
    id vx(x1?{`QBARMASSLESS'}, x2?{`QMASSLESS'}, p1?, p2?, idx1?, idx2?) = (1/1) * (-1) * i_ * ((4/3)*gs^2/16/pi^2) * (1/ep) * d_(colF[idx1], colF[idx2]);

*   * the finite part needs to be checked, also because the factor 4/3 on the pole of the mass correction is pure fudge for now.
*   * dZ massive quark
    id vx(x1?{`QBARMASSIVE'}, x2?{`QMASSIVE'}, p1?, p2?, idx1?, idx2?) = (1/1) * (-1) * i_ * ((4/3)*gs^2/16/pi^2) * d_(colF[idx1], colF[idx2]);

*   * dZ gluon

*   * The version below is for contributions to the gluon wavefunction from g, gh and down quark only, so it is good for e+ e- > j j j / u c s b t
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

#endprocedure

********************************************
* subprocedures for FeynmanRulesMomentum() *
********************************************

* do the spin sum external particles
#procedure SpinSumMomentum()

    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHO', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
    repeat id prop(`GLU', in, p?, idx1?)*prop(`GLU', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
    repeat id prop(`Z', in, p?, idx1?)*prop(`Z', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
    repeat id prop(x?{`L'}, in, p?, idx1?)*prop(x?{`L',}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
    repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
    repeat id prop(x?{`Q'}, in, p?, idx1?)*prop(x?{`Q'}, out, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) + masses(x)*gamma(dirac[idx1], dirac[idx2]);
    repeat id prop(x?{`LBAR'}, out, p?, idx1?)*prop(x?{`LBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);
    repeat id prop(x?{`QBAR'}, out, p?, idx1?)*prop(x?{`QBAR'}, in, p?, idx2?) = gamma(dirac[idx1], p, dirac[idx2]) - masses(x)*gamma(dirac[idx1], dirac[idx2]);

#endprocedure

* virtual edges 
#procedure VirtualEdgesMomentum()

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

#endprocedure

**************************************************
* START SE prop Lorentz Feynman rules
**************************************************
#procedure SEPropLorentzFeynmanRules()

*   * The original spin-sum can be used:
*   *repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
*   * or one can include a projector, like below
    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOPRIME', out, p?, idx2?) = -d_(lorentz[idx1], lorentz[idx2]) + energyselector(lorentz[idx1]) * energyselector(lorentz[idx2]);
*   * same for all repeat ID below:
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
    id prop(`GLUPRIME', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
    id prop(x?{`GHOPRIME',`GHOPRIMEBAR'}, virtual, p?, idx1?, idx2?) = 1;

#endprocedure
**************************************************
* END SE prop Lorentz Feynman rules
**************************************************

**************************************************
* START amp prop Lorentz Feynman rules
**************************************************
#procedure AmpPropLorentzFeynmanRules()

    id prop(`PHOAMPPRIME', virtual, p?, idx1?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);
    repeat id prop(`PHO', in, p?, idx1?)*prop(`PHOAMPPRIME', out, p?, idx2?) = d_(lorentz[idx1], lorentz[idx2]);

#endprocedure
**************************************************
* END amp prop Lorentz Feynman rules
**************************************************

**************************************************
* START SE vx Lorentz Feynman rules
**************************************************
#procedure SEVxLorentzFeynmanRules()

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
    id vx(`GLU', `GLUPRIME', `SDUMMY', p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(lorentz[idx1], lorentz[idx2]);
    id vx(`GHO', `SDUMMY', `GHOPRIME', p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(lorentz[idx1], lorentz[idx3]);
    id vx(`GHOBAR', `SDUMMY', `GHOPRIMEBAR', p1?, p2?, p3?, idx1?, idx2?, idx3?) = d_(lorentz[idx1], lorentz[idx3]);

#endprocedure 
**************************************************
* END SE vx Lorentz Feynman rules
**************************************************

**************************************************
* START amp vx Lorentz Feynman rules
**************************************************


#procedure AmpVxLorentzFeynmanRules()
    id vx(x1?{`QBAR'}, `PHOAMPPRIME', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
    id vx(x1?{`LBAR'}, `PHOAMPPRIME', x2?{`L'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
    id vx(`PHOAMPPRIME', `PHOAMPPRIME', `PHO', `PHO', p4?, p3?, p2?, p1?, idx4?, idx3?, idx2?, idx1?) = 
    ( APHOAMPFFSTU(p1,p2,p3) ) * (d_(lorentz[idx1], lorentz[idx2]) - (2*FFS(p1, p2, p3)^-1)*p1(lorentz[idx2])*p2(lorentz[idx1])  ) * ( d_(lorentz[idx3], lorentz[idx4]) + (2*FFS(p1, p2, p3)^-1)*p1(lorentz[idx3])*p3(lorentz[idx4]) + (2*FFS(p1, p2, p3)^-1)*p2(lorentz[idx3])*p3(lorentz[idx4]) + (2*FFS(p1, p2, p3)^-1)*p3(lorentz[idx3])*p3(lorentz[idx4]))
+ ( APHOAMPFFTSU(p1,p2,p3) ) * (d_(lorentz[idx3], lorentz[idx2]) - (2*FFT(p1, p2, p3)^-1)*p3(lorentz[idx2])*p2(lorentz[idx3])  ) * ( d_(lorentz[idx1], lorentz[idx4]) + (2*FFT(p1, p2, p3)^-1)*p3(lorentz[idx1])*p1(lorentz[idx4]) + (2*FFT(p1, p2, p3)^-1)*p2(lorentz[idx1])*p1(lorentz[idx4]) + (2*FFT(p1, p2, p3)^-1)*p1(lorentz[idx1])*p1(lorentz[idx4]))
+ ( APHOAMPFFUST(p1,p2,p3) ) * (d_(lorentz[idx3], lorentz[idx1]) - (2*FFU(p1, p2, p3)^-1)*p3(lorentz[idx1])*p1(lorentz[idx3])  ) * ( d_(lorentz[idx2], lorentz[idx4]) + (2*FFU(p1, p2, p3)^-1)*p3(lorentz[idx2])*p2(lorentz[idx4]) + (2*FFU(p1, p2, p3)^-1)*p1(lorentz[idx2])*p2(lorentz[idx4]) + (2*FFU(p1, p2, p3)^-1)*p2(lorentz[idx2])*p2(lorentz[idx4]))
+ ( BPHOAMPFFSTU(p1,p2,p3) ) * (( - ( d_(lorentz[idx1], lorentz[idx2]) - (2*FFS(p1, p2, p3)^-1)*p1(lorentz[idx2])*p2(lorentz[idx1]) ) * ( p1(lorentz[idx3]) - (FFU(p1, p2, p3)*FFT(p1, p2, p3)^-1) * p2(lorentz[idx3]) ) * ( p2(lorentz[idx4]) - (FFU(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p3(lorentz[idx4]) ) ) 
+ (p3(lorentz[idx1]) - (FFU(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p2(lorentz[idx1]) )* (p1(lorentz[idx2]) - (FFS(p1, p2, p3)*FFT(p1, p2, p3)^-1) * p3(lorentz[idx2])) * (d_(lorentz[idx3], lorentz[idx4]) + (2*FFT(p1, p2, p3)^-1) * p2(lorentz[idx3]) * p1(lorentz[idx4]) + (2*FFU(p1, p2, p3)^-1)*p1(lorentz[idx3]) * p2(lorentz[idx4]) 
+ (2*FFS(p1, p2, p3)^-1) * (p1(lorentz[idx3]) + p2(lorentz[idx3]) + p3(lorentz[idx3])) * (p1(lorentz[idx4]) + p2(lorentz[idx4]) + p3(lorentz[idx4])) + p3(lorentz[idx3])*p3(lorentz[idx4])) ) 
+ ( BPHOAMPFFTSU(p1,p2,p3) ) * (( - ( d_(lorentz[idx3], lorentz[idx2]) - (2*FFT(p1, p2, p3)^-1)*p3(lorentz[idx2])*p2(lorentz[idx3]) ) * ( p3(lorentz[idx1]) - (FFU(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p2(lorentz[idx1]) ) * ( p2(lorentz[idx4]) - (FFU(p1, p2, p3)*FFT(p1, p2, p3)^-1) * p1(lorentz[idx4]) ) ) 
+ (p1(lorentz[idx3]) - (FFU(p1, p2, p3)*FFT(p1, p2, p3)^-1) * p2(lorentz[idx3]) )* (p3(lorentz[idx2]) - (FFT(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p1(lorentz[idx2])) * (d_(lorentz[idx1], lorentz[idx4]) + (2*FFS(p1, p2, p3)^-1) * p2(lorentz[idx1]) * p3(lorentz[idx4]) + (2*FFU(p1, p2, p3)^-1)*p3(lorentz[idx1]) * p2(lorentz[idx4]) 
+ (2*FFT(p1, p2, p3)^-1) * (p3(lorentz[idx1]) + p2(lorentz[idx1]) + p1(lorentz[idx1])) * (p3(lorentz[idx4]) + p2(lorentz[idx4]) + p1(lorentz[idx4])) + p1(lorentz[idx1])*p1(lorentz[idx4])) ) 
+ ( BPHOAMPFFUST(p1,p2,p3) ) * (( - ( d_(lorentz[idx3], lorentz[idx1]) - (2*FFU(p1, p2, p3)^-1)*p3(lorentz[idx1])*p1(lorentz[idx3]) ) * ( p3(lorentz[idx2]) - (FFT(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p1(lorentz[idx2]) ) * ( p1(lorentz[idx4]) - (FFT(p1, p2, p3)*FFU(p1, p2, p3)^-1) * p2(lorentz[idx4]) ) ) 
+ (p2(lorentz[idx3]) - (FFT(p1, p2, p3)*FFU(p1, p2, p3)^-1) * p1(lorentz[idx3]) )* (p3(lorentz[idx1]) - (FFU(p1, p2, p3)*FFS(p1, p2, p3)^-1) * p2(lorentz[idx1])) * (d_(lorentz[idx2], lorentz[idx4]) + (2*FFS(p1, p2, p3)^-1) * p1(lorentz[idx2]) * p3(lorentz[idx4]) + (2*FFT(p1, p2, p3)^-1)*p3(lorentz[idx2]) * p1(lorentz[idx4]) 
+ (2*FFU(p1, p2, p3)^-1) * (p3(lorentz[idx2]) + p1(lorentz[idx2]) + p2(lorentz[idx2])) * (p3(lorentz[idx4]) + p1(lorentz[idx4]) + p2(lorentz[idx4])) + p2(lorentz[idx2])*p2(lorentz[idx4])) ) 
+ ( CPHOAMPFFSTU(p1,p2,p3)*FFS(p1, p2, p3)^-1*FFT(p1, p2, p3)^-1*FFT(p1, p2, p3)^-1*FFU(p1, p2, p3)^-1) * ( FFS(p1, p2, p3) * p3(lorentz[idx1]) - FFU(p1, p2, p3) * p2(lorentz[idx1]) ) * ( FFS(p1, p2, p3) * p3(lorentz[idx2]) - FFT(p1, p2, p3) * p1(lorentz[idx2]) ) * ( FFU(p1, p2, p3) * p2(lorentz[idx3]) - FFT(p1, p2, p3) * p1(lorentz[idx3]) ) * ( FFS(p1, p2, p3) * p2(lorentz[idx4]) - FFU(p1, p2, p3) * p3(lorentz[idx4]) ) ;
#endprocedure
**************************************************
* END amp vx Lorentz Feynman rules
**************************************************

**************************************************
* START amp FF substitution
**************************************************
#procedure AmpFFSubstitution()
    CF FFSINV, FFTINV, FFUINV;
    id FFS(p1?, p2?, p3?) = 2*p1.p2;
    id FFT(p1?, p2?, p3?) = 2*p2.p3;
    id FFU(p1?, p2?, p3?) = 2*p1.p3;
    id FFS(p1?, p2?, p3?)^-1 = FFSINV(0,0,p1.p2);
    id FFT(p1?, p2?, p3?)^-1 = FFTINV(0,0,p2.p3);
    id FFU(p1?, p2?, p3?)^-1 = FFUINV(0,0,p1.p3);
    id APHOAMPFFSTU(p1?, p2?, p3?) = APHOAMPFFSTU(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id APHOAMPFFTSU(p1?, p2?, p3?) = APHOAMPFFTSU(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id APHOAMPFFUST(p1?, p2?, p3?) = APHOAMPFFUST(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id BPHOAMPFFSTU(p1?, p2?, p3?) = BPHOAMPFFSTU(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id BPHOAMPFFTSU(p1?, p2?, p3?) = BPHOAMPFFTSU(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id BPHOAMPFFUST(p1?, p2?, p3?) = BPHOAMPFFUST(0, 0, 0, p1.p2, p1.p3, p2.p3);
    id CPHOAMPFFSTU(p1?, p2?, p3?) = CPHOAMPFFSTU(0, 0, 0, p1.p2, p1.p3, p2.p3);
#endprocedure
**************************************************
* END amp FF substitution
**************************************************

* Replace these functions by 1, so the multiplicy check does not crash
#procedure MultiplictyPatch()
id APHOAMPFFSTU(?p) = 1;
id APHOAMPFFTSU(?p) = 1;
id APHOAMPFFUST(?p) = 1;
id BPHOAMPFFSTU(?p) = 1;
id BPHOAMPFFTSU(?p) = 1;
id BPHOAMPFFUST(?p) = 1;
id CPHOAMPFFSTU(?p) = 1;

id FFS(?p) = 1;
id FFT(?p)= 1;
id FFU(?p) = 1;

id FFSINV(?p) = 1;
id FFTINV(?p) = 1;
id FFUINV(?p) = 1; 
#endprocedure

* vertices

#procedure VerticesMomentum()

id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) ;
* The version below is the text-book ghost Feynman rules, but we prefer the symmetrized version
*id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = p3(lorentz[idx2]);
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

#endprocedure

*remaining gluon vertices
#procedure GluonVerticesMomentum()
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
#endprocedure