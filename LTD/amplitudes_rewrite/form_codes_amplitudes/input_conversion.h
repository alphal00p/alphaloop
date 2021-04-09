#procedure scalar-prop-to-clTD-prop
*replace propagator
    id sprop(p?,y2?) = prop(penergy(p), sqrt(spatial(p,p) + y2^2) );
*chain out energies
    argument prop;
        id penergy(p?!vector_) = p(mu)*penergy(mu);
        id p?(mu) * penergy(mu) = penergy(p);
    endargument;
*chain out spatial dots
    argument prop;
        argument sqrt;
            id spatial(p?!vector_,p1?) = p(mu1)*spatial(mu1,p1);
            id spatial(p1?,p?!vector_) = p(mu2)*spatial(mu2,p1);
            id spatial(p?!vector_,k?!vector_) = p(mu3)*k(mu4)*spatial(mu3,mu4);        
            id p?(mu1?)*k?(mu2?)* spatial(mu1?,mu2?) = spatial(p,k);
            id p?(mu?) * spatial(mu?,p1?) = spatial(p,p1);
            symmetrize spatial;
        endargument;
    endargument;
    B+ sp;
    .sort:treat-spatial-dots;
*chain out scalar products
    Keep brackets;
    repeat;
        id, once, ifmatch -> simpl, sp(p?!vector_,p1?) = p(mud1)*sp(mud1,p1);
        id, once, ifmatch -> simpl, sp(p1?,p?!vector_) = p(mud2)*sp(mud2,p1);
        id, once, ifmatch -> simpl, sp(p?!vector_,k?!vector_) = p(mud3)*k(mud4)*sp(mud3,mud4);        
        label simpl;
        id p?(mud1?)*k?(mud2?)*sp(mud1?,mud2?) = sp(p,k);
        id p?(mud?)*sp(mud?,p1?) = sp(p,p1);
    endrepeat;
#endprocedure

#procedure translate-inds-and-ext-data
* transform indices
    Argument; 
        id indS(aa?) = indS(aa+`INDSHIFT');
        id muL(aa?)  = muL(aa+`INDSHIFT');
    EndArgument;
    Argument;
        id indS(aa?) = dirac[aa];
        id muL(aa?)  = lorentz[aa];
    EndArgument;

* transform external data
    #do i = 1,`NPOL'
        id pol(`i',?aa,ind1?) = lVec(eps`i',ind1);    
    #enddo
    #do i = 1,`NCPOL'
        id cpol(`i',?aa,ind1?) = lVec(ceps`i',ind1);
    #enddo
    #do i = 1,`NSPINU'
        id uSpinor(`i',?aa,ind1?) = spinor(sU`i',ind1);
    #enddo
    #do i = 1,`NSPINV'
        id vSpinor(`i',?aa,ind1?) = spinor(sV`i',ind1);
    #enddo
    #do i = 1,`NSPINUBAR'
        id ubarSpinor(`i',?aa,ind1?) = spinor(sUbar`i',ind1);
    #enddo
    #do i = 1,`NSPINVBAR'
        id vbarSpinor(`i',?aa,ind1?) = spinor(sVbar`i',ind1);
    #enddo
    .sort:replace-ext-data;
    if (count(pol,1,cpol,1,uSpinor,1,vSpinor,1,ubarSpinor,1,vbarSpinor,1));
        Print "Unsubstituted polarization: %t";
        exit "Critical error";
    endif;

*transform functions:
    id lVec(p?,mu?) = p(mu);
    id gMetric(mu1?,mu2?) = d_(mu1,mu2);
    id deltaS(s1?,s2?) = d_(s1,s2);
    argument;
        id sp(p1?,p2?) = p1.p2;
    endargument;
#endprocedure