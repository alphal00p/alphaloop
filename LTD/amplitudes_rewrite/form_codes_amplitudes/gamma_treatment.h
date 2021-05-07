#procedure gamma - traces - chains
* chain out momenta gam(...,lVec(p1+p2),...) = gam(...,p1,...)+gam(...,p2,...)
    repeat id gam(?aa,lVec(p?vector_),?bb) = gam(?aa,p,?bb);
    repeat;
        id once gam(?aa,lVec(p?!vector_),?bb) = p(mu)*gam(?aa,mu,?bb);
        id p?(mu?) * gam(?aa, mu?,?bb) = gam(?aa, p,?bb);
    endrepeat;
    repeat id gam(?aa,s1?)*gam(s1?,?bb) = gam(?aa,?bb);
    id gam(s1?,?aa,s1?) = gammatrace(?aa);

    B + gammatrace;
    .sort : chain - out - momenta;

*perform traces Keep brackets;
    repeat;
        id once gammatrace(?aa) = g_(1,?aa);
        trace4, 1;
        endrepeat;
        B + gam;
        .sort:gamma-traces;
*simplify gamma - strings Keep brackets;
* these are all the chains which are not simplified
    id gam(?aa) = gammatensor(?aa);
#call Gstring4D(gammatensor, 0)
    id gammatensor(?aa) = gamma(?aa);
#endprocedure

#procedure expand-gamma-chains
***********expand gamma chains **********************CF gammaAll;
    S iter, intSym;

    Table explSpinor(1:4,p?);
    Fill explSpinor(1) = penergy(p);
    Fill explSpinor(2) = spatialComp(p, 1);
    Fill explSpinor(3) = spatialComp(p, 2);
    Fill explSpinor(4) = spatialComp(p, 3);

    .sort
* spinor expansion
    repeat id once gamma(?aa,mud1?,?bb)*gamma(?cc,mud1?,?dd) = sum_(y,0,3,gamma(?aa,y,?bb)*gamma(?cc,y,?dd));
    repeat id once spinor(p?,sd?diracdummy)*gamma(sd?diracdummy,?aa) = sum_(y,1,4,spinor(p,y)*gamma(y,?aa));
    repeat id once spinor(p?,sd?diracdummy)*gamma(?aa,sd?diracdummy) = sum_(y,1,4,spinor(p,y)*gamma(?aa,y));
    repeat id once spinor(p?,s1?dirac)*gamma(s1?dirac,?aa) = sum_(y,1,4,spinor(p,y)*gamma(y,?aa));
    repeat id once spinor(p?,s1?dirac)*gamma(?aa,s1?dirac) = sum_(y,1,4,spinor(p,y)*gamma(?aa,y));

    id spinor(p?,intSym?int_) = explSpinor(intSym,p);
    .sort

*    id gamma(?aa) = gammaAll(gamma(?aa));
    B + gamma;
    .sort

    #include- definition_gamma_explicit.h
***********************gamma expansion *********************************
***************************keep brackets;
*gamma(1, p1 + q1, 2) expansion 
    repeat;
        id once gamma(xx?number_,?aa,p?!vector_,aa?number_) = p(mu)*gamma(xx,?aa,mu,aa);
        id p?(mu?)*gamma(xx?number_, mu?, aa?number_) = gamma(xx, p, aa);
    endrepeat;
        B + gamma;
    .sort:p-expand;
        keep brackets;
    repeat;
        id gamma(x?int_, p?vector_,  y?int_) = slash(x,y,p);
        id gamma(x?int_, intSym?int_, y?int_) = gamtab(x,y,intSym);
        id gamma(x?int_,?aa, p?vector_,  y?int_) = sum_(iter,1,4,gamma(x,?aa,iter)*slash(iter,y,p));
        id gamma(x?int_,?aa, intSym?int_, y?int_) = sum_(iter,1,4,gamma(x,?aa,iter)*gamtab(iter,y,intSym));
    endrepeat;
    .sort:gamma-chain-explicit;

*sanity check 
    if (count(gamma, 1) || count(gam, 1));
        print "Some gammas are not replaced: %t";
        exit "Critical ERROR";
    endif;
#endprocedure
