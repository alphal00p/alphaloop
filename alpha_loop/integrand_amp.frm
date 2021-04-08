#-
* I expect the following input:
* OPTIMISATIONSTRATEGY=CSEgreedy 
* SGID=0 
* NINITIALMOMENTA=2 
* really means final independent final momenta!
* NFINALMOMENTA=3 -D 
* really means loop-momenta!
* NLOOPMOMENTA
* NPOL=2 
* NCPOL=0 
* NSPINV=0 
* NSPINU=0 
* NSPINVBAR=0 
* NSPINUBAR=0 
* OPTIMITERATIONS=1000 
* OPTIMLVL=4 -M -l -C numerator_0.log
* INDSHIFT
#define NINITIALMOMENTA "2" 
#define NFINALMOMENTA "2" 
#define NPOL "2"
#define NCPOL "0"
#define NSPINV "0"
#define NSPINU "0"
#define NSPINVBAR "0" 
#define NLOOPMOMENTA "1"
#define NSPINUBAR "0"
#define INDSHIFT "4"
#define SGID "0"


* defintions for the computation
V p1,...,p40,k1,...,k40,c1,...,c40; * force this internal ordering in FORM
Auto V p,k,eps,ceps, sV,sVbar,sU,sUbar;
Auto S lm;

AutoDeclare Index s=4, mu=4, ind;
Set dirac: s1,...,s80;
Set diracdummy: sd1,...,sd80;
Set lorentz: mu1,...,mu80;
Set lorentzdummy: mud1,...,mud80;
Auto S x;
* for the integrand
CF hermconjugate, pol, cpol, uSpinor, ubarSpinor, vSpinor, vbarSpinor, sp(s);
CF gamma, gam, spinorU, spinorUbar,spinorV,spinorVbar,lVec, gMetric(s),deltaS(s),muL,indS;
CF spinor;
CT gammatensor;
CF GGstring, NN, gammatrace(c);

CF ampDenom;

* for replacements
CF LmbToCmbSubs;

* for components
CF penergy, spatial(s),spatialComp,spatialCompTmp;

*for clTD
CF sprop, prop, sqrt, penergy, spatial, onshellenergy, toCLTD;
Auto S E, m, energyk, y, sqrt;

S ii,aa,bb,cc,dd,m,n,y,z;

S yT,mT;


#define NCUTMOMENTA "`NFINALMOMENTA'+`NLOOPMOMENTA'"


* code for the computation

* Load the diagrams
*#include- input_`SGID'.h
L F = sprop(-k1, mT)*sprop(-k1 + p1, mT)*sprop(-k1 - p2, mT)*sprop(-k1 + p1 - k2, mT)*sprop(-k1 - p2 + k3, mT)*LmbToCmbSubs(k2,c1)*LmbToCmbSubs(k3,c2)*LmbToCmbSubs(k1,c3)*(1/4*pol(1,muL(-1))*pol(2,muL(-3))*ii^10*(gam(indS(1),lVec(-k1),indS(2))+deltaS(indS(1),indS(2))*mT)*(gam(indS(3),lVec(-k1+p1),indS(4))+deltaS(indS(3),indS(4))*mT)*(gam(indS(5),lVec(-k1-p2),indS(6))+deltaS(indS(5),indS(6))*mT)*(gam(indS(7),lVec(-k1+p1-k2),indS(8))+deltaS(indS(7),indS(8))*mT)*(gam(indS(9),lVec(-k1-p2+k3),indS(10))+deltaS(indS(9),indS(10))*mT)*1^2*gam(indS(4),muL(-1),indS(1))*gam(indS(2),muL(-3),indS(5))*yT^3*2^(1/2)*deltaS(indS(8),indS(3))*deltaS(indS(6),indS(9))*deltaS(indS(10),indS(7)))*(hermconjugate(1));




*********************************************************************************
*****************       PART I: INPUT CONVERSION       **************************
*********************************************************************************


* We only allow for the numerator of the effective vertex to be 1 :> no squaring/interferences
id hermconjugate(1) = 1;
if ( count(hermconjugate,1)>0 ); 
    Print "Only hermconjugat(1) is allowed: %t";
    exit "Critical error";
endif;
id denom_(?aa) = ampDenom(?aa);
if ( count(ampDenom,1)>0 ); 
    Print "I have not yet implemented denominators properly: %t";
    exit "Critical error";
endif;
* apply mapping of external and loop-momenta to the cmb
* from form manual  "one should not use more than a single one at the same time inside a term"
#do i=0,1
    id once LmbToCmbSubs(?aa) = replace_(?aa);
    if ( count(LmbToCmbSubs,1)>0 ) redefine i "0";
    .sort        
#enddo
.sort:lmb-to-cmb;

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
        id spatial(p1?,p?!vector_) = p(mu2)*spatial(mu3,p1);
        id spatial(p?!vector_,k?!vector_) = p(mu3)*k(mu4)*spatial(mu3,mu4);
        id p?(mu?) * spatial(mu?,p1?) = spatial(p,p1);
        id p?(mu1?)*k?(mu2?)* spatial(mu1?,mu2?) = spatial(p,k);
    endargument;
endargument;
.sort 


* this is for maple input:
id ii = i_;
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
id sp(p1?,p2?) = p1.p2;
B+ gam;
.sort




*********************************************************************************
***************** PART II: GAMMA-CHAINS AND TRACES ******************************
*********************************************************************************

* chain out momenta gam(...,lVec(p1+p2),...) = gam(...,p1,...)+gam(...,p2,...)
Keep brackets;
repeat id gam(?aa,lVec(p?vector_),?bb) = gam(?aa,p,?bb);
repeat; 
    id once gam(?aa,lVec(p?!vector_),?bb) = p(mu)*gam(?aa,mu,?bb);
    id p?(mu?)*gam(?aa,mu?,?bb) = gam(?aa,p,?bb);
endrepeat;
repeat id gam(?aa,s1?)*gam(s1?,?bb) = gam(?aa,?bb);
id gam(s1?,?aa,s1?) = gammatrace(?aa);

B+ gammatrace;
.sort:chain-out-momenta;

*perform traces
Keep brackets;
repeat;
    id once gammatrace(?aa) = g_(1,?aa);
    trace4,1;
endrepeat;
B+ gam;
.sort:gamma-traces;



* simplify gamma-strings
Keep brackets;
* these are all the chains which are not simplified
id gam(?aa) = gammatensor(?aa);
#call Gstring4D(gammatensor,0)
id gammatensor(?aa) = gamma(?aa);
.sort:simplify-gamma-chains;


*********** expand gamma chains **********************
CF gammaAll;
S iter, intSym ;

Table explSpinor(1:4,p?);
Fill explSpinor(1) = penergy(p);
Fill explSpinor(2) = spatialComp(p,1);
Fill explSpinor(3) = spatialComp(p,2);
Fill explSpinor(4) = spatialComp(p,3);

.sort
repeat id once gamma(?aa,mud1?,?bb)*gamma(?cc,mud1?,?dd) = sum_(y,0,3,gamma(?aa,y,?bb)*gamma(?cc,y,?dd));
repeat id once spinor(p?,sd?diracdummy)*gamma(sd?diracdummy,?aa) = sum_(y,1,4,spinor(p,y)*gamma(y,?aa));
repeat id once spinor(p?,sd?diracdummy)*gamma(?aa,sd?diracdummy) = sum_(y,1,4,spinor(p,y)*gamma(?aa,y));
repeat id once spinor(p?,s1?dirac)*gamma(s1?dirac,?aa) = sum_(y,1,4,spinor(p,y)*gamma(y,?aa));
repeat id once spinor(p?,s1?dirac)*gamma(?aa,s1?dirac) = sum_(y,1,4,spinor(p,y)*gamma(?aa,y));

id spinor(p?,intSym?int_) = explSpinor(intSym,p);
.sort

id gamma(?aa) = gammaAll(gamma(?aa));
B+ gammaAll;
.sort
nTable slash(1:4,1:4,p?);
nTable gamtab(1:4,1:4,0:3);

#include- definition_gamma_components.h
*********************** EXPAND GAMMA FUNCTIONS ************************************************************
keep brackets;

argument gammaAll;
    #do i=1,1
        label retry4;  
*       gamma(1,p1+q1,2) expansion              
        id p?(mu?)*gamma(xx?number_,mu?,aa?number_) = gamma(xx,p,aa);
        id once ifmatch->retry4 gamma(xx?number_,?aa,p?!vector_,aa?number_) = p(mu)*gamma(xx,?aa,mu,aa);
        id once ifmatch->retry4  gamma(x?int_, p?vector_,  y?int_) = slash(x,y,p);
        id once ifmatch->retry4  gamma(x?int_, intSym?int_,  y?int_) = gamma(x,y,intSym);    
        id once ifmatch->retry4  gamma(x?int_,?aa, p?vector_,  y?int_) = gamma(x,?aa,1)*slash(1,y,p)+gamma(x,?aa,2)*slash(2,y,p)+gamma(x,?aa,3)*slash(3,y,p)+gamma(x,?aa,4)*slash(4,y,p);
*        sum_(iter,1,4,gamma(x,?aa,iter)*slash(iter,y,p));
        id once ifmatch->retry4  gamma(x?int_,?aa, intSym?int_, y?int_) = sum_(iter,1,4,gamma(x,?aa,iter)*gamtab(iter,y,intSym));   
*       if ( count(gamma,1) ) redefine i "0";
    #enddo
endargument;
.sort
#do i=1,1
    id once gammaAll(aa?) =aa;
    if ( count(gammaAll,1) > 0 ) redefine i "0";
    .sort
#enddo

.sort:gamma-chain-explicit;
* sanity check
if (count(gamma,1) || count(gam,1));
        print "Some gammas are not replaced: %t";
        exit "Critical ERROR";
endif;




*********************************************************************************
***************** PART III: PREPARATION FOR CLTD ********************************
*********************************************************************************


* REPLACE EVERYTHING APART FROM LOOP-MOMENTA BY lm's
* replace scalar products involving loop-momenta
#do i=`NFINALMOMENTA'+1,`NCUTMOMENTA'
    id c`i'.p1? = penergy(c`i')*penergy(p1) - spatial(c`i',p1);
    argument;
        id c`i'.p1? = penergy(c`i')*penergy(p1) - spatial(c`i',p1);
    endargument;
#enddo
.sort:expand-sps-loop-mom;
***************** TRANSLATION OF FUNCTIONS INVOLVING ONLY MOMENTA (as in standard LTD)
* lm-replacements: we do not replace loop-energies, since they are needed for the PF-routine
* Convert the dot products and energies to a symbol
#$MAXK = `NCUTMOMENTA';
#$MAXP = `NINITIALMOMENTA';
#$OFFSET = 0;
#do i=1,`$MAXP'
    id penergy(p`i') = lm`$OFFSET';
    #$OFFSET = $OFFSET + 1;
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
        argument;
            id p`i'.p`j' = lm`$OFFSET';
            argument;
                id p`i'.p`j' = lm`$OFFSET';
            endargument;
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`i', p`j') = lm`$OFFSET';
        argument;
            id spatial(p`i', p`j') = lm`$OFFSET';
            argument;
                id spatial(p`i', p`j') = lm`$OFFSET';
            endargument;
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXK'
    #if (`i'<= `NFINALMOMENTA');
        id penergy(c`i') = lm`$OFFSET';
        argument;
            id penergy(c`i') = lm`$OFFSET';
        endargument;
    #endif
    #$OFFSET = $OFFSET + 1;
    #do j=1,`$MAXP'
* there sould not exist any, which involve loop-momenta, because sp's involving loop-momenta are expanded beforehand
        id c`i'.p`j' = lm`$OFFSET';
        argument;
            id c`i'.p`j' = lm`$OFFSET';
            argument;
                id c`i'.p`j' = lm`$OFFSET';
            endargument;            
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(p`j', c`i') = lm`$OFFSET';
        argument;
            id spatial(p`j', c`i') = lm`$OFFSET';
            argument;
                id spatial(p`j', c`i') = lm`$OFFSET';
            endargument;
        endargument;
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=`i',`$MAXK'
        id c`i'.c`j' = lm`$OFFSET';
        argument;
            id c`i'.c`j' = lm`$OFFSET';
            argument;
                id c`i'.c`j' = lm`$OFFSET';
            endargument;
        endargument;
        #$OFFSET = $OFFSET + 1;
        id spatial(c`i', c`j') = lm`$OFFSET';
        argument;
            id spatial(c`i', c`j') = lm`$OFFSET';
            argument;
                id spatial(c`i', c`j') = lm`$OFFSET';
            endargument;
        endargument;
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
argument prop;
    id sqrt(x?) = (x)^(1/2);
endargument;
.sort:introduce-lms;

*** replace on-shell energies by energy symbols (E) for the cltd-code
*** TODO: check that there are no E variables ***
ExtraSymbols, underscore, E;
#redefine oldextrasymbols "`extrasymbols_'"
argtoextrasymbol prop, 2;
id onshellenergy(y?) = y;
.sort:onshell-enegies;
multiply replace_(<E{`oldextrasymbols'+1}_,E{`oldextrasymbols'+1}>\
	                  ,...,<E`extrasymbols_'_,E`extrasymbols_'>);
.sort
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out`SGID'.txt> "\tE{`i'-`oldextrasymbols'} = %$;" $y
#enddo
#write<out`SGID'.txt> ""
B+ prop;
.sort:replace-energies;
**** REPLACE cut-energies by cLTD loop energy symbols (energyk`i')
#define LOOPS "`NLOOPMOMENTA'"
#define LoopSymbol "c";
#define LoopSymbolStart "`NFINALMOMENTA'+1";
Keep brackets;
#do i=0,`NLOOPMOMENTA'-1
    id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = prop(energyk`i');
    argument prop;
        id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = energyk`i';
    endargument;
#enddo
*notice: we also put energies in prop
B+ prop;
.sort-prep-cltd-final;


*********************************************************************************
*****************       PART IV: CLTD            ********************************
*********************************************************************************
Keep brackets;
* replace back energies
id prop(x?) = x;
.sort;
Keep brackets;
*call cLTD
#include ./pf_complete.frm 

ExtraSymbols, underscore, Z;
Format C;
Format O1,stats=on, saIter=1000;
#Optimize F
.sort
#write<out`SGID'.txt> "%O"
#write<out`SGID'.txt> "\treturn %e" F;
.sort
.end