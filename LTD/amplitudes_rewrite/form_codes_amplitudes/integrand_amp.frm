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
#define NCUTMOMENTA "`NFINALMOMENTA'+`NLOOPMOMENTA'"
#define OPTIMISATIONSTRATEGY "CSEgreedy"
#define OPTIMLVL "4"
#define OPTIMITERATIONS "1000"

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
CF ampDenom;

* for components
CF penergy, spatial(s), spatialComp, spatialCompTmp;

* for gamma-algebra
CF GGstring, NN, gammatrace(c);
CT gammatensor;

* for replacements
CF LmbToCmbSubs;


*for clTD
CF sprop, prop, sqrt, penergy, onshellenergy, toCLTD;
Auto S E, m, energyk, y, sqrt;

S ii,aa,bb,cc,dd,m,n,y,z;

S yT,mT;









* code for the computation

* Load the diagrams
*#include- input_`SGID'.h
L F = 1/(sp(k3,p1)+sp(k2,p1))^3*sprop(-k1, mT)*sprop(-k1 + p1, mT)*sprop(-k1 - p2, mT)*sprop(-k1 + p1 - k2, mT)*sprop(-k1 - p2 + k3, mT)*LmbToCmbSubs(k2,c1)*LmbToCmbSubs(k3,c2)*LmbToCmbSubs(k1,c3)*(1/4*pol(1,muL(-1))*pol(2,muL(-3))*ii^10*(gam(indS(1),lVec(-k1),indS(2))+deltaS(indS(1),indS(2))*mT)*(gam(indS(3),lVec(-k1+p1),indS(4))+deltaS(indS(3),indS(4))*mT)*(gam(indS(5),lVec(-k1-p2),indS(6))+deltaS(indS(5),indS(6))*mT)*(gam(indS(7),lVec(-k1+p1-k2),indS(8))+deltaS(indS(7),indS(8))*mT)*(gam(indS(9),lVec(-k1-p2+k3),indS(10))+deltaS(indS(9),indS(10))*mT)*1^2*gam(indS(4),muL(-1),indS(1))*gam(indS(2),muL(-3),indS(5))*yT^3*2^(1/2)*deltaS(indS(8),indS(3))*deltaS(indS(6),indS(9))*deltaS(indS(10),indS(7)))*(hermconjugate(1));
.sort



*********************************************************************************
*****************       PART I: INPUT CONVERSION       **************************
*********************************************************************************
* all the procedures in PART I are defined in 
#include input_conversion.h

* this is for maple input:
id ii = i_;
* We only allow for the numerator of the effective vertex to be 1 :> no squaring/interferences
id hermconjugate(1) = 1;
if ( count(hermconjugate,1)>0 ); 
    Print "Only hermconjugat(1) is allowed: %t";
    exit "Critical error";
endif;
* overall denominators of external kinematics
id denom_(?aa) = ampDenom(?aa);
if ( count(ampDenom,1)==0 ); 
    #define TREATAMDENOM "0"
else;
    #define TREATAMDENOM "1"
endif;

* apply mapping of external and loop-momenta to the cmb
* from form manual  "one should not use more than a single one at the same time inside a term"
#do i=0,1
    id once LmbToCmbSubs(?aa) = replace_(?aa);
    if ( count(LmbToCmbSubs,1)>0 ) redefine i "0";
    .sort        
#enddo

.sort:lmb-to-cmb;
#call scalar-prop-to-clTD-prop
.sort:bilinear-sp; 
#call translate-inds-and-ext-data




*********************************************************************************
***************** PART II: GAMMA-CHAINS AND TRACES ******************************
*********************************************************************************

* load relevant procedures
#include gamma_treatment.h
B+ gam;
.sort:translate-fcts;
Keep brackets;

#call gamma-traces-chains
.sort:simplify-gamma-chains;
#call expand-gamma-chains
* I have to ask ben why that is not done by default
symmetrize spatial;
.sort:gamma-chain-final;



*********************************************************************************
***************** PART III: PREPARATION FOR CLTD ********************************
*********************************************************************************

* we replace everything apart from loop-momenta energies by lms (variables for c-code)
* only loop-energies are needed for cLTD
* I have to ask ben why that is not done by default
* load relevant procedures
#include translate_to_lms.h
.sort
#call introduce-lms
.sort:introduce-lm;
print;
.sort
* treatement of overall denominators: optimize and export
#if `TREATAMDENOM'
    #redefine oldextrasymbols "`extrasymbols_'"
    argtoextrasymbol tonumber ampDenom;
    .sort:define-ampDenom;
    
    #do i={`oldextrasymbols'+1},`extrasymbols_'
        L globalDenom`i' = extrasymbol_(`i');
    #enddo
    .sort
    Hide F;
    print;
    .sort
    ExtraSymbols, underscore, Z;
    #do i={`oldextrasymbols'+1},`extrasymbols_'    
        Format C;
        Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
        #Optimize globalDenom`i'
        .sort
        #write<out`SGID'.txt> "ampDenom`i'\n"
        #write<out`SGID'.txt> "%O"
        #write<out`SGID'.txt> "\treturn %e" globalDenom`i';
        #write<out`SGID'.txt> "\n"
        #clearoptimize
        Drop globalDenom`i';
        .sort
    #enddo
    delete  extrasymbols>`oldextrasymbols';
    .sort
    Unhide F;
#endif
.sort

*** replace on-shell energies by energy symbols (E) for the cltd-code
ExtraSymbols, underscore, E;
#redefine oldextrasymbols "`extrasymbols_'"
argtoextrasymbol prop, 2;
id onshellenergy(y?) = y;
.sort:onshell-enegies;
multiply replace_(<E{`oldextrasymbols'+1}_,E{`oldextrasymbols'+1}>\
	                  ,...,<E`extrasymbols_'_,E`extrasymbols_'>);
.sort
#write<out`SGID'.txt> "energies \n"
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out`SGID'.txt> "\tE{`i'-`oldextrasymbols'} = %$;" $y
#enddo
#write<out`SGID'.txt> ""
delete  extrasymbols>`oldextrasymbols';
.sort:replace-energies;

**** REPLACE cut-energies by cLTD loop energy symbols (energyk`i')
#define LOOPS "`NLOOPMOMENTA'"
#define LoopSymbol "c";
#define LoopSymbolStart "`NFINALMOMENTA'+1";

#do i=0,`NLOOPMOMENTA'-1
    id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = prop(energyk`i');
    argument prop;
        id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = energyk`i';
    endargument;
#enddo
*notice: we also put energies in prop
B+ prop;
.sort-prep-cltd-final;

CF diag;
Keep brackets;

multiply diag(1);
repeat id prop(?aa)*diag(bb?) = diag(prop(?aa)*bb);
argument diag;
* replace back energies
    id prop(x?) = x;
endargument;
#redefine oldextrasymbols "`extrasymbols_'"
argtoextrasymbol tonumber diag;
.sort:diag-replacements;
#do i={`oldextrasymbols'+1},`extrasymbols_'
    L diag`i' = extrasymbol_(`i');
#enddo
.sort
#define NDIAGS "`extrasymbols_'"
delete  extrasymbols>`oldextrasymbols';
Hide F;
.sort

*********************************************************************************
*****************       PART IV: CLTD            ********************************
*********************************************************************************
#redefine oldextrasymbols "`extrasymbols_'"
#include pf_complete.h 
Keep brackets;

* Get cLTD expression
#call partial-fractioning

* Expand numerator
ExtraSymbols, underscore, den;
argtoextrasymbol den;
id den(y?) = y;
.sort:den;
off statistics;
#call unfold-numerator;
on statistics; 

* Store inverse denominators
multiply replace_(<den{`oldextrasymbols'+1}_,invden{1}>
                  ,...,<den`extrasymbols_'_,invden{`extrasymbols_'-`oldextrasymbols'}>);
#write<out`SGID'.txt> "denoms \n"
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out`SGID'.txt> "\tinvden{`i'-`oldextrasymbols'} = 1/(%$);" $y
#enddo
.sort:end-numerator;
delete  extrasymbols>`oldextrasymbols';
.sort
*********************************************************************************
*****************       PART V: EXPORT            *******************************
*********************************************************************************


* Optimize and export the diags1,...,diagsN
ExtraSymbols, underscore, Z;
#do i=1,`NDIAGS'    
    Format C;
    Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
    #Optimize diag`i'
    .sort
    #write<out`SGID'.txt> "diag`i'\n"
    #write<out`SGID'.txt> "%O"
    #write<out`SGID'.txt> "\treturn %e" diag`i';
    #write<out`SGID'.txt> "\n"
    #clearoptimize
    Drop diag`i';
    .sort
#enddo
.sort
Unhide F;
.sort
* Optimize and export complete expression
ExtraSymbols, underscore, Z;
Format C;
Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
#Optimize F
.sort
#write<out`SGID'.txt> "integrand\n"
#write<out`SGID'.txt> "%O"
#write<out`SGID'.txt> "\treturn %e" F;
.end