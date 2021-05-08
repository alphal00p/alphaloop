#-
Off statistics;

* I expect the following input:
*
*#define NINITIALMOMENTA "2" (number of INDEPENDENT outgoing momenta)
*#define NFINALMOMENTA "2" 
*#define NPOL "2"
*#define NCPOL "0"
*#define NSPINV "0"
*#define NSPINU "0"
*#define NSPINVBAR "0" 
*#define NLOOPMOMENTA "1"
*#define NSPINUBAR "0"
*#define INDSHIFT "4"
*#define SGID "0"
*#define NCUTMOMENTA "`NFINALMOMENTA'+`NLOOPMOMENTA'"
*#define OPTIMISATIONSTRATEGY "CSEgreedy"
*#define OPTIMLVL "4"
*#define OPTIMITERATIONS "1000"
*#define DEBUGLVL "1"
*#define NUMERICGCHAINS "1"
*#define NUMERICGCHAINS "1"
* defintions for the computation
V p1,...,p40,k1,...,k40,c1,...,c40; * force this internal ordering in FORM
Auto V p,k,eps,ceps, sV,sVbar,sU,sUbar,c;
Auto S lm;

AutoDeclare Index s=4, mu=4, ind;
Set dirac: s1,...,s80;
Set diracdummy: sd1,...,sd80;
Set lorentz: mu1,...,mu80;
Set lorentzdummy: mud1,...,mud80;
Auto S x;
CF repl;
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




*for clTD
CF sprop, prop, sqrt, penergy, onshellenergy, toCLTD;
Auto S E, m, energyk, y, sqrt;

S ii,aa,bb,cc,dd,m,n,y,z;




* code for the computation
#define NCUTMOMENTA "`NFINALMOMENTA'+`NLOOPMOMENTA'"
* Load the diagrams
#include- input_`SGID'.h

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
id ampDenom(1) = 1;
.sort

* treat overall numerator
#define TREATAMDENOM "0"
if ( count(ampDenom,1) ); 
    #redefine TREATAMDENOM "1"
endif;
.sort

* apply mapping of external and loop-momenta to the cmb
#do i=1,`NFINALMOMENTA'
    multiply replace_(p{`i'+`NINITIALMOMENTA'},c`i');
    .sort        
#enddo
#do i=1,`NLOOPMOMENTA'
    multiply replace_(k{`i'},c{`i'+`NFINALMOMENTA'});
    multiply (2*pi_*i_);
    .sort        
#enddo
.sort:lmb-to-cmb;
Format Mathematica;
#if (`DEBUGLVL'>1)
    #write<debug_stage_0_diag_`SGID'.m> "{diag`SGID'->%E}" F
#endif
.sort-debug-output;

#call scalar-prop-to-clTD-prop
.sort:bilinear-sp; 
#call translate-inds-and-ext-data
.sort
Format Mathematica;
#if (`DEBUGLVL'>1)
    #write<debug_stage_1_diag_`SGID'.m> "{diag`SGID'->%E}" F
#endif
.sort



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
#if (`NUMERICGCHAINS'<=0)
    #call expand-gamma-chains
#endif
symmetrize spatial;
.sort:gamma-chain-final;


#if (`NUMERICGCHAINS'>0)
    V envec, spatvec1, spatvec2 spatvec3;
    CF gchain;
    Auto S spatialc;

    id spinor(p?,sd1?)*gamma(sd1?,?aa) =gamma(p,?aa);
    id spinor(p?,sd1?)*gamma(?aa,sd1?) = gamma(?aa,p);
    id spinor(p?,s1?)*gamma(s1?,?aa) =gamma(p,?aa);
    id spinor(p?,s1?)*gamma(?aa,s1?) = gamma(?aa,p);
    B+ gamma;
    .sort
    Keep brackets;
* split off energies only
    #do i=`NFINALMOMENTA'+1, `NFINALMOMENTA'+`NLOOPMOMENTA'
        repeat;
            id once gamma(?aa,c`i',?bb) = penergy(c`i')*gamma(?aa,envec,?bb) + gamma(?aa,spatialc`i',?bb); 
        endrepeat;
    #enddo
    
* if that becomes a problem the gamma chain has to be re-tought    
    repeat;
        id once gamma(?aa,mu?,?bb)*gamma(?cc,mu?,?dd) = gamma(?aa,envec,?bb)*gamma(?cc,envec,?dd)
            - gamma(?aa,spatvec1,?bb)*gamma(?cc,spatvec1,?dd)
            - gamma(?aa,spatvec2,?bb)*gamma(?cc,spatvec2,?dd)
            - gamma(?aa,spatvec3,?bb)*gamma(?cc,spatvec3,?dd);
    endrepeat;
* put into gamma-chain format
    id gamma(?aa) = gamma(nargs_(?aa),gchain(?aa));
    #redefine oldextrasymbols "`extrasymbols_'"
    argtoextrasymbol tonumber gamma 2;
    .sort
    #write<out_integrand_PF_`SGID'.proto_c> "// dummy vectors";
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s envec[4] = {1.,0.,0.,0.};";
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s spatvec1[4] = {0.,1.,0.,0.};";
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s spatvec2[4] = {0.,0.,1.,0.};";
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s spatvec3[4] = {0.,0.,0.,1.};\n";
#endif
symmetrize spatial;
.sort:gamma-chain-final;
.sort


* I have to ask ben why that is not done by default


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
repeat;
    id once repl(?aa) = replace_(?aa);
endrepeat;
id replace_(?aa) = 1;
.sort:repl;
* write out g-chains

#if (`NUMERICGCHAINS'>0)  
    Format C;
    #write<out_integrand_PF_`SGID'.proto_c> "//gamma chains\n"
    #do i={`oldextrasymbols'+1},`extrasymbols_'
        #$y = extrasymbol_(`i');
        #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s *gchain`i'[] = %$;" $y;
    #enddo
    delete  extrasymbols>`oldextrasymbols';
    
#endif
#write<out_integrand_PF_`SGID'.proto_c> "\n";

.sort
* declare each gamma_chain to a symbol
CF gamcollector, gc;
multiply gamcollector(1);
repeat id gamma(?aa)*gamcollector(?bb) = gamcollector(gamma(?aa),?bb);
id gamcollector(?aa,1) =gamcollector(?aa);
#redefine oldextrasymbols "`extrasymbols_'"
argtoextrasymbol tonumber gamcollector;
.sort
#if (`NUMERICGCHAINS'>0)  
    Format C;
    #write<out_integrand_PF_`SGID'.proto_c> "//gamma chain computation\n"
    #do i={`oldextrasymbols'+1},`extrasymbols_'
        #$y = extrasymbol_(`i');
        #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s gc`i'= %$;" $y;
    #enddo
    delete  extrasymbols>`oldextrasymbols';  
#endif
#write<out_integrand_PF_`SGID'.proto_c> "\n";
b+ gamcollector;
.sort:collect-gammas;
Keep brackets;
repeat id gamcollector(aa?int_,?cc) = gamcollector(gc(aa))*gamcollector(?cc);
id gamcollector =1;
repeat id gamcollector(aa?)*gamcollector(bb?) = gamcollector(aa*bb);

* factorize full chain
#redefine oldextrasymbols "`extrasymbols_'"
argtoextrasymbol tonumber gamcollector;
.sort
#if (`NUMERICGCHAINS'>0)  
    Format C;
    #write<out_integrand_PF_`SGID'.proto_c> "//gamma chain multiplications \n"
    #do i={`oldextrasymbols'+1},`extrasymbols_'
        #$y = extrasymbol_(`i');
        #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s gcmult`i'= %$;" $y;
    #enddo
    delete  extrasymbols>`oldextrasymbols';  
#endif
#write<out_integrand_PF_`SGID'.proto_c> "\n";
b+ gamcollector;

.sort;
CF gcmult;
Keep brackets;
id gamcollector(aa?int_) = gcmult(aa);
.sort:gamma-factorization-final;
* treatement of overall denominators: optimize and export
#if `TREATAMDENOM'
    

    #redefine oldextrasymbols "`extrasymbols_'"
    argtoextrasymbol tonumber ampDenom;
    .sort
    Hide F;
    .sort:define-ampDenom;
    
    #do i={`oldextrasymbols'+1},`extrasymbols_'
        L globalDenom`i' = extrasymbol_(`i');
    #enddo
    .sort
    #if (`DEBUGLVL'>1)
        #write<debug_overall_denoms_diag_`SGID'.m> "{"
    #endif

    #do i={`oldextrasymbols'+1},`extrasymbols_'
        #if (`DEBUGLVL'>1)
            Format Mathematica;
            #write<debug_overall_denoms_diag_`SGID'.m> "ampDenom[`i'] = 1/(%E)," globalDenom`i'
            #write<debug_overall_denoms_diag_`SGID'.m> "\n"
        #endif
        .sort-debug;   
        ExtraSymbols, underscore, AD`i'Z; 
        Format C;
        Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
        #Optimize globalDenom`i'
        .sort
        #write<out_integrand_PF_`SGID'.proto_c> "//ampDenom`i'\n"
        #write<out_integrand_PF_`SGID'.proto_c> "%O"
        #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s ampDenom`i' = 1/(%E);" globalDenom`i'
        #write<out_integrand_PF_`SGID'.proto_c> "\n"
        #clearoptimize
        Drop globalDenom`i';
        .sort
    #enddo
    delete  extrasymbols>`oldextrasymbols';
    .sort
    #if (`DEBUGLVL'>1)
        #write<debug_overall_denoms_diag_`SGID'.m> "1->1}"
    #endif
    Unhide F;
#endif
.sort
#if (`DEBUGLVL'>1)
    Format Mathematica;
    #write<debug_stage_2_diag_`SGID'.m> "{diag`SGID'->%E}" F
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


#write<out_integrand_PF_`SGID'.proto_c> "//energies \n"
#if (`DEBUGLVL'>1)
    #write<debug_energies_diag_`SGID'.m> "energies -> {"
#endif
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #if (`DEBUGLVL'>1)
            Format Mathematica;
            #write<debug_energies_diag_`SGID'.m> "E{`i'-`oldextrasymbols'} -> %$," $y
            #write<debug_energies_diag_`SGID'.m> "\n"
    #endif
    .sort:debug_output;
    Format C;
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s E{`i'-`oldextrasymbols'} = %$;" $y
#enddo
#write<out_integrand_PF_`SGID'.proto_c> ""
#if (`DEBUGLVL'>1)
    #write<debug_energies_diag_`SGID'.m> "1->1}"
#endif
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
.sort:replace_loop_energies;

* safety check
*at this point, there should be no vectors. If there are still vectors, it means the input is not valid
if (match(p?));
    Print "Vector without replacement: %t";
    exit "Critical error";
endif;

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
#write<out_integrand_PF_`SGID'.proto_c> "//denoms cLTD \n"
#if (`DEBUGLVL'>1)
    #write<debug_cltd_diag_`SGID'.m> "{"
#endif

#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #if (`DEBUGLVL'>1)
            Format Mathematica;
            #write<debug_cltd_diag_`SGID'.m> "invden{`i'-`oldextrasymbols'} -> 1/(%$)," $y
            #write<debug_cltd_diag_`SGID'.m> "\n"
    #endif
    .sort:debug_output;

    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s invden{`i'-`oldextrasymbols'} = 1/(%$);" $y
#enddo
.sort:end-numerator;
delete  extrasymbols>`oldextrasymbols';
.sort
*********************************************************************************
*****************       PART V: EXPORT            *******************************
*********************************************************************************


* Optimize and export the diags1,...,diagsN

#do i=1,`NDIAGS'
    #if (`DEBUGLVL'>1)
            Format Mathematica;
            #write<debug_cltd_diag_`SGID'.m> "diag[`i'] -> %E," diag`i'
            #write<debug_cltd_diag_`SGID'.m> "\n"
    #endif
    .sort:debug_output;    
    ExtraSymbols, underscore, DIA`i'Z;
    Format C;
    Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
    #Optimize diag`i'
    .sort
    #write<out_integrand_PF_`SGID'.proto_c> "\n//energy config `i'\n"
    #write<out_integrand_PF_`SGID'.proto_c> "%O"
    #write<out_integrand_PF_`SGID'.proto_c> "%%(numbertype)s diag`i' = %e" diag`i';
    #write<out_integrand_PF_`SGID'.proto_c> "\n"
    #clearoptimize
    Drop diag`i';
    .sort
#enddo
#if (`DEBUGLVL'>1)
    #write<debug_cltd_diag_`SGID'.m> "1 -> 1}" 
#endif
.sort
Unhide F;
.sort
#if (`DEBUGLVL'>1)
    Format Mathematica;
    #write<debug_stage_3_diag_`SGID'.m> "{diag`SGID'->%E}" F
#endif
.sort
* Optimize and export complete expression
ExtraSymbols, underscore, Z;
Format C;
Format O`OPTIMLVL',method=`OPTIMISATIONSTRATEGY',stats=on,saIter=`OPTIMITERATIONS';
#Optimize F
.sort
#write<out_integrand_PF_`SGID'.proto_c> "//integrand\n"
#write<out_integrand_PF_`SGID'.proto_c> "%O"
#write<out_integrand_PF_`SGID'.proto_c> "\t*out = %e" F;
#write<success_`SGID'.proto_c> "1";
.sort

#if (`DEBUGLVL'>1)
    Format Mathematica;
    #write<debug_full_diag_`SGID'.m> "(*This code debugs diagram `SGID'*)\n\n"
    #write<debug_full_diag_`SGID'.m> "(*useful definitions*)"
    #write<debug_full_diag_`SGID'.m> "sp[a_List,b_List]:=a[[1]]*b[[1]]-a[[2;;]].b[[2;;]];"
    #write<debug_full_diag_`SGID'.m> "spatial[a_List,b_List]:=a[[2;;]].b[[2;;]];"
    #write<debug_full_diag_`SGID'.m> "penergy[a_List]:=a[[1]];"
    #write<debug_full_diag_`SGID'.m> "spatialComp[a_List,b_Integer]:=a[[b+1]];"
    #write<debug_full_diag_`SGID'.m> "importFORM[x_String]:=ToExpression@(StringReplace[{\"\\n\"->\"\",\" \"->\"\",\"i_\"->\"I\",\"pi_\"->\"Pi\"}][Import[x,\"Text\"]]); \n\n"


    #write<debug_full_diag_`SGID'.m> "(* Replacement rules for dummy functions and variables *)\n"
    #write<debug_full_diag_`SGID'.m> "(* overall denominators wrapped in ampDenom wrapper *)"
    #write<debug_full_diag_`SGID'.m> "overallDenoms = importFORM[\"debug_overall_denoms_diag_`SGID'.m\"];"
    #write<debug_full_diag_`SGID'.m> "(* lms *)"
    #write<debug_full_diag_`SGID'.m> "lmRepl=Reverse/@importFORM[\"debug_lm_diag_`SGID'.m\"];"
    #write<debug_full_diag_`SGID'.m> "(* replacements from cLTD: energies *)"
    #write<debug_full_diag_`SGID'.m> "cLTDEnergies=importFORM[\"debug_energies_diag_`SGID'.m\"];"
    #write<debug_full_diag_`SGID'.m> "(* replacements from cLTD: denoms and configs in wrapper diag[] *)"
    #write<debug_full_diag_`SGID'.m> "cLTDRepl=importFORM[\"debug_cltd_diag_`SGID'.m\"];\n"
    



    #write<debug_full_diag_`SGID'.m> "(* diagram as imported mapped to the cmb*)"
    #write<debug_full_diag_`SGID'.m> "diagOringinal=diag0 /. importFORM[\"debug_stage_0_diag_`SGID'.m\"] ;\n"
    #write<debug_full_diag_`SGID'.m> "(* diagram mapped propagators and expanded SPs*)"
    #write<debug_full_diag_`SGID'.m> "diagMod1=diag0 /. importFORM[\"debug_stage_1_diag_`SGID'.m\"] ;\n"
    #write<debug_full_diag_`SGID'.m> "(* diagram gamma algebra performed and introduced lms*)"
    #write<debug_full_diag_`SGID'.m> "diagMod2=diag0 /. importFORM[\"debug_stage_2_diag_`SGID'.m\"] ;\n"
    #write<debug_full_diag_`SGID'.m> "(* diagram final version: CLTD is performed*)"
    #write<debug_full_diag_`SGID'.m> "diagMod3=diag0 /. importFORM[\"debug_stage_3_diag_`SGID'.m\"] ;\n"

#endif
.sort
.end