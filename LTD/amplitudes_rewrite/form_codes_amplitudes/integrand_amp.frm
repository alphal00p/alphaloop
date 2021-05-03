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