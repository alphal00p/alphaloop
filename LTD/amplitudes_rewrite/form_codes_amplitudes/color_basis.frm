#-
Off statistics;

*--#[ setup :
#-
Off statistics;

* I expect the following input:
*
*#define INDSHIFT "4" (some int>0)
*#define SGID "0" (some int>0)
*#define CPERGRAPH "0" (1 or 0)

* defintions for the computation
V p1,...,p40,k1,...,k40,c1,...,c40; * force this internal ordering in FORM
Auto V p,k,eps,ceps, sV,sVbar,sU,sUbar,c;
Auto S x;

* for the integrand
CF sprop,hermconjugate, pol, cpol, uSpinor, ubarSpinor, vSpinor, vbarSpinor, sp(s);
CF gamma, gam, spinorU, spinorUbar,spinorV,spinorVbar,lVec, gMetric(s),deltaS(s),muL,indS;
CF spinor;

#include- coltracebased.h

I ind1, ind2;
CF  deltaA, deltaF, Tc , colA,colF, colorfct, fC, ffC,colTemp;
CF trCol;
S a,x,m,n,Ca,Cf;


S ii,aa,bb,cc,dd,m,n,y,z;


* Load the diagrams
#include- input_`SGID'.h
id sprop(?aa) = 1;
.sort

* Transform indices
Argument; 
    id colA(aa?) = colA(aa+`INDSHIFT');
    id colF(aa?)  = colF(aa+`INDSHIFT');
EndArgument;
Argument;
    id colA(aa?) = ccolF[aa];
    id colF(aa?)  = ccolA[aa];
EndArgument;
id deltaA(cOlj1?,cOlj2?) = d_(cOlj1,cOlj2);
id deltaF(cOli1?,cOli2?) = d_(cOli1,cOli2);
id fC(cOlj1?,cOlj2?,cOlj3?) = fCol(cOlj1,cOlj2,cOlj3);
id ffC(?a) = ffCol(?a);
id Tc(?a) = TCol(?a);
.sort:translate-inds-and-fcts;


******************
* Color evaluation
******************
#call tracebased
* set the SU(3) values
Multiply replace_(cOlNF,3,cOlNA,8,Tf,1/2);
.sort:color;
*Manipulate color string
Splitarg,color;
repeat id color(x?,y?,?a)= color(x)+color(y,?a);
FactArg,color;
repeat id color(x?,y?,?a)= color(x)*color(y,?a);
id color(x?number_) = x;
.sort:color-final;

* Factor out rational numbers and i
repeat id  color(y?*x?) = color(x)*color(y);
id color(x?number_) = x;
id color(i_) = i_;
.sort:factor-out-rationals;
* in case there is no color structure
Multiply color(1);
repeat id  color(x?)*color(y?) = color(x*y);
*id p1?.p2? = sp(p1,p2);
.sort
* factorize color per graph
#if `CPERGRAPH';
    AB+ color;
#endif;
.sort
#if `CPERGRAPH';
    collect numtemp;
    FactArg numtemp;
    repeat id numtemp(?x,y?number_,?z) = y*numtemp(?x,?z);
    argument numtemp;
        id color(x?) = x;
    endargument;
    repeat id numtemp(?x) = color(?x);
#endif;
.sort
B+ color;
.sort


* split off every color configuration into a new expression
id color(?a) = color(color(?a));
argtoextrasymbol tonumber,color,1;
#redefine oldextrasymbols "`extrasymbols_'"
B+ color;
.sort:color-1;
Hide F;
#redefine colorsymbolstart "`extrasymbols_'"
#write<SG_`SGID'_color_decomp.txt> "{\n"
#do cc={`oldextrasymbols'+1}, `colorsymbolstart'
    delete extrasymbols>`colorsymbolstart'; * clear all extra symbols from the last configuration
    #$color = extrasymbol_(`cc');
    #write<SG_`SGID'_color_decomp.txt> "\"%$\":", $color;
    #$struc = F[color(`cc')];
    #write<SG_`SGID'_color_decomp.txt> "\"%$\",\n", $struc;
#enddo
    #write<SG_`SGID'_color_decomp.txt> "}"
.end