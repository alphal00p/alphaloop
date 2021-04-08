#-
#include ./pf_complete.frm 
CF sprop, prop, sqrt, penergy, spatial, onshellenergy;
Auto S E, m, energyk, y;
Auto V k,p, q, mom, armin;
Auto I mu;

L F = penergy(k1)*sprop(k1+q1,m1)*sprop(k1+k2+q5,m5)*sprop(k2+p4,m4);
print +s;

id sprop(mom?,y2?) = prop(penergy(mom), sqrt(mom.mom + y2^2));
argument prop;
    id penergy(mom?!vector_) = mom(mu)*penergy(mu);
    id mom?(mu) * penergy(mu) = penergy(mom);
endargument;
.sort 

* Extract energies
#redefine oldextrasymbols "`extrasymbols_'"
ExtraSymbols, underscore, E;
argtoextrasymbol prop, 2;
id onshellenergy(y?) = y;
.sort:onshell-enegies;
multiply replace_(<E{`oldextrasymbols'+1}_,E{`oldextrasymbols'+1}>\
	                  ,...,<E`extrasymbols_'_,E`extrasymbols_'>);
.sort
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out.txt> "\tE{`i'-`oldextrasymbols'} = %$;" $y
#enddo
#write<out.txt> ""


* Map loop momenta energies
#define LOOPS "2";
#define LoopSymbol "k";
#define LoopSymbolStart "1";
#do i=0,4
    id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = energyk`i';
    argument prop;
        id penergy(`LoopSymbol'{`LoopSymbolStart'+`i'}) = energyk`i';
    endargument;
#enddo
.sort


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
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out.txt> "\tinvden{`i'-`oldextrasymbols'} = 1/(%$);" $y
#enddo
.sort:end-numerator;

* Optimize
ExtraSymbols, underscore, Z;
Format C;
Format O1,stats=on, saIter=1000;
#Optimize F
.sort
#write<out.txt> "%O"
#write<out.txt> "\treturn %e" F;

.sort
*#write "%X"
.end
