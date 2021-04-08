#-
Auto S n,y,E,energyp,energyk;
CF ncmd, num, Eres, Echain;
CF den, prop, norm, error;
CF polex, polexbar;
CF a;

set onshellenergies:E0,...,E1000;

*--#[ partial_fractioning:
multiply norm(1)*num();
repeat id norm(E?)*prop(y0?,y1?) = norm(2*E*y1)*(den(y0-y1)-den(y0+y1)); 
.sort

* NOTE:
*    This procedure follows what presented in the paper 2009.05509
*    The only difference is that the numerator function is that the
*    divided difference is definded with a minus sing.
*    As a result we have that each time a ncmd gain an extra term 
*    we need to compensate with an additional minus sing.
#do i=0,{`LOOPS'-1}
* Normalise by the contour orientation 
    id norm(n?) = norm(-n);
  
* Activate numerator
    id num(?y) = num(1,?y,ncmd());

* Identify poles in the loop momentum energy 
    splitarg (energyk`i') den;
    id den(?a,energyp?) = den(energyp,?a);
    factarg (-1) den 1;
    id den(energyp?,y?,energyp1?) = polex(energyp1/y,0,0)/y;
    id den(energyp?,y?) = den(y*energyp);
    .sort:get-poles;

* Split the poles by their position in the complex plane
    splitarg polex; 
    .sort
    id polex(?y1,+E?onshellenergies,?y2,E1?,E2?) = polex(?y1,?y2,E1+E,E2);
    id polex(?y1,-E?onshellenergies,?y2,E1?,E2?) = polex(?y1,?y2,E1,E2+E);

* check for errors
    id polex(?y,0,0) = error(polex(?y,0,0));
    id polex(?y,E1?!{,0},E2?!{,0}) = error(polex(?y, E1,E2));
    if (count(error,1));
         print "[energyk`i'] Cannot find energy signature %t";
         exit "Critical ERROR";
    endif;

* re-absorb 
    id polex(?y,0,E?!{,0}) = polex(?y,-E);
    id polex(?y,E?!{,0},0) = polexbar(?y,E);
    transform, polex,polexbar , addargs(1,last);
    id polex(y?) = polex(-y);
    id polexbar(y?) = polexbar(-y);
*    print +s;

* Remove elements with zero residue
    if (count(polex,1) == 0) discard;
    if ( count(polexbar,1) == 0 ) 
        repeat id num(?n,ncmd(?n1))*polex(y1?) = -num(?n,ncmd(?n1,y1));
    .sort:x-xbar-`i';

    
* Start unfolding the expression 
    repeat;
        if ( count(polex,1) == 1 );
	    repeat id polex(y1?)*polexbar(y2?) = polex(y1)*den(y1-y2);
	    id num(0,?n,ncmd(?n1))*polex(y1?) = num(?n,ncmd(?n1));
	    id num(1,?n,ncmd(?n1))*polex(y1?) = -num(?n,ncmd(?n1,y1));
        endif;


        id once num(0,?n,ncmd(?n1))*polex(y1?)*polexbar(y2?) = 
                     -num(0,?n,ncmd(?n1))*Eres(y1,y2)*polexbar(y2);
    
        id once num(1,?n,ncmd(?n1))*polex(y1?)*polexbar(y2?) = 
                     +num(0,?n,ncmd(?n1,y1))*Eres(y1,y2)*polexbar(y2)
                     -num(1,?n,ncmd(?n1,y1))*polexbar(y2);

        repeat id Eres(y1?,y2?)*polexbar(y2?)*polexbar(y3?)= 
                      den(y1-y2)*polexbar(y3)*(polexbar(y2)+Eres(y1,y3));
        id once Eres(y1?,y2?)*polexbar(y2?)= den(y1-y2)*polexbar(y2);
    endrepeat;
    id num(n0?{0,1},?n) = num(?n);
    .sort:recursive-expansion-`i';

* Sanity check
    if (count(polexbar,1) || count(polex,1));
         print "[energyk`i'] Some roots were not absorbed: %t";
         exit "Critical ERROR";
    endif;

*    repeat id Echain(y0?, y1?,?y)= -den(y0-y1)*Echain(y0,?y);
*    repeat id Echain(y0?, y1?,?y)= den(y0-y1)*Echain(y0,?y);
*    id Echain(y0?)= 1;
*    .sort:terminate-k`i';
#enddo
*--#] partial_fractioning:

*id norm(?y) = 1;
*multiply k0^3;
*id num(ncmd(y1?,y2?,?y)) = 0;
*print +s ;
.sort 

* Read out the numerator instructions
#procedure NumUnfold
Auto S z,r,s,invden;
* Start unfoding the numerator instructions loop by loop
#do i = 0,{`LOOPS'-1}
	id once num(ncmd(?z),?y) =ncmd(?z,0)*num(?y);
	B+ ncmd, energyk`i';
	.sort:collect-ncmd;
	keep brackets;
	
	id energyk`i'^r?*ncmd(?z,z1?,0) = sum_(s,0, r - nargs_(?z) + 0, a(r,s,?z)*z1^s);
	B+ a;
	.sort:collect-a;
	keep brackets;

	repeat id a(r?,s?,?z,z1?) = -sum_(energyk, s+1,r-nargs_(?z) + 0, a(r,energyk,?z)*z1^(energyk-s-1));
	id a(r?,s?) = delta_(r,s);
	.sort:energy-k`i';
#enddo
id num=1;

* check that the substitution is complete
if (count(ncmd, 1,num,1,a, 1));
    Print "Unsubstituted ncmd: %t";
    exit "Critical error";
endif;
.sort:rm-num;
#endprocedure

* Expand numerator
#message "Start expanding the numerator"
#redefine oldextrasymbols "`extrasymbols_'"
ExtraSymbols, underscore, den;
argtoextrasymbol den;
id den(y?) = y;

.sort:den;
off statistics;
#call NumUnfold;
on statistics; 

*multiply replace_(<den{`oldextrasymbols'+1}_,den(extrasymbol_({`oldextrasymbols'+1}))>\
*                  ,...,<den`extrasymbols_'_,den(extrasymbol_(`extrasymbols_'))>);
multiply replace_(<den{`oldextrasymbols'+1}_,invden{1}>
                  ,...,<den`extrasymbols_'_,invden{`extrasymbols_'-`oldextrasymbols'}>);
#do i={`oldextrasymbols'+1},`extrasymbols_'
    #$y = extrasymbol_(`i');
    #write<out.txt> "\tinvden{`i'-`oldextrasymbols'} = %$;" $y
#enddo

.sort:end-numerator;

**Format 255;
*#clearoptimize;
*ExtraSymbols, underscore, Z;
*Format C;
*Format O1,stats=on, saIter=1000;
*#Optimize F
*#write "%O"
*.sort
*print +s;
*.sort
*#write<pf_form.out> "      F = %e",F

**print +s;
*.end
*

