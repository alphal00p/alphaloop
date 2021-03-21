#-
Format nospaces;
Format 100; 

*--#[ partial_fractioning_vars :
Auto S n, y,E,p,k;
CF ncmd, num, Eres, Echain;
CF den, prop, norm, error;
CF x, xbar;

set energies:E0,...,E1000;
set shifts:p0,...,p1000;
*--#] partial_fractioning_vars :

* The propagators are defined has prop(q0,qE):=q0^2-qE^2;

off statistics;
**************
*** 1 LOOP ***
**************
*Bubble
L F2 = prop(k0+p0,E0)*prop(k0+p1,E1);
*Triangle
L F3 = prop(k0+p0,E0)*prop(k0+p1,E1)*prop(k0+p2,E2);
*Box
*L F4 = prop(k0+p0,E0)*prop(k0+p1,E1)*prop(k0+p2,E2)*prop(k0+p3,E3);
L F4 = prop(k0-p1,E1)*prop(k0-p2,E2)*prop(k0-p3,E3)*prop(k0-p4,E4);
*L F4 = prop(k0,E1)*prop(k0-p1,E2)*prop(k0-p1-p2,E3)*prop(k0+p1,E4);
*Box
*L F5 = prop(k0+p0,E0)*prop(k0+p1,E1)*prop(k0+p2,E2)*prop(k0+p3,E3)*prop(k0+p4,E4);
L F5 = prop(k0-p1,E1)*prop(k0-p2,E2)*prop(k0-p3,E3)*prop(k0-p4,E4)*prop(k0-p5,E5);
*Decagon
L F10 = prop(k0+p0,E0)*prop(k0+p1,E1)*prop(k0+p2,E2)*prop(k0+p3,E3)*
        prop(k0+p4,E4)*prop(k0+p5,E5)*prop(k0+p6,E6)*prop(k0+p7,E7)*
        prop(k0+p8,E8)*prop(k0+p9,E9);
*Endecagon
L F11 = prop(k0+p0,E0)*prop(k0+p1,E1)*prop(k0+p2,E2)*prop(k0+p3,E3)*
        prop(k0+p4,E4)*prop(k0+p5,E5)*prop(k0+p6,E6)*prop(k0+p7,E7)*
        prop(k0+p8,E8)*prop(k0+p9,E9)*prop(k0+p10,E10);

**************
*** 2 LOOP ***
**************
*Sunrise
L F2a = prop(k0+p0,E0)*prop(k0-k1+p1,E1)*prop(k1+p2,E2);
L F2b = prop(k0+p0,E0)*prop(-k0-k1+p1,E1)*prop(k1+p2,E2);
*PentaBox
L F2c = prop(k0+p0,E0)*
        prop(k0-k1+p1,E1)*prop(k0-k1+p2,E2)*prop(k0-k1+p3,E3)*
        prop(k1+p4,E4)*prop(k1+p5,E5)*prop(k1+p6,E6)*prop(k1+p7,E7);
*UV topology
L Ftest = prop(k0,E0)*prop(k0-p0,E1)*prop(k0,E2)^6*
          prop(k1,E3)*prop(k1,E4)^4*
          prop(k0+k1-p0,E5)*prop(k0+k1,E6)^5;

**************
*** 4 LOOP ***
**************
* Fishnet 2x2
L F4a =	prop(k0+p0,E0)*prop(k0+p1,E1)*
        prop(k1+p2,E2)*prop(k1+p3,E3)*
        prop(k0-k1+p4,E4)*
        prop(-k0-k2+p5,E5)*
        prop(-k1-k3+p6,E6)*
        prop(+k2+p7,E7)*prop(+k2+p8,E8)*
        prop(-k2+k3+p9,E9)*
        prop(-k3+p10,E10)*prop(-k3+p11,E11);
** Kite
L F4b =	prop(k0+p0,E0)*
        prop(k1+p1,E1)*
        prop(k0-k1+p2,E2)*
        prop(-k0-k2+p3,E3)*
        prop(-k1-k3+p4,E4)*
        prop(+k2+p5,E5)*
        prop(-k2+k3+p6,E6)*
        prop(-k3+p7,E7);
*** 
L F4c =	prop(k0+p0,E0)*prop(k0+p1,E1)*
        prop(k1+p2,E2)*
        prop(k0-k1+p3,E3)*
        prop(-k0-k2+p4,E4)*
        prop(-k1-k3+p5,E5)*
        prop(+k2+p6,E6)*
        prop(-k2+k3+p7,E7)*
        prop(-k3+p8,E8);
*** 
L F4d =	prop(k0+p0,E0)*prop(k0+p1,E1)*
        prop(k1+p2,E2)*prop(k1+p3,E3)*
        prop(k0-k1+p4,E4)*
        prop(-k0-k2+p5,E5)*
        prop(-k1-k3+p6,E6)*
        prop(+k2+p7,E7)*
        prop(-k2+k3+p8,E8)*
        prop(-k3+p9,E9);
*** 
L F4e =	prop(k0+p0,E0)*prop(k0+p1,E1)*
        prop(k1+p2,E2)*prop(k1+p3,E3)*
        prop(k0-k1+p4,E4)*
        prop(-k0-k2+p5,E5)*
        prop(-k1-k3+p6,E6)*
        prop(+k2+p7,E7)*prop(+k2+p8,E8)*
        prop(-k2+k3+p9,E9)*
        prop(-k3+p10,E10);
*** 
L F4f =	prop(k0+p0,E0)*prop(k0+p1,E1)*
        prop(k1+p2,E2)*prop(k1+p3,E3)*
        prop(k0-k1+p4,E4)*
        prop(-k0-k2+p5,E5)*
        prop(-k1-k3+p6,E6)*
        prop(+k2+p7,E7)*prop(+k2+p8,E8)*
        prop(-k2+k3+p9,E9)*
        prop(-k3+p10,E10)*prop(-k3+p11,E11);
.sort
on statistics;
hide;

* Select ONE

#ifndef `topoID'
	#define topoID "2";
#endif
#ifndef `LOOPS'
	#define LOOPS "1";
#endif
L F = F`topoID';
*L F = Ftest;
.sort
Drop;
ndrop F;

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
    splitarg (k`i') den;
    id den(?a,p?) = den(p,?a);
    factarg (-1) den 1;
    id den(p?,y?,p1?) = x(p1/y,0,0)/y;
    id den(p?,y?) = den(y*p);
    .sort:get-poles;

* Split the poles by their position in the complex plane
    splitarg x; 
    .sort
    id x(?y1,+E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1+E,E2);
    id x(?y1,-E?energies,?y2,E1?,E2?) = x(?y1,?y2,E1,E2+E);

* check for errors
    id x(?y,0,0) = error(x(?y,0,0));
    id x(?y,E1?!{,0},E2?!{,0}) = error(x(?y, E1,E2));
    if (count(error,1));
         print "[k`i'] Cannot find energy signature %t";
         exit "Critical ERROR";
    endif;

* re-absorb 
    id x(?y,0,E?!{,0}) = x(?y,-E);
    id x(?y,E?!{,0},0) = xbar(?y,E);
    transform, x,xbar , addargs(1,last);
    id x(y?) = x(-y);
    id xbar(y?) = xbar(-y);
*    print +s;

* Remove elements with zero residue
    if (count(x,1) == 0) discard;
    if ( count(xbar,1) == 0 ) 
        repeat id num(?n,ncmd(?n1))*x(y1?) = -num(?n,ncmd(?n1,y1));
    .sort:x-xbar;

    
* Start unfolding the expression 
    repeat;
        if ( count(x,1) == 1 );
*	    repeat id xbar(?y1)*xbar(?y2) = xbar(?y1,?y2);
*	    id num(0,?n,ncmd(?n1))*x(y1?)*xbar(?y2) = num(?n,ncmd(?n1))*Echain(y1,?y2);
*	    id num(1,?n,ncmd(?n1))*x(y1?)*xbar(?y2) = -num(?n,ncmd(?n1,y1))*Echain(y1,?y2);
	    repeat id x(y1?)*xbar(y2?) = x(y1)*den(y1-y2);
	    id num(0,?n,ncmd(?n1))*x(y1?) = num(?n,ncmd(?n1));
	    id num(1,?n,ncmd(?n1))*x(y1?) = -num(?n,ncmd(?n1,y1));
        endif;


        id once num(0,?n,ncmd(?n1))*x(y1?)*xbar(y2?) = 
                     -num(0,?n,ncmd(?n1))*Eres(y1,y2)*xbar(y2);
    
        id once num(1,?n,ncmd(?n1))*x(y1?)*xbar(y2?) = 
                     +num(0,?n,ncmd(?n1,y1))*Eres(y1,y2)*xbar(y2)
                     -num(1,?n,ncmd(?n1,y1))*xbar(y2);

        repeat id Eres(y1?,y2?)*xbar(y2?)*xbar(y3?)= 
                      den(y1-y2)*xbar(y3)*(xbar(y2)+Eres(y1,y3));
        id once Eres(y1?,y2?)*xbar(y2?)= den(y1-y2)*xbar(y2);
    endrepeat;
    id num(n0?{0,1},?n) = num(?n);
    .sort:recursive-expansion;

* Sanity check
    if (count(xbar,1) || count(x,1));
         print "[k`i'] Some roots were not absorbed: %t";
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

* Expand numerator
*#define numEXPANSION "0";
#ifdef `numEXPANSION'
	#message "Start expanding the numerator"
	#redefine oldextrasymbols "`extrasymbols_'"
	ExtraSymbols, underscore, den;
	argtoextrasymbol den;
	id den(y?) = y;
	
	.sort:den;
	off statistics;
	#include pf_numerator.frm # pf_num 
	on statistics; 
	
	multiply replace_(<den{`oldextrasymbols'+1}_,den(extrasymbol_({`oldextrasymbols'+1}))>\
	                  ,...,<den`extrasymbols_'_,den(extrasymbol_(`extrasymbols_'))>);
	
    Format mathematica;
	.sort:num-expansion;
#endif


#write<form_pf_`topoID'.out> "%e", F;

*
** Compare
*#message "Compare against reference expression"
*factarg (-1) norm 1;
*id norm(y?,n?) = sign_(n);
*.sort
*#include reference.frm
*L diff = F - refF;
*B+ num;
**print +s diff;
*.sort
*
*#write<diff.out> "L F = %e", diff;
*
.end
