#-
* PDGs:
*d,u,c,s,b,t -> 1 to 6
*d~,u~,c~,s~,b~,t~ -> -1 to -6
*g, photon -> 21, 22
*e+ e- > -11, 11
*mu+, mu-, ta+ ,ta- > -12, 12, -13, 13
*z w+ w- z -> 23, 24, -24
*H 25

#define GLU "21"
#define PHO "22"
#define EP "-11"
#define EM "11"
#define H "25"
#define GHO "82"
#define GHOBAR "-82"
#define FERM "-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,-11,11,-12,12,-13,13"
#define Q "1,2,3,4,5,6"
#define QBAR "-1,-2,-3,-4,-5,-6"

Auto S mass;
CTable masses(-30:30);

#ifndef `OPTIMITERATIONS'
    #define OPTIMITERATIONS "100"
#endif

#ifndef `HEAVYFERMIONS'
Fill masses(1) = 0;
Fill masses(2) = 0;
Fill masses(3) = 0;
Fill masses(4) = 0;
Fill masses(5) = 0;
Fill masses(-1) = 0;
Fill masses(-2) = 0;
Fill masses(-3) = 0;
Fill masses(-4) = 0;
Fill masses(-5) = 0;
Fill masses(11) = 0;
Fill masses(12) = 0;
Fill masses(13) = 0;
Fill masses(-11) = 0;
Fill masses(-12) = 0;
Fill masses(-13) = 0;
#else
Fill masses(1) = mass_d;
Fill masses(2) = mass_u;
Fill masses(3) = mass_c;
Fill masses(4) = mass_s;
Fill masses(5) = mass_b;
Fill masses(-1) = mass_d;
Fill masses(-2) = mass_u;
Fill masses(-3) = mass_c;
Fill masses(-4) = mass_s;
Fill masses(-5) = mass_b;
Fill masses(11) = mass_e;
Fill masses(12) = mass_mu;
Fill masses(13) = mass_tau;
Fill masses(-11) = mass_e;
Fill masses(-12) = mass_mu;
Fill masses(-13) = mass_tau;
#endif

Fill masses(6) = mass_t;
Fill masses(-6) = mass_t;
Fill masses(25) = mass_h;

Auto V p,k;
Auto S lm,ext;
Auto I mu=4,s=4;
Symbol ge, gs, gy, type, in, out, virtual;
Auto S x, idx;

Set dirac: s1,...,s40;
Set lorentz: mu1,...,mu40;
Set lorentzdummy: mud1,...,mud40;

CF gamma, vector,g(s),delta(s),T, counter,color, prop;
CF u,ubar,v,vbar,f,vx, vec;
Symbol ca,cf,nf,[dabc^2/n],[d4RR/n],[d4RA/n],[d4AA/n];

S  i, m, n;

#include- diacolor.h
Set colF: cOli1,...,cOli40;
Set colA: cOlj1,...,cOlj40;
Set colAdum: cOljj1,...,cOljj40;

* Load the diaggrams
#include- input.h

************************************************
* Substitute the Feynman rules for the numerator
************************************************

* Fix a quirk where 0 does not match to a vector
* The only 0 in a propagator or vertex is when a momentum is 0
* All indices and pdgs are non-zero
repeat id prop(?a, 0, ?b) = prop(?a, pzero, ?b);
repeat id vx(?a, 0, ?b) = vx(?a, pzero, ?b);

* incoming edges
id prop(`EM', in, p?, idx?) = u(p, masses(`EM'), dirac[idx]);
id prop(x?{`Q'}, in, p?, idx?) = u(p, masses(x), dirac[idx]);
id prop(`EP', in, p?, idx?) = vbar(p, masses(`EP'), dirac[idx]);
id prop(x?{`QBAR'}, in, p?, idx?) = vbar(p, masses(x), dirac[idx]);

* outgoing edges
id prop(`EM', out, p?, idx?) = ubar(p, masses(`EM'), dirac[idx]);
id prop(`EP', out, p?, idx?) = v(p, masses(`EP'), dirac[idx]);
id prop(x?{`Q'}, out, p?, idx?) = ubar(p, masses(x), dirac[idx]);
id prop(x?{`Q'}, out, p?, idx?) = v(p, masses(x), dirac[idx]);

* virtual edges
id prop(`GLU', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]) * d_(colA[idx1], colA[idx2]);
id prop(x?{`GHO',`GHOBAR'}, virtual, p?, idx1?, idx2?) = i_ * d_(colA[idx1], colA[idx2]);
id prop(`PHO', virtual, p?, idx1?, idx2?) = - i_ * d_(lorentz[idx1], lorentz[idx2]);
id prop(x?{`FERM'}, virtual, p?, idx1?, idx2?) = - i_ * (gamma(dirac[idx1], p, dirac[idx2]) + masses(x) * gamma(dirac[idx1], dirac[idx2])) * d_(colF[idx1], colF[idx2]);
id prop(`H', virtual, p?, idx1?, idx2?) = -i_;

.sort:feynman-rules-edges;

* vertices
id vx(x1?{`QBAR'}, `GLU', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = gs * i_ * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * T(colF[idx1], colA[idx2], colF[idx3]);
id vx(`GHOBAR', `GLU', `GHO', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gs * cOlf(colA[idx3], colA[idx1], colA[idx2]) * p3(lorentz[idx2]);
id vx(x1?{`QBAR'}, `PHO', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = 2/3 * ge * i_* gamma(dirac[idx1], lorentz[idx2], dirac[idx3]) * d_(colF[idx1], colF[idx3]);
id vx(`EP', `PHO', `EM', p1?, p2?, p3?, idx1?, idx2?, idx3?) = -ge * i_ * gamma(dirac[idx1], lorentz[idx2], dirac[idx3]);
id vx(x1?{`QBAR'}, `H', x2?{`Q'}, p1?, p2?, p3?, idx1?, idx2?, idx3?) = -gy * i_ * d_(dirac[idx1], dirac[idx3]) * d_(colF[idx1], colF[idx3]);

* TODO: use momentum conservation to reduce the number of different terms
id vx(`GLU', `GLU', `GLU', p1?, p2?, p3?, idx1?, idx2?, idx3?) = gs * cOlf(colA[idx1], colA[idx2], colA[idx3]) *(
    - d_(lorentz[idx1], lorentz[idx3]) * p1(lorentz[idx2])
    + d_(lorentz[idx1], lorentz[idx2]) * p1(lorentz[idx3])
    + d_(lorentz[idx2], lorentz[idx3]) * p2(lorentz[idx1])
    - d_(lorentz[idx1], lorentz[idx2]) * p2(lorentz[idx3])
    - d_(lorentz[idx2], lorentz[idx3]) * p3(lorentz[idx1])
    + d_(lorentz[idx1], lorentz[idx3]) * p3(lorentz[idx2]));

* For the quartic gluon vertex we need an extra dummy index
Multiply counter(1);
repeat id vx(`GLU', `GLU', `GLU', `GLU', p1?, p2?, p3?, p4?, idx1?, idx2?, idx3?, idx4?)*counter(i?) = counter(i + 1) * gs^2 * i_ *(
    + cOlf(colA[idx1], colA[idx2], colAdum[i]) * cOlf(colA[idx3], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]))
    + cOlf(colA[idx1], colA[idx3], colAdum[i]) * cOlf(colA[idx2], colA[idx4], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx4]) * d_(lorentz[idx2], lorentz[idx3]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4]))
    + cOlf(colA[idx1], colA[idx4], colAdum[i]) * cOlf(colA[idx2], colA[idx3], colAdum[i])
        * (d_(lorentz[idx1], lorentz[idx3]) * d_(lorentz[idx2], lorentz[idx4]) - d_(lorentz[idx1], lorentz[idx2]) * d_(lorentz[idx3], lorentz[idx4])));
id counter(x?) = 1;

if (count(vx, 1, prop, 1));
    Print "Unsubstituted propagator or vertex: %t";
    exit "Critical error";
endif;

.sort:feynman-rules-vertices;

******************
* Color evaluation
******************
repeat id T(cOli1?,?a,cOli2?)*T(cOli2?,?b,cOli3?) = T(cOli1,?a,?b,cOli3); * collect the colour string
id  T(cOli1?, ?a, cOli1?) = cOlTr(?a);
id  cOlTr(cOli1?) = 0;
id  cOlTr = cOlNR;
Multiply color(1);
repeat id cOlTr(?a)*color(x?) = color(x * cOlTr(?a));
repeat id cOlf(cOlj1?,cOlj2?,cOlj3?)*color(x?) = color(x * cOlf(cOlj1,cOlj2,cOlj3));
B+ color;
.sort:color-prep;
Keep brackets;
* TODO: do we need to set the dimension?
* Only evaluate this part per unique color by bracketing
Argument color;
    #call color
    #call simpli
    id  cOlI2R = cOlcR*cOlNR/cOlNA;
    id  cOlNR/cOlNA*cOlcR = cOlI2R;
    id  cOld33(cOlpR1,cOlpR2) = [dabc^2/n]*cOlNR;
    id  cOlNR/cOlNA = nf/cf/2;
    id  cOlcR = cf;
    id  cOlcA = ca;
    id  cOlI2R = nf/2;
	id	cOld44(cOlpA1,cOlpA2) = [d4AA/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpR2) = [d4RR/n]*cOlNA;
	id	cOld44(cOlpR1,cOlpA1) = [d4RA/n]*cOlNA;

* set the SU(3) values
    id cOlNR = 3;
    id cOlNA = 8;
    id cf = 4 / 3;
    id ca = 3;
EndArgument;
.sort:color;

* set the SU(3) values
id cOlNR = 3;

* Do the spin sum
repeat id u(p?,m?,s1?)*ubar(p?,m?,s2?) = gamma(s1, p, s2) + m*gamma(s1, s2);
repeat id v(p?,m?,s1?)*vbar(p?,m?,s2?) = gamma(s1, p, s2) - m*gamma(s1, s2);

* construct gamma string
repeat id gamma(s1?,?a,s2?)*gamma(s2?,?b,s3?) = gamma(s1,?a,?b,s3);

Multiply counter(1);
repeat id gamma(mu?,?a,p?,?b,mu?)*counter(i?) = vec(p,lorentzdummy[i])*gamma(mu,?a,lorentzdummy[i],?b,mu)*counter(i+1);
id counter(n?) = 1;

#do i=1,10
    id once gamma(mu?,?a,mu?) = g_(`i',?a);
#enddo

id vec(p?,mu?) = p(mu);

#do i=1,10
    trace4 `i';
#enddo

if (count(gamma, 1));
    Print "Unsubstituted gamma string: %t";
    exit "Critical error";
endif;

.sort:gamma-traces;

id color(x?) = x;

* Set all external momenta on-shell
#do i=1,10
    id p`i'.p`i' = 0;
#enddo

id pzero = 0; * Substitute the 0-momentum by 0
.sort:feynman-rules-final;

*********************************************
* Construction of optimized numerator C code
*********************************************

* TODO: split off the energy components of the LTD loop momenta
* and optimize a monomial in the energies

* Translate all dot products to a linear representation
* We could also write out the dot products, at the cost of a drastic increase in variables, but does it help?
* It could do k1.p1+k1.p2 => k1.(p1+p2)
#$MAXK = 0;
#do i=10,0,-1
    if (count(k`i', 1));
        if ($MAXK < `i');
            $MAXK = `i';
        endif;
        goto donek;
    endif;
#enddo
label donek;

#$MAXP = 0;
#do i=10,0,-1
    if (count(p`i', 1));
        if ($MAXP < `i');
            $MAXP = `i';
        endif;
        goto donep;
    endif;
#enddo
label donep;
.sort

#$OFFSET = 0;
#do i=1,`$MAXP'
    #do j=`i',`$MAXP'
        id p`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

#do i=1,`$MAXK'
    #do j=1,`$MAXP'
        id k`i'.p`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo

    #do j=`i',`$MAXK'
        id k`i'.k`j' = lm`$OFFSET';
        #$OFFSET = $OFFSET + 1;
    #enddo
#enddo

Format float 16; * print all constants as floats in the C output
Format C;
.sort
Format O4,stats=on,saIter=`OPTIMITERATIONS';

#Optimize F
#write<out.proto_c> "%O\n\treturn %e",F
.end