#procedure Gstring(G,par)
*
*	When par == 0, we just manipulate a string of gamma matrices
*	               G should be just a tensor.
*	When par == 1, we take the D-dimensional trace.
*	               G should be a cyclesymmetric tensor.
*	               We also need an auxiliary function GGstring.
*	We assume dimension D = 4-2*ep
*	Indices are mu1,...,mun, but they are assumed to have been summed over.
*	p should be a vector
*	There is a polyratfun rat.
*	Other assumption: the indices mu occur only in the Gs and are contracted.
*	We also need a function NN, a symbol n and symbols x1,x2.
*
*	Strategy: we mimic the tracen algorithms (par == 1) or just
*	the simplifying parts of them when par == 0.
*
#switch `par'
#case 1
*--#[ Trace :
*
*	Here we do the trace.
*	We make sure that in each module one pair of identical indices
*	or identical vectors is eliminated. Eventually the string has
*	only different objects and we have to use the tracen statement.
*
#$pass = 0;
#do iGstring = 1,1
#$pass = $pass+1;
#$catch = 0;
#$num = 0;
*
label easycatch;
*
*	Distance zero and a match that generates only a single term
*
id,ifmatch->easycatch,`G'(p?,p?,?b)=`G'(?b)*p.p;
id,ifmatch->easycatch,`G'(mu?,mu?,?b)=`G'(?b)*rat(4-2*ep,1);
id,ifmatch->easycatch,`G'(mu?,mu1?,mu?,?b)=rat(-2+2*ep,1)*`G'(mu1,?b);
*
*	Now the matches that generate two terms
*
id ifmatch->catch `G'(p?,mu1?,p?,?b)= 2*`G'(p,?b)*p(mu1)-p.p*`G'(mu1,?b);
id ifmatch->catch `G'(mu?,mu1?,mu2?,mu?,?b)=
			+`G'(?b)*4*d_(mu1,mu2)
			+rat(-2*ep,1)*`G'(mu1,mu2,?b);
id ifmatch->catch `G'(mu?,mu1?,mu2?,mu3?,mu?,?b)=
			-2*`G'(mu3,mu2,mu1,?b)
			+rat(2*ep,1)*`G'(mu1,mu2,mu3,?b);
*
*	Three terms
*
id ifmatch->catch `G'(p?,mu1?,mu2?,p?,?b)=
			+2*p(mu2)*`G'(p,mu1,?b)
			-2*p(mu1)*`G'(p,mu2,?b)
			+p.p*`G'(mu1,mu2,?b);
*
*	Intermezzo, to determine the number of arguments in `G'
id	`G'(?a) = `G'(?a)*NN(nargs_(?a));
id	NN(n?$num) = 1;
if ( $num >= 10 )
id ifmatch->catch `G'(mu?,mu1?,mu2?,mu3?,mu4?,mu?,?b) =
		+2*`G'(mu4,mu1,mu2,mu3,?b)
        +2*`G'(mu3,mu2,mu1,mu4,?b)
		+rat(-2*ep,1)*`G'(mu1,mu2,mu3,mu4,?b);
*
*	Four terms
*
if ( $num >= 8 )
id ifmatch->catch `G'(p?,mu1?,mu2?,mu3?,p?,?b)=
		+2*p(mu3)*`G'(p,mu1,mu2,?b)
		-2*p(mu2)*`G'(p,mu1,mu3,?b)
		+2*p(mu1)*`G'(p,mu2,mu3,?b)
		-p.p*`G'(mu1,mu2,mu3,?b);
if ( $num >= 12 )
id ifmatch->catch `G'(mu?,mu1?,...,mu5?,mu?,?b) =
		+2*`G'(mu5,mu1,...,mu4,?b)
		-2*`G'(mu4,mu1,...,mu3,mu5,?b)
		-2*`G'(mu3,mu2,mu1,mu4,mu5,?b)
		+rat(2*ep,1)*`G'(mu1,...,mu5,?b);
*
*	Five terms
*
if ( $num >= 10 )
id ifmatch->catch `G'(p?,mu1?,...,mu4?,p?,?b)=
		+2*p(mu4)*`G'(p,mu1,...,mu3,?b)
		-2*p(mu3)*`G'(p,mu1,mu2,mu4,?b)
		+2*p(mu2)*`G'(p,mu1,mu3,mu4,?b)
		-2*p(mu1)*`G'(p,mu2,mu3,mu4,?b)
		+p.p*`G'(mu1,mu2,mu3,mu4,?b);
if ( $num >= 14 )
id ifmatch->catch `G'(mu?,mu1?,...,mu6?,mu?,?b) =
		+2*`G'(mu6,mu1,...,mu5,?b)
		-2*`G'(mu5,mu1,...,mu4,mu6,?b)
		+2*`G'(mu4,mu1,...,mu3,mu5,mu6,?b)
		+2*`G'(mu3,mu2,mu1,mu4,mu5,mu6,?b)
		+rat(-2*ep,1)*`G'(mu1,...,mu6,?b);
*
*	Six terms
*
if ( $num >= 12 )
id ifmatch->catch `G'(p?,mu1?,...,mu5?,p?,?b)=
		+2*p(mu5)*`G'(p,mu1,...,mu4,?b)
		-2*p(mu4)*`G'(p,mu1,...,mu3,mu5,?b)
		+2*p(mu3)*`G'(p,mu1,mu2,mu4,mu5,?b)
		-2*p(mu2)*`G'(p,mu1,mu3,mu4,mu5,?b)
		+2*p(mu1)*`G'(p,mu2,mu3,mu4,mu5,?b)
		-p.p*`G'(mu1,mu2,mu3,mu4,mu5,?b);
if ( $num >= 16 )
id ifmatch->catch `G'(mu?,mu1?,...,mu7?,mu?,?b) =
		+2*`G'(mu7,mu1,...,mu6,?b)
		-2*`G'(mu6,mu1,...,mu5,mu7,?b)
		+2*`G'(mu5,mu1,...,mu4,mu6,mu7,?b)
		-2*`G'(mu4,mu1,...,mu3,mu5,mu6,mu7,?b)
		-2*`G'(mu3,mu2,mu1,mu4,mu5,mu6,mu7,?b)
		+rat(2*ep,1)*`G'(mu1,...,mu7,?b);
*
*	Seven terms
*
if ( $num >= 14 )
    id ifmatch->catch `G'(p?,mu1?,...,mu6?,p?,?b)=
		+2*p(mu6)*`G'(p,mu1,...,mu5,?b)
		-2*p(mu5)*`G'(p,mu1,...,mu4,mu6,?b)
		+2*p(mu4)*`G'(p,mu1,...,mu3,mu5,mu6,?b)
		-2*p(mu3)*`G'(p,mu1,mu2,mu4,mu5,mu6,?b)
		+2*p(mu2)*`G'(p,mu1,mu3,mu4,mu5,mu6,?b)
		-2*p(mu1)*`G'(p,mu2,mu3,mu4,mu5,mu6,?b)
		+p.p*`G'(mu1,mu2,mu3,mu4,mu5,mu6,?b);
*
*	Eight terms
*
if ( $num >= 16 )
    id ifmatch->catch `G'(p?,mu1?,...,mu7?,p?,?b)=
		+2*p(mu7)*`G'(p,mu1,...,mu6,?b)
		-2*p(mu6)*`G'(p,mu1,...,mu5,mu7,?b)
		+2*p(mu5)*`G'(p,mu1,...,mu4,mu6,mu7,?b)
		-2*p(mu4)*`G'(p,mu1,...,mu3,mu5,mu6,mu7,?b)
		+2*p(mu3)*`G'(p,mu1,mu2,mu4,mu5,mu6,mu7,?b)
		-2*p(mu2)*`G'(p,mu1,mu3,mu4,mu5,mu6,mu7,?b)
		+2*p(mu1)*`G'(p,mu2,mu3,mu4,mu5,mu6,mu7,?b)
		-p.p*`G'(mu1,mu2,mu3,mu4,mu5,mu6,mu7,?b);
*
*	If the above did not catch anything, but there are still repeated
*	indices or vectors things are a bit more complicated than with the
*	strings, because in principle we could get an infinite loop, due to
*	the cyclic property of G.
*	We need the shortest distance!
*
if ( $num > 16 );
  if ( match(`G'(mu?,?a,mu?,?b)) );
	id,once,`G'(mu?,?a,mu?,?b) =
			GGstring(nargs_(?a),mu,?a)*GGstring(nargs_(?b),mu,?b);
	id,ifmatch->catch,GGstring(x1?,mu?,?a,mu1?)*GGstring(x2?,mu?,?b) =
			2*`G'(mu1,?a,?b)-`G'(mu,?a,mu,mu1,?b);
  elseif ( match(`G'(p?,?a,p?,?b)) );
	id,once,`G'(p?,?a,p?,?b) =
			GGstring(nargs_(?a),p,?a)*GGstring(nargs_(?b),p,?b);
	id,ifmatch->catch,GGstring(x1?,p?,?a,mu1?)*GGstring(x2?,p?,?b) =
			2*p(mu1)*`G'(p,?a,?b)-`G'(p,?a,p,mu1,?b);
  endif;
endif;
*
goto nocatch;
*
label catch;
	$catch = 1;
	rcyclesymmetrize `G';
label nocatch;
B+	`G';
moduleoption maximum,$catch;
moduleoption local,$num;
.sort: trace pass `$pass';
#if ( `$catch' > 0 )
	#redefine iGstring "0"
#endif
Keep Brackets;
#enddo
*
*	Now cleanup and after that the remaining trace
*
rcyclesymmetrize `G';
renumber;
B+ `G';
.sort:simplefilter;
Keep brackets;
id	`G' = 4;
repeat;
  id,once,`G'(?a) = g_(1,?a);
  tracen,1;
  id  D = rat(4-2*ep,1);
endrepeat;
.sort:trace;
*
*	Done
*
*--#] Trace : 
#break
#case 0
*--#[ String :
*
*	Here we do just simplify strings.
*	This case is more complicated than the trace, because we do not have
*	cyclic symmetry and hence the identical objects can be further away
*	from each other. Hence we have a generic code at the end. It may not
*	be super fast, but at least it will work.
*
#$pass = 0;
#do iGstring = 1,1
#$pass = $pass+1;
#$catch = 0;
#$num = 0;
*
label easycatch;
*
*	Distance zero and a match that generates only a single term
*
id,ifmatch->easycatch,`G'(?a,p?,p?,?b)=`G'(?a,?b)*p.p;
id,ifmatch->easycatch,`G'(?a,mu?,mu?,?b)=`G'(?a,?b)*rat(4-2*ep,1);
id,ifmatch->easycatch,`G'(?a,mu?,mu1?,mu?,?b)=rat(-2+2*ep,1)*`G'(?a,mu1,?b);
*
*	Now the matches that generate two terms
*
id ifmatch->catch `G'(?a,p?,mu1?,p?,?b)= 2*`G'(?a,p,?b)*p(mu1)-p.p*`G'(?a,mu1,?b);
id ifmatch->catch `G'(?a,mu?,mu1?,mu2?,mu?,?b)=
			+`G'(?a,?b)*4*d_(mu1,mu2)
			+rat(-2*ep,1)*`G'(?a,mu1,mu2,?b);
id ifmatch->catch `G'(?a,mu?,mu1?,mu2?,mu3?,mu?,?b)=
			-2*`G'(?a,mu3,mu2,mu1,?b)
			+rat(2*ep,1)*`G'(?a,mu1,mu2,mu3,?b);
*
*	Three terms
*
id ifmatch->catch `G'(?a,p?,mu1?,mu2?,p?,?b)=
			+2*p(mu2)*`G'(?a,p,mu1,?b)
			-2*p(mu1)*`G'(?a,p,mu2,?b)
			+p.p*`G'(?a,mu1,mu2,?b);
*
*	Intermezzo, to determine the number of arguments in `G'
id	`G'(?a) = `G'(?a)*NN(nargs_(?a));
id	NN(n?$num) = 1;
if ( $num >= 6 )
id ifmatch->catch `G'(?a,mu?,mu1?,mu2?,mu3?,mu4?,mu?,?b) =
		+2*`G'(?a,mu4,mu1,mu2,mu3,?b)
        +2*`G'(?a,mu3,mu2,mu1,mu4,?b)
		+rat(-2*ep,1)*`G'(?a,mu1,mu2,mu3,mu4,?b);
*
*	Four terms
*
if ( $num >= 5 )
id ifmatch->catch `G'(?a,p?,mu1?,mu2?,mu3?,p?,?b)=
		+2*p(mu3)*`G'(?a,p,mu1,mu2,?b)
		-2*p(mu2)*`G'(?a,p,mu1,mu3,?b)
		+2*p(mu1)*`G'(?a,p,mu2,mu3,?b)
		-p.p*`G'(?a,mu1,mu2,mu3,?b);
if ( $num >= 7 )
id ifmatch->catch `G'(?a,mu?,mu1?,...,mu5?,mu?,?b) =
		+2*`G'(?a,mu5,mu1,...,mu4,?b)
		-2*`G'(?a,mu4,mu1,...,mu3,mu5,?b)
		-2*`G'(?a,mu3,mu2,mu1,mu4,mu5,?b)
		+rat(2*ep,1)*`G'(?a,mu1,...,mu5,?b);
*
*	Five terms
*
if ( $num >= 6 )
id ifmatch->catch `G'(?a,p?,mu1?,...,mu4?,p?,?b)=
		+2*p(mu4)*`G'(?a,p,mu1,...,mu3,?b)
		-2*p(mu3)*`G'(?a,p,mu1,mu2,mu4,?b)
		+2*p(mu2)*`G'(?a,p,mu1,mu3,mu4,?b)
		-2*p(mu1)*`G'(?a,p,mu2,mu3,mu4,?b)
		+p.p*`G'(?a,mu1,mu2,mu3,mu4,?b);
if ( $num >= 8 )
id ifmatch->catch `G'(?a,mu?,mu1?,...,mu6?,mu?,?b) =
		+2*`G'(?a,mu6,mu1,...,mu5,?b)
		-2*`G'(?a,mu5,mu1,...,mu4,mu6,?b)
		+2*`G'(?a,mu4,mu1,...,mu3,mu5,mu6,?b)
		+2*`G'(?a,mu3,mu2,mu1,mu4,mu5,mu6,?b)
		+rat(-2*ep,1)*`G'(?a,mu1,...,mu6,?b);
*
*	Six terms
*
if ( $num >= 7 )
id ifmatch->catch `G'(?a,p?,mu1?,...,mu5?,p?,?b)=
		+2*p(mu5)*`G'(?a,p,mu1,...,mu4,?b)
		-2*p(mu4)*`G'(?a,p,mu1,...,mu3,mu5,?b)
		+2*p(mu3)*`G'(?a,p,mu1,mu2,mu4,mu5,?b)
		-2*p(mu2)*`G'(?a,p,mu1,mu3,mu4,mu5,?b)
		+2*p(mu1)*`G'(?a,p,mu2,mu3,mu4,mu5,?b)
		-p.p*`G'(?a,mu1,mu2,mu3,mu4,mu5,?b);
if ( $num >= 9 )
id ifmatch->catch `G'(?a,mu?,mu1?,...,mu7?,mu?,?b) =
		+2*`G'(?a,mu7,mu1,...,mu6,?b)
		-2*`G'(?a,mu6,mu1,...,mu5,mu7,?b)
		+2*`G'(?a,mu5,mu1,...,mu4,mu6,mu7,?b)
		-2*`G'(?a,mu4,mu1,...,mu3,mu5,mu6,mu7,?b)
		-2*`G'(?a,mu3,mu2,mu1,mu4,mu5,mu6,mu7,?b)
		+rat(2*ep,1)*`G'(?a,mu1,...,mu7,?b);
*
*	Seven terms
*
if ( $num >= 8 )
    id ifmatch->catch `G'(?a,p?,mu1?,...,mu6?,p?,?b) =
		+2*p(mu6)*`G'(?a,p,mu1,...,mu5,?b)
		-2*p(mu5)*`G'(?a,p,mu1,...,mu4,mu6,?b)
		+2*p(mu4)*`G'(?a,p,mu1,...,mu3,mu5,mu6,?b)
		-2*p(mu3)*`G'(?a,p,mu1,mu2,mu4,mu5,mu6,?b)
		+2*p(mu2)*`G'(?a,p,mu1,mu3,mu4,mu5,mu6,?b)
		-2*p(mu1)*`G'(?a,p,mu2,mu3,mu4,mu5,mu6,?b)
		+p.p*`G'(?a,mu1,mu2,mu3,mu4,mu5,mu6,?b);
*
*	Eight terms
*
if ( $num >= 9 )
    id ifmatch->catch `G'(?a,p?,mu1?,...,mu7?,p?,?b) =
		+2*p(mu7)*`G'(?a,p,mu1,...,mu6,?b)
		-2*p(mu6)*`G'(?a,p,mu1,...,mu5,mu7,?b)
		+2*p(mu5)*`G'(?a,p,mu1,...,mu4,mu6,mu7,?b)
		-2*p(mu4)*`G'(?a,p,mu1,...,mu3,mu5,mu6,mu7,?b)
		+2*p(mu3)*`G'(?a,p,mu1,mu2,mu4,mu5,mu6,mu7,?b)
		-2*p(mu2)*`G'(?a,p,mu1,mu3,mu4,mu5,mu6,mu7,?b)
		+2*p(mu1)*`G'(?a,p,mu2,mu3,mu4,mu5,mu6,mu7,?b)
		-p.p*`G'(?a,mu1,mu2,mu3,mu4,mu5,mu6,mu7,?b);
*
*	Generic reduction
*
id ifmatch->catch `G'(?a,mu?,?b,mu1?,mu?,?c) =
		+2*`G'(?a,mu1,?b,?c) - `G'(?a,mu,?b,mu,mu1,?c);
id ifmatch->catch `G'(?a,p?,?b,mu?,p?,?c) =
		+2*p(mu)*`G'(?a,p,?b,?c) - `G'(?a,p,?b,p,mu,?c);
goto nocatch;
*
label catch;
	$catch = 1;
label nocatch;
B+	`G';
moduleoption maximum,$catch;
moduleoption local,$num;
.sort: string pass `$pass';
#if ( `$catch' > 0 )
	#redefine iGstring "0"
#endif
Keep Brackets;
#enddo
*
*	Now cleanup and after that the remaining trace
*
.sort:strings;
*
*	Done
*
*--#] String : 
#break
#default
#message The case par = `par' has not been implemented in Gstring.
#message Program aborting.....
#terminate -1
#break
#endswitch
#endprocedure
