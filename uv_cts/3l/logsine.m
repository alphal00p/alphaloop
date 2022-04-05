(* ::Package:: *)

BeginPackage["LsToLi`"];


Print["LsToLi: evaluating log-sine integrals in polylogarithmic terms
	accompanying the paper \"Special values of generalized log-sine integrals\"
	-- Jonathan M. Borwein, University of Newcastle
	-- Armin Straub, Tulane University
	-- Version 2.0 (2013/04/03)"];


Ls::usage=
	"Ls[n,k,s] represents the generalized log-sine integral.
	Usage:
		N[Ls[5,1,Pi/3]]";
Lsh::usage=
	"Lsh[n,k,s] represents the generalized log-sinh integral.";


Li::usage=
	"Li[{a1,..,an},x] represents the multiple polylogarithm.";
Cl::usage=
	"Cl[{a1,..,an},t] represents the multiple Clausen function.";
Gl::usage=
	"Gl[{a1,..,an},t] represents the multiple Glaisher function.";


LsExpand::usage=
	"Normalize the arguments of generalized log-sine integrals to lie in the interval [0,Pi].";
LsToLi::usage=
	"Evaluate generalized log-sine integral in polylogarithmic terms.";
LiReduce::usage=
	"Reduce polylogarithmic terms, if possible, to lower dimensional ones.";
NielsenS1::usage=
	"The Nielsen polylogarithm at argument 1, reduced to zeta values.";


Options[LsToLi]={
	UseReductionTable->True};


Options[LiReduce]={
	UseReductionTable->True};


Begin["Private`"];


(* ::Section:: *)
(*Top level*)


LsToLi[expr_,OptionsPattern[]]:=Module[{x=expr},
	(* Reduce arguments to range [0,Pi]. *)
	x=LsExpand[x];
	(* The next step is not used anymore since the resulting polylogs at -1 need
		a different approach for reduction. Try e.g. Ls[5,1,Pi]. *)
	(* Reduce values at Pi. *)
	(* x=Expand[x/.Ls[n_,k_,Pi]:>LsPiToLi[n,k]]; *)
	(* Reduce values at other arguments. *)
	x=Expand[x/.Ls[n_,k_,t_/;Simplify[0<=t<=2Pi]]:>LsGeneralToLi[n,k,t]];
	x=LiReduce[x];
	If[OptionValue[UseReductionTable],
		x=x/.LiReductionTable
	];
	Expand[x]
];


LsExpand[expr_]:=Module[{x=expr},
	(* Treat negative arguments. *)
	x=x/.Ls[n_,k_,s_/;s<0]:>(-1)^(k+1)Ls[n,k,-s];
	(* Reduce arguments to range [0,Pi] and multiples of 2Pi. *)
	x=x/.Ls[n_,k_,s_/;s>Pi&&!Element[s/(2Pi),Integers]]:>LsReduceArgument[n,k,s];
	(* Reduce values at multiples of 2Pi. *)
	x=x/.Ls[n_,k_,Pi a_?EvenQ/;a>2]:>LsMult2PiTo2Pi[n,k,a/2];
	x=x/.Ls[n_,k_,2Pi]:>Ls2Pi[n,k];
	x=x/.Ls[n_,0,Pi]:>LsPi0[n];
	(* Do trivial evaluations. *)
	x=Expand[x/.{Ls[n0_,k0_,t_]/;(n0==k0+1):>-1/n0(t)^n0}];
	x
];


(* ::Section:: *)
(*Basic numerics*)


Ls/:N[Ls[n_,k_,s_],prec_?NumericQ]:=-NIntegrate[t^k Log[Abs[2Sin[t/2]]]^(n-k-1),{t,0,s},WorkingPrecision->prec];
Lsh/:N[Lsh[n_,k_,s_],prec__?NumericQ]:=-NIntegrate[t^k Log[Abs[2Sinh[t/2]]]^(n-k-1),{t,0,s},WorkingPrecision->prec];
Ls/:N[Ls[n_,k_,s_]]:=N[Ls[n,k,s],$MachinePrecision];
Lsh/:N[Lsh[n_,k_,s_]]:=N[Lsh[n,k,s],$MachinePrecision];


NielsenQ[L_]:=Length[L]>=1&&Complement[Rest[L],{1}]==={};
N[Li[L_List,x_],prec_:$MachinePrecision]:=N[PolyLog[L[[1]]-1,Length[L],x],prec]/;NielsenQ[L];
N[Cl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Im[p],Re[p]]]/;NielsenQ[L];
N[Gl[L_List,x_],prec_:$MachinePrecision]:=Module[{p=N[PolyLog[L[[1]]-1,Length[L],Exp[I x]],prec]},If[EvenQ[Total[L]],Re[p],Im[p]]]/;NielsenQ[L];


(* ::Section:: *)
(*At pi and 2pi*)


LsPi0[n_/;n>=3]:=LsPi0[n]=Module[{a,k},
	a[m_]:=(1-2^(1-m))Zeta[m];
	(-1)^(n)(n-2)!(Pi a[n-1]+Sum[(-1)^k/(k+1)!a[n-2-k]LsPi0[k+2],{k,1,n-4}])//Expand];
LsPi0[2]:=0;
LsPi0[1]:=-Pi;


Ls2Pi[n_,k_]:=Ls2Pi[n,k]=Module[{m=n-k-1,a,expa0,idxcoeff,expa,x,ls},
	a[nn_]:=-2Sum[(-1/2)^r/r!(nn+r-1)!Zeta[nn+r]x^r,{r,0,m}]; (*TYPO in Lewin!*)
	a[1]=I Pi;
	(* The below can surely be optimized by extracting only the needed coefficient. *)
	expa0=1-1/Pi Sum[LsPi0[r+1]/r!x^r,{r,2,m}];
	idxcoeff[ii_]:=1/Times@@((#[[2]]!)&/@Tally[ii]);
	expa=Sum[Sum[idxcoeff[ii]Product[a[i]/i!,{i,ii}],{ii,IntegerPartitions[k,{j},Join[{1},Table[r,{r,2,k,2}]]]}],{j,0,k}];
	ls=Coefficient[expa0 expa,x,m];
	ls=Expand[-2Pi (-I)^k k!m!ls];
	ls
];


LsPiToLi[n_,k_]:=Module[{m=n-k-1,ls,L},
	If[k==0,Return[LsPi0[n]]];
	(* Contribution from the first term *)
	ls=If[Mod[m,2]==0,Pi^(m+k+1)/2^m Beta[m+1,k+1](-1)^Floor[(m+k)/2],0];
	(* Polylogarithmic contribution from the other terms *)
	L[a_,b_,x_]:=Li[Join[{a},Table[1,{b}]],x];
	ls+=(-1)^m m!(
		Sum[(-1)^c/2^c Pochhammer[c+1,k-d]Binomial[k,d](-1)^Ceiling[d/2]Pi^d L[c+k-d+2,m-c-1,-1],{c,0,m-1},{d,1+Mod[k,2],k,2}]
		-Sum[(-1)^c/2^c Pochhammer[c+1,k]/b!(-1)^Ceiling[b/2](Pi/2)^b NielsenS1[c+k+1,m-c-b],{c,0,m-2},{b,1+Mod[k,2],m-c-1,2}]
		-If[Mod[k,2]==0,0,Sum[(-1)^c/2^c Pochhammer[c+1,k](L[c+k+2,m-c-1,1]-L[c+k+2,m-c-1,-1]),{c,0,m-1}]]);
	(* Adjust for the fact the egf has the second variable imaginary *)
	Expand[-(-1)^Floor[k/2]ls]
];


(* ::Section:: *)
(*Quasiperiodicity*)


LsMult2PiTo2Pi[n_,k_,m_]:=Sum[(-1)^(k-j)(2Pi)^j Binomial[k,j]Ls[n-j,k-j,2Pi]HarmonicNumber[m,-j],{j,0,k}];


LsReduceArgument[n_,k_,s_]:=Module[{m,s0,pm=1},
	(* Write s = 2m*Pi +- s0 *)
	m=Floor[s/(2Pi)];
	s0=s-2m Pi;
	If[s0>Pi,s0=2Pi-s0;pm=-1;m++];
	Ls[n,k,2m Pi]+Sum[pm^(k-j+1)(2m Pi)^j Binomial[k,j]Ls[n-j,k-j,s0],{j,0,k}]
];


(* ::Section:: *)
(*At general values*)


LsToLiStep[n_,k_,t_]:=Module[{ls},
	ls=(-1)^(n+k)(n-k-1)!k!Re[(I)^(k+1)]NielsenS1[n-k-1,k+1]
		(* Observe that on the next line Cl is only created with odd depth, and Gl only with even depth.
			This means that no general polylogarithmic reductions are necessary. *)
		-(-1)^(n+k)(n-k-1)!k!Sum[(-t)^r/r!Re[(I)^(k+r+1) Li[Join[{2+k-r},Table[1,{n-k-2}]],Exp[I t]]],{r,0,k}]
		-Sum[Binomial[n-k-1,r]Binomial[r,m](I/2)^r(-Pi)^(r-m)Ls[n-r+m,k+m,t],{r,2,n-k-1,2},{m,0,r}];
	ls
];


LsGeneralToLi[n_,k_,t_]:=Module[{ls,r},
	ls=Ls[n,k,t];
	(* Reduce the log-sines level by level. *)
	For[r=0,r<n-k-1,r+=2,
		ls=Expand[ls/.{Ls[n0_,k0_,t]/;(r==(n-k)-(n0-k0)):>LsToLiStep[n0,k0,t]}];
	];
	(* Possibly remaining log-sines on the last level are of the form Ls[m,m-1]. These evaluate trivially. *)
	ls=Expand[ls/.{Ls[n0_,k0_,t]/;(n0==k0+1):>-1/n0(t)^n0}];
	ls
];


(* ::Section:: *)
(*Reduction of polylogarithms*)


LiReduce[expr_,OptionsPattern[]]:=Module[{x=expr},
	x=FixedPoint[DoLiReduce,Expand[x]];
	If[OptionValue[UseReductionTable],
		x=x/.LiReductionTable;
		Expand[x]];		
	x
];
DoLiReduce[expr_]:=Module[{x=expr},
	x=x/.Li[L_?NielsenQ,1]:>NielsenS1[L[[1]]-1,Length[L]];
	x=x/.Li[{n_},-1]:>-(1-2^(1-n))Zeta[n];
	x=x/.Li[L_?NielsenQ,-1]/;(EvenQ[L[[1]]]):>NielsenLiMReduced[L[[1]]-1,Length[L]];
	x=x/.Gl[{n_},x_/;Simplify[Element[x,Reals]]]:>Expand[-(-1)^Floor[n/2]/n!2^(n-1)Pi^n BernoulliB[n,x/(2Pi)]];
	x=x/.Gl[L_?NielsenQ,Pi/3]/;(L[[1]]==2):>(-1)^(Length[L]+Floor[(Length[L]+2)/2])(Pi/3)^(Length[L]+1)/(2(Length[L]+1)!);
	x=x/.Cl[L_?NielsenQ,Pi/3]/;(L[[1]]==2&&OddQ[Length[L]]):>With[{n=Length[L]},-Sum[(-1)^Floor[(m-1)/2](Pi/3)^m/m!Cl[{n-m+1}, Pi/3],{m,0,n-1}]];
	x=x/.Cl[{n_?OddQ},Pi/3]:>1/2(1-2^(1-n))(1-3^(1-n))Zeta[n];
	x=x/.Cl[{n_?OddQ},2Pi/3]:>1/(1/2^(n-1)+(-1)^n)*1/2(1-2^(1-n))(1-3^(1-n))Zeta[n];
	x=x/.Cl[{n_?OddQ},Pi/2]:>2(1-2^(n-1))/4^n Zeta[n];
	x=NielsenClGlReduce[x];
	x=NielsenClGlReduceMu[x];
	Expand[x]
];


Li/:Re[Li[L_,x_/;Element[x,Reals]&&Abs[x]<=1]]:=Li[L,x];
Li/:Im[Li[L_,x_/;Element[x,Reals]&&Abs[x]<=1]]:=0;
Li/:Re[Li[L_,Exp[x_]/;Simplify[0<=x/I<=2Pi]]]:=If[EvenQ[Total[L]],Gl[L,x/I],Cl[L,x/I]];
Li/:Im[Li[L_,Exp[x_]/;Simplify[0<=x/I<=2Pi]]]:=If[EvenQ[Total[L]],Cl[L,x/I],Gl[L,x/I]];
Li/:Re[Li[L_,x_]]:=If[EvenQ[Total[L]],Gl[L,-I Log[x]],Cl[L,-I Log[x]]];
Li/:Im[Li[L_,x_]]:=If[EvenQ[Total[L]],Cl[L,-I Log[x]],Gl[L,-I Log[x]]];


(* ::Subsection:: *)
(*Reduction at 1*)


(* ::Text:: *)
(*We use K. S. Koelbig's expression (1982) for Nielsen polylogarithms at 1.*)


(* ::Text:: *)
(*The Taylor coefficients of Gamma(1+x) and 1/Gamma(1+x) respectively:*)


ga[k_]:=ga[k]=1/k Sum[(-1)^(m+1)Zeta[m]ga[k-m],{m,2,k}];
ga[0]=1;
gb[k_]:=gb[k]=-1/k Sum[(-1)^(m+1)Zeta[m]gb[k-m],{m,2,k}];
gb[0]=1;


NielsenS1[n_,p_]:=Expand[(-1)^(n+p-1)Sum[gb[n-m-1]Sum[Binomial[m+r+1,r]gb[p-r]ga[m+r+1],{r,0,p}],{m,0,n-1}]];
(*NielsenS1[n_,p_]:=Global`MZV[Join[{n+1},Table[1,{p-1}]]];*)


(* ::Subsection:: *)
(*Reduction of Cl and Gl for even/odd depth*)


NielsenClGlReduced[n_,p_,t_]:=Module[{r,s},1/2(-Sum[(t-Pi)^r/r!Binomial[n+s-r-1,s-r](-1)^(s+If[EvenQ[n],1,0]r+Floor[r/2])If[EvenQ[p],Cl,Gl][Join[{n+s-r+1},Table[1,{p-s-1}]],t],{s,1,p-1},{r,0,s}]+Sum[If[EvenQ[n+r+1],1,0](-1)^(p+Floor[r/2])(t-Pi)^r/r!CC[n-r,p],{r,0,n-1}]-If[EvenQ[p],0,(-1)^Floor[(n+p)/2](t-Pi)^(n+p)/(n+p)!])];
CC[n_,p_]:=(-1)^(n+p-1)Sum[(-1)^k(1-(-1)^n If[k==0,1,0])Binomial[n+k-1,k](-1)^j Pi^(2j)/(2j)!NielsenS1[n+k-2j,p-k],{k,0,p-1},{j,0,(n+k-1)/2}]+(-1)^(n+Floor[(n+p)/2])If[OddQ[n+p],0,1]Pi^(n+p)/(n+p)/(n-1)!/p!;


NielsenLiMReduced[n_,p_]:=1/2(-Sum[(-1)^s Binomial[n+s-1,s]Li[Join[{n+s+1},Table[1,{p-s-1}]],-1],{s,1,p-1}]+(-1)^p CC[n,p]);


NielsenClGlReduce[expr_]:=Module[{x=expr},
	x=x/.Cl[L_?NielsenQ,t_]/;(EvenQ[Length[L]]):>NielsenClGlReduced[L[[1]]-1,Length[L],t];
	x=x/.Gl[L_?NielsenQ,t_]/;(OddQ[Length[L]]):>NielsenClGlReduced[L[[1]]-1,Length[L],t];
	Expand[x]
];


(* ::Subsection:: *)
(*Reduction at pi/3*)


(* ::Text:: *)
(*Very preliminary (rough and dirty) implementation:*)


KoelbigIdentity[a_,b_]:=NielsenS1[a+1,b+1]+(-1)^a(I Pi/3)^(a+b+2)/(a+1)!/(b+1)!-Sum[(-I Pi/3)^r/r!Li[Join[{a+2-r},Table[1,{b}]],Exp[I Pi/3]],{r,0,a}]-Conjugate[Sum[(-I Pi/3)^r/r!Li[Join[{b+2-r},Table[1,{a}]],Exp[I Pi/3]],{r,0,b}]];


NielsenClMuReduced[n_,p_]:=Module[{li},
	li=Cl[Join[{n+1},Table[1,{p-1}]],Pi/3];
	li/.Solve[If[EvenQ[n+p],Im,Re][KoelbigIdentity[n-1,p-1]]==0,li]//First];
NielsenGlMuReduced[n_,p_]:=Module[{li},
	li=Gl[Join[{n+1},Table[1,{p-1}]],Pi/3];
	li/.Solve[If[OddQ[n+p],Im,Re][KoelbigIdentity[n-1,p-1]]==0,li]//First];


NielsenClGlReduceMu[expr_]:=Module[{x=expr},
	x=x/.Cl[L_?NielsenQ,Pi/3]/;(Length[L]>L[[1]]-1):>NielsenClMuReduced[L[[1]]-1,Length[L]];
	x=x/.Gl[L_?NielsenQ,Pi/3]/;(Length[L]>=L[[1]]-1):>NielsenGlMuReduced[L[[1]]-1,Length[L]];
	Expand[x]
];


(* ::Text:: *)
(*At pi/3 there are extra relations. The following is a table of these relations up to weight 16.*)


MCVReductions={Cl[{4,1,1},\[Pi]/3]->-(1/18) \[Pi]^2 Cl[{4},\[Pi]/3]+3 Cl[{6},\[Pi]/3]-11/324 \[Pi]^3 Zeta[3]-29/108 \[Pi] Zeta[5],Cl[{5,1,1},\[Pi]/3]->-(2/3) \[Pi] Cl[{6},\[Pi]/3]-(17 \[Pi]^4 Zeta[3])/4860-25/972 \[Pi]^2 Zeta[5]+(3229 Zeta[7])/1296,Cl[{7,1,1},\[Pi]/3]->-\[Pi] Cl[{8},\[Pi]/3]-(5423 \[Pi]^6 Zeta[3])/7348320+Zeta[3]^3/18-(53 \[Pi]^4 Zeta[5])/6480-(1285 \[Pi]^2 Zeta[7])/23328+(37729 Zeta[9])/7776,Cl[{6,1,1,1,1},\[Pi]/3]->-((5 \[Pi]^4 Cl[{6},\[Pi]/3])/1944)+2/9 \[Pi]^2 Cl[{8},\[Pi]/3]-21 Cl[{10},\[Pi]/3]-1/18 \[Pi]^2 Cl[{6,1,1},\[Pi]/3]+4 Cl[{8,1,1},\[Pi]/3]+(10441 \[Pi]^7 Zeta[3])/7348320-5/54 \[Pi] Zeta[3]^3+(67 \[Pi]^5 Zeta[5])/3888+(41417 \[Pi]^3 Zeta[7])/209952+(15241 \[Pi] Zeta[9])/209952,Cl[{7,1,1,1,1},\[Pi]/3]->-(1/27) \[Pi]^3 Cl[{8},\[Pi]/3]+8 \[Pi] Cl[{10},\[Pi]/3]-\[Pi] Cl[{8,1,1},\[Pi]/3]+5/2 Cl[{9,1,1},\[Pi]/3]-(20143 \[Pi]^8 Zeta[3])/44089920-5/648 \[Pi]^2 Zeta[3]^3-(2851 \[Pi]^6 Zeta[5])/489888+5/4 Zeta[3]^2 Zeta[5]-(202835 \[Pi]^4 Zeta[7])/2519424-(325997 \[Pi]^2 Zeta[9])/629856+(422407 Zeta[11])/279936,Cl[{8,1,1,1,1},\[Pi]/3]->-((5 \[Pi]^4 Cl[{8},\[Pi]/3])/1944)+2 \[Pi]^2 Cl[{10},\[Pi]/3]-75 Cl[{12},\[Pi]/3]-1/18 \[Pi]^2 Cl[{8,1,1},\[Pi]/3]+7/6 \[Pi] Cl[{9,1,1},\[Pi]/3]+8 Cl[{10,1,1},\[Pi]/3]+(280843 \[Pi]^9 Zeta[3])/558013050-(31 \[Pi]^3 Zeta[3]^3)/1944+(59 \[Pi]^7 Zeta[5])/9720-7/12 \[Pi] Zeta[3]^2 Zeta[5]+(479 \[Pi]^5 Zeta[7])/7290+(7106723 \[Pi]^3 Zeta[9])/11337408-(40263227 \[Pi] Zeta[11])/5038848,Cl[{9,1,1,1,1},\[Pi]/3]->-(10/81) \[Pi]^3 Cl[{10},\[Pi]/3]+65/3 \[Pi] Cl[{12},\[Pi]/3]-1/18 \[Pi]^2 Cl[{9,1,1},\[Pi]/3]-4/3 \[Pi] Cl[{10,1,1},\[Pi]/3]+7 Cl[{11,1,1},\[Pi]/3]-(205917611 \[Pi]^10 Zeta[3])/4714094246400-(25 \[Pi]^4 Zeta[3]^3)/11664-(444637 \[Pi]^8 Zeta[5])/881798400-1/36 \[Pi]^2 Zeta[3]^2 Zeta[5]+3/2 Zeta[3] Zeta[5]^2-(49183 \[Pi]^6 Zeta[7])/7348320+3/2 Zeta[3]^2 Zeta[7]-(35468449 \[Pi]^4 Zeta[9])/340122240-(20263021 \[Pi]^2 Zeta[11])/90699264-(1828525379 Zeta[13])/60466176,Cl[{11,1,1,1,1},\[Pi]/3]->-(25/162) \[Pi]^3 Cl[{12},\[Pi]/3]+130/3 \[Pi] Cl[{14},\[Pi]/3]-1/18 \[Pi]^2 Cl[{11,1,1},\[Pi]/3]-5/3 \[Pi] Cl[{12,1,1},\[Pi]/3]+25/2 Cl[{13,1,1},\[Pi]/3]-(455303841241 \[Pi]^12 Zeta[3])/231650591268096000-(2585 \[Pi]^6 Zeta[3]^3)/8817984+Zeta[3]^5/360-(8585987 \[Pi]^10 Zeta[5])/523788249600-(341 \[Pi]^4 Zeta[3]^2 Zeta[5])/38880-1/24 \[Pi]^2 Zeta[3] Zeta[5]^2+(7 Zeta[5]^3)/12-(611083 \[Pi]^8 Zeta[7])/2645395200-1/24 \[Pi]^2 Zeta[3]^2 Zeta[7]+7/2 Zeta[3] Zeta[5] Zeta[7]-(55637 \[Pi]^6 Zeta[9])/11022480+65/36 Zeta[3]^2 Zeta[9]-(1100975773 \[Pi]^4 Zeta[11])/9795520512-(3937699 \[Pi]^2 Zeta[13])/40310784-(390843481583 Zeta[15])/4081466880,Cl[{8,1,1,1,1,1,1},\[Pi]/3]->(7 \[Pi]^6 Cl[{8},\[Pi]/3])/262440-37/972 \[Pi]^4 Cl[{10},\[Pi]/3]-106/9 \[Pi]^2 Cl[{12},\[Pi]/3]+693 Cl[{14},\[Pi]/3]+(\[Pi]^4 Cl[{8,1,1},\[Pi]/3])/1944-7/324 \[Pi]^3 Cl[{9,1,1},\[Pi]/3]-1/9 \[Pi]^2 Cl[{10,1,1},\[Pi]/3]-47/6 \[Pi] Cl[{11,1,1},\[Pi]/3]-52 Cl[{12,1,1},\[Pi]/3]+6 Cl[{10,1,1,1,1},\[Pi]/3]-(11626062137 \[Pi]^11 Zeta[3])/28284565478400+(2927 \[Pi]^5 Zeta[3]^3)/349920-(232817111 \[Pi]^9 Zeta[5])/47617113600+22/81 \[Pi]^3 Zeta[3]^2 Zeta[5]+11/4 \[Pi] Zeta[3] Zeta[5]^2-(2379829 \[Pi]^7 Zeta[7])/44089920+11/4 \[Pi] Zeta[3]^2 Zeta[7]-(554756773 \[Pi]^5 Zeta[9])/1020366720-(826373575 \[Pi]^3 Zeta[11])/181398528+(26764061837 \[Pi] Zeta[13])/362797056,Cl[{9,1,1,1,1,1,1},\[Pi]/3]->-((11 \[Pi]^5 Cl[{10},\[Pi]/3])/7290)+55/18 \[Pi]^3 Cl[{12},\[Pi]/3]-118 \[Pi] Cl[{14},\[Pi]/3]+(\[Pi]^4 Cl[{9,1,1},\[Pi]/3])/1944-4/81 \[Pi]^3 Cl[{10,1,1},\[Pi]/3]+29/18 \[Pi]^2 Cl[{11,1,1},\[Pi]/3]+10 \[Pi] Cl[{12,1,1},\[Pi]/3]+7 Cl[{13,1,1},\[Pi]/3]-4/3 \[Pi] Cl[{10,1,1,1,1},\[Pi]/3]+(11541670393 \[Pi]^12 Zeta[3])/188028077328000-(829 \[Pi]^6 Zeta[3]^3)/314928+(7 Zeta[3]^5)/180+(1216275073 \[Pi]^10 Zeta[5])/1571364748800-(1873 \[Pi]^4 Zeta[3]^2 Zeta[5])/19440-\[Pi]^2 Zeta[3] Zeta[5]^2+(7 Zeta[5]^3)/2+(150755 \[Pi]^8 Zeta[7])/17635968-\[Pi]^2 Zeta[3]^2 Zeta[7]+21 Zeta[3] Zeta[5] Zeta[7]+(706963291 \[Pi]^6 Zeta[9])/9183300480+203/18 Zeta[3]^2 Zeta[9]+(3647612837 \[Pi]^4 Zeta[11])/8162933760-(23048449031 \[Pi]^2 Zeta[13])/1088391168+(1136684170913 Zeta[15])/10883911680,Cl[{10,1,1,1,1,1,1},\[Pi]/3]->-((61 \[Pi]^6 Cl[{10},\[Pi]/3])/524880)+(235 \[Pi]^4 Cl[{12},\[Pi]/3])/1944-347/9 \[Pi]^2 Cl[{14},\[Pi]/3]+2730 Cl[{16},\[Pi]/3]-(5 \[Pi]^4 Cl[{10,1,1},\[Pi]/3])/1944+1/18 \[Pi]^3 Cl[{11,1,1},\[Pi]/3]+11/18 \[Pi]^2 Cl[{12,1,1},\[Pi]/3]-91/6 \[Pi] Cl[{13,1,1},\[Pi]/3]-148 Cl[{14,1,1},\[Pi]/3]-1/18 \[Pi]^2 Cl[{10,1,1,1,1},\[Pi]/3]+11 Cl[{12,1,1,1,1},\[Pi]/3]-(11932186349231 \[Pi]^13 Zeta[3])/99278824829184000+(7789 \[Pi]^7 Zeta[3]^3)/4898880-1/120 \[Pi] Zeta[3]^5-(237224761 \[Pi]^11 Zeta[5])/168360508800+(5929 \[Pi]^5 Zeta[3]^2 Zeta[5])/116640+13/24 \[Pi]^3 Zeta[3] Zeta[5]^2+73/36 \[Pi] Zeta[5]^3-(53208533 \[Pi]^9 Zeta[7])/3401222400+13/24 \[Pi]^3 Zeta[3]^2 Zeta[7]+73/6 \[Pi] Zeta[3] Zeta[5] Zeta[7]-(171641 \[Pi]^7 Zeta[9])/1049760+71/12 \[Pi] Zeta[3]^2 Zeta[9]-(1462574791 \[Pi]^5 Zeta[11])/906992640-(19008737965 \[Pi]^3 Zeta[13])/1224440064+(4559169750899 \[Pi] Zeta[15])/24488801280};
(* MCV: weight 8 *)
MGVReductions={Gl[{5,1},\[Pi]/3]->(209 \[Pi]^6)/918540-Zeta[3]^2/6,Gl[{6,1,1,1},\[Pi]/3]->-((645259 \[Pi]^9)/4081466880)-1/18 \[Pi]^2 Gl[{6,1},\[Pi]/3]+5/6 \[Pi] Gl[{7,1},\[Pi]/3]+5 Gl[{8,1},\[Pi]/3]+7/216 \[Pi]^3 Zeta[3]^2+5/6 \[Pi] Zeta[3] Zeta[5],Gl[{7,1,1,1},\[Pi]/3]->(174605509 \[Pi]^10)/4714094246400-1/18 \[Pi]^2 Gl[{7,1},\[Pi]/3]-\[Pi] Gl[{8,1},\[Pi]/3]+7/2 Gl[{9,1},\[Pi]/3]+(53 \[Pi]^4 Zeta[3]^2)/12960+1/36 \[Pi]^2 Zeta[3] Zeta[5]-(5 Zeta[5]^2)/4-5/2 Zeta[3] Zeta[7],Gl[{9,1,1,1},\[Pi]/3]->(17014205987 \[Pi]^12)/4454819062848000-1/18 \[Pi]^2 Gl[{9,1},\[Pi]/3]-4/3 \[Pi] Gl[{10,1},\[Pi]/3]+8 Gl[{11,1},\[Pi]/3]+(1529 \[Pi]^6 Zeta[3]^2)/2449440-Zeta[3]^4/72+(25 \[Pi]^4 Zeta[3] Zeta[5])/1944+1/36 \[Pi]^2 Zeta[5]^2+1/18 \[Pi]^2 Zeta[3] Zeta[7]-3 Zeta[5] Zeta[7]-28/9 Zeta[3] Zeta[9],Gl[{8,1,1,1,1,1},\[Pi]/3]->(260777767519 \[Pi]^13)/7942305986334720-(5 \[Pi]^4 Gl[{8,1},\[Pi]/3])/1944+7/162 \[Pi]^3 Gl[{9,1},\[Pi]/3]+7/18 \[Pi]^2 Gl[{10,1},\[Pi]/3]-49/6 \[Pi] Gl[{11,1},\[Pi]/3]-63 Gl[{12,1},\[Pi]/3]-1/18 \[Pi]^2 Gl[{8,1,1,1},\[Pi]/3]+7 Gl[{10,1,1,1},\[Pi]/3]-(13619 \[Pi]^7 Zeta[3]^2)/6298560+7/216 \[Pi] Zeta[3]^4-(929 \[Pi]^5 Zeta[3] Zeta[5])/19440-55/216 \[Pi]^3 Zeta[5]^2-55/108 \[Pi]^3 Zeta[3] Zeta[7]-35/6 \[Pi] Zeta[5] Zeta[7]-301/54 \[Pi] Zeta[3] Zeta[9],Gl[{9,1,1,1,1,1},\[Pi]/3]->-((417553336429 \[Pi]^14)/77216863756032000)+(\[Pi]^4 Gl[{9,1},\[Pi]/3])/1944-4/81 \[Pi]^3 Gl[{10,1},\[Pi]/3]+17/9 \[Pi]^2 Gl[{11,1},\[Pi]/3]+20 \[Pi] Gl[{12,1},\[Pi]/3]-66 Gl[{13,1},\[Pi]/3]-4/3 \[Pi] Gl[{10,1,1,1},\[Pi]/3]+6 Gl[{11,1,1,1},\[Pi]/3]+(644621 \[Pi]^8 Zeta[3]^2)/2645395200+1/324 \[Pi]^2 Zeta[3]^4+(1703 \[Pi]^6 Zeta[3] Zeta[5])/262440-1/2 Zeta[3]^3 Zeta[5]+(65 \[Pi]^4 Zeta[5]^2)/1296+65/648 \[Pi]^4 Zeta[3] Zeta[7]+14/9 \[Pi]^2 Zeta[5] Zeta[7]-(3 Zeta[7]^2)/2+128/81 \[Pi]^2 Zeta[3] Zeta[9]-4 Zeta[5] Zeta[9]-6 Zeta[3] Zeta[11],Gl[{10,1,1,1,1,1},\[Pi]/3]->(831867323171 \[Pi]^15)/85289535875980800-(5 \[Pi]^4 Gl[{10,1},\[Pi]/3])/1944+5/36 \[Pi]^3 Gl[{11,1},\[Pi]/3]+19/6 \[Pi]^2 Gl[{12,1},\[Pi]/3]-143/4 \[Pi] Gl[{13,1},\[Pi]/3]-165 Gl[{14,1},\[Pi]/3]-1/18 \[Pi]^2 Gl[{10,1,1,1},\[Pi]/3]+3/2 \[Pi] Gl[{11,1,1,1},\[Pi]/3]+12 Gl[{12,1,1,1},\[Pi]/3]-(159470827 \[Pi]^9 Zeta[3]^2)/285702681600+(41 \[Pi]^3 Zeta[3]^4)/7776-(549671 \[Pi]^7 Zeta[3] Zeta[5])/44089920+1/4 \[Pi] Zeta[3]^3 Zeta[5]-(5089 \[Pi]^5 Zeta[5]^2)/77760-(5089 \[Pi]^5 Zeta[3] Zeta[7])/38880-803/648 \[Pi]^3 Zeta[5] Zeta[7]-29/8 \[Pi] Zeta[7]^2-(2327 \[Pi]^3 Zeta[3] Zeta[9])/1944-27/4 \[Pi] Zeta[5] Zeta[9]-23/4 \[Pi] Zeta[3] Zeta[11],Gl[{11,1,1,1,1,1},\[Pi]/3]->-((319670954888762597 \[Pi]^16)/306223549609121464320000)-(5 \[Pi]^4 Gl[{11,1},\[Pi]/3])/1944-25/162 \[Pi]^3 Gl[{12,1},\[Pi]/3]+133/36 \[Pi]^2 Gl[{13,1},\[Pi]/3]+124/3 \[Pi] Gl[{14,1},\[Pi]/3]-182 Gl[{15,1},\[Pi]/3]-1/18 \[Pi]^2 Gl[{11,1,1,1},\[Pi]/3]-5/3 \[Pi] Gl[{12,1,1,1},\[Pi]/3]+23/2 Gl[{13,1,1,1},\[Pi]/3]+(51208819 \[Pi]^10 Zeta[3]^2)/4714094246400+(341 \[Pi]^4 Zeta[3]^4)/466560+(374947 \[Pi]^8 Zeta[3] Zeta[5])/1322697600+1/72 \[Pi]^2 Zeta[3]^3 Zeta[5]+(6089 \[Pi]^6 Zeta[5]^2)/2449440-7/8 Zeta[3]^2 Zeta[5]^2+(6089 \[Pi]^6 Zeta[3] Zeta[7])/1224720-7/12 Zeta[3]^3 Zeta[7]+(191 \[Pi]^4 Zeta[5] Zeta[7])/1944+11/12 \[Pi]^2 Zeta[7]^2+(6071 \[Pi]^4 Zeta[3] Zeta[9])/58320+67/36 \[Pi]^2 Zeta[5] Zeta[9]+35/6 Zeta[7] Zeta[9]+23/12 \[Pi]^2 Zeta[3] Zeta[11]+7/2 Zeta[5] Zeta[11]};
(* MGV: weight 10 *)
LiReductionTable=Join[MCVReductions,MGVReductions];


End[];


EndPackage[];
