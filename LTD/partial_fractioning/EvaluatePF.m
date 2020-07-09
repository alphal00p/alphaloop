(* ::Package:: *)

(* ::Section:: *)
(*Evaluate PF Expression*)


(* ::Subsection::Closed:: *)
(*Path*)


$PFPATH = NotebookDirectory[];


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Numerator*)


ClearAll[NumVarUnfold,num]
NumVarUnfold[NumFun_,y_List,yvar_]:=Module[{newy,y0},
	If[Length[y] == 1,
		NumFun/.yvar->y[[1]],
		newy=Drop[y,-2];
		(NumVarUnfold[NumFun,Append[newy,y[[-1]]],yvar]-NumVarUnfold[NumFun,Append[newy,y[[-2]]],yvar])/(y[[-2]]-y[[-1]])]//Cancel
]


num[ys__List,Num_,vars_List]:=Module[{xs,x,vi,numtmp},
	(*Check inputs*)
	If[Length[{ys}]!=Length[vars],Print["ys list must be of same length as vars"];Abort[]];
	
	(*Top numerator N0 initialized with the provided Num function*)
	numtmp=Num@@vars;
	(*Apply steps*)
	Table[
		xs = vars[[i+1;;]]/.List->Sequence;
		numtmp=NumVarUnfold[numtmp,{ys}[[i]],vars[[i]]];
		(*num[i][xs]*),
	{i,Length[{ys}]}];
	(*Return Final expression*)
	numtmp
]



(*MyFun[x_,y_]:=x^2
MyFun[x_,y_]:=f[x,y]
num[{x1+y,x2,x3},{1},MyFun,{x,y}]*)


(* ::Subsubsection::Closed:: *)
(*Evaluate Element*)


ClearAll[evaluate]
evaluate[{factor_,dens_List,zs_List},NumFunction_,vars_]:=Module[{
		denominator = Times@@dens,
		numerator=num @@Join[zs,{NumFunction,vars}]
	},
	factor *numerator/denominator
]


(* ::Subsection:: *)
(*Import and Evaluate*)


(* ::Subsubsection::Closed:: *)
(*Box*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$PFPATH <> "pf_1l_box.m"]/.a_Real:>Rationalize[a];
(*Define Numerator*)
MyNum[x_]:=x^3
(*MyNum[x_]:=1*)
(*Evaluate*)
pf2=evaluate[#,MyNum,{k0}]&/@pf//Simplify//Cancel//MatrixForm


(* ::Subsubsection::Closed:: *)
(*Sunrise*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$PFPATH <> "pf_2l_sunrise.m"]/.a_Real:>Rationalize[a];
(*Define Numerator*)
MyNum[x_,y_]:=c[1]x+c[2]y^2x^2
(*MyNum[x_,y_]:=1*)
(*Evaluate*)
evaluate[#,MyNum,{k0,k1}]&/@pf//MatrixForm//Simplify


(* ::Input:: *)
(**)


(* ::Subsubsection:: *)
(*PentaBox*)


ClearAll[pf, MyNum]
pf = Import[$PFPATH <> "pf_2l_pentabox.m"]/.a_Real:>Rationalize[a]; 
MyNum[x_, y_] := x^2y^2
(*MyNum[x_, y_] := 1*)
Print["Evaluate:"]
Monitor[pf2=Table[evaluate[pf[[ipf]], MyNum, {k0, k1}], {ipf, Length[pf]}], PercentForm[N[ipf/Length[pf]]]]


blocks=Select[Chop[pf2]/.a_Real:>Round[a],#=!=0&]//Cancel
dens=Union[Flatten[Chop[pf[[All,2]]]/.a_Real:>Round[a]]][[2;;]]//Simplify;
subDen=Table[dens[[i]]->Den[i],{i,Length[dens]}];


Plus@@(blocks/.subDen)


(* ::Subsubsection::Closed:: *)
(*2x2 Fishnet*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import[$PFPATH <> "pf_4l_2x2fishnet.m"]/.a_Real:>Rationalize[a];
(*Define Numerator*)
MyNum[x1_,x2_,x3_,x4_]:=1
(*Evaluate*)
Print["Evaluate:"]
Monitor[ipf=0;Table[evaluate[pf[[ipf]],MyNum,{k0,k1,k2,k3}],{ipf,Length[pf]}],PercentForm[N[ipf/Length[pf]]]]



