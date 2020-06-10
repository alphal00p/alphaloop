(* ::Package:: *)

(* ::Section:: *)
(*Evaluate PF Expression*)


(* ::Subsection:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*Numerator*)


ClearAll[num]
num[n_Integer][y_List,vars___?(Head[#]=!=List&)]:=Module[{
newy},
If[Length[y] == 1,
num[n][y[[1]],vars]
,
newy=Drop[y,-2];
(num[n][Append[newy,y[[-1]]],vars]-num[n][Append[newy,y[[-2]]],vars])/(y[[-2]]-y[[-1]])]//Cancel
]

num[ys__List,Num_,vars_List]:=Module[{
xs,x},
If[Length[{ys}]!=Length[vars],Print["ys list must be of same length as vars"];Abort[]];
num[0][x__?(Head[#]=!=List&)]:=Num[x];
Table[
xs = vars[[i+1;;]]/.List->Sequence;
num[i][x___?(Head[#]=!=List&)]=num[i-1][{ys}[[i]],x];
(*num[i][xs]*),{i,Length[{ys}]}];
num[Length[vars]][xs]
]


(* ::Subsubsection::Closed:: *)
(*Evaluate Element*)


ClearAll[evaluate]
evaluate[{factor_,dens_List,zs_List},NumFunction_,vars_]:=Module[{
denominator = Times@@dens,
numerator=num @@Join[zs,{NumFunction,vars}]},
factor *numerator/denominator
]


(* ::Subsection:: *)
(*Import and Evaluate*)


(* ::Subsubsection::Closed:: *)
(*Box*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import["pf_1l_box.m"];
(*Define Numerator*)
MyNum[x_]:=c[0]+c[1]x +c[2]x^2
MyNum[x_]:=1
(*Evaluate*)
evaluate[#,MyNum,{x}]&/@pf//MatrixForm//Simplify


(* ::Subsubsection::Closed:: *)
(*Sunrise*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import["pf_2l_sunrise.m"];
(*Define Numerator*)
MyNum[x_,y_]:=c[1]x^2+c[2]y^3x +c[3]x^2y^2 
MyNum[x_,y_]:=1
(*Evaluate*)
evaluate[#,MyNum,{x,y}]&/@pf//MatrixForm//Simplify


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*PentaBox*)


ClearAll[pf, MyNum]
pf = Import["pf_2l_pentabox.m"]; 
(*MyNum[x_, y_] := c[1]*x^2 + c[2]*y^3*x + c[3]*x^2*y^2*)
MyNum[x_, y_] := 1
Print["Evaluate:"]
Monitor[Table[evaluate[pf[[ipf]], MyNum, {x, y}], {ipf, Length[pf]}], PercentForm[N[ipf/Length[pf]]]]


(* ::Subsubsection::Closed:: *)
(*2x2 Fishnet*)


ClearAll[pf,MyNum]
(*Import instructions*)
pf= Import["pf_4l_2x2fishnet.m"];
(*Define Numerator*)
MyNum[x1_,x2_,x3_,x4_]:=1
(*Evaluate*)
Print["Evaluate:"]
Monitor[ipf=0;Table[evaluate[pf[[ipf]],MyNum,{x1,x2,x3,x4}],{ipf,Length[pf]}],
PercentForm[N[ipf/Length[pf]]]]
