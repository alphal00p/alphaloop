(* ::Package:: *)

(* ::Section:: *)
(*Box in 4 dimensions with subtractions*)


(* ::Subsection::Closed:: *)
(*Functions*)


(* ::Input:: *)
(*ClearAll[AnalyticBox,AnalyticBub,AnalyticTri]*)
(*AnalyticBox[s_,t_]:=2 Gamma[1+\[Epsilon]]Gamma[1-\[Epsilon]]^2/(\[Epsilon]^2Gamma[1-2\[Epsilon]] t s )( *)
(*(-t)^(-\[Epsilon])Hypergeometric2F1[1,-\[Epsilon],1-\[Epsilon],1+t/s-I \[Delta]]*)
(*+(-s-I \[Delta])^(-\[Epsilon])Hypergeometric2F1[1,-\[Epsilon],1-\[Epsilon],1+s/t+I \[Delta]]*)
(*)*)
(*AnalyticBub[s_]:=(-s)^(-\[Epsilon])Gamma[\[Epsilon]]Gamma[1-\[Epsilon]]^2/(Gamma[2-2\[Epsilon]] )*)
(*AnalyticTri[s_]:=-(-s-I \[Delta])^(-1-\[Epsilon])Gamma[1+\[Epsilon]]Gamma[1-\[Epsilon]]^2/(Gamma[1-2\[Epsilon]] )/\[Epsilon]^2*)


(* ::Input:: *)
(*((1-2 \[Epsilon])AnalyticBub[s])/(s \[Epsilon])(*Triangle!!!!*);*)


(* ::Input:: *)
(*msp[p1_,p2_]:= p1[[1]]*p2[[1]]-Sum[p1[[i]]p2[[i]],{i,2,4}]*)
(*msq[p1_]:=msp[p1,p1]*)


(* ::Input:: *)
(*<<HypExp`*)


(* ::Subsection:: *)
(*Seeds*)


(* ::Subsubsection::Closed:: *)
(*33*)


(* ::Input:: *)
(*(*SEED 33*)*)
(*p1={   -3.9756551766525951`17*10^+02    ,0 , 0 , -3.9756551766525951`17*10^+02 };p2={   -3.9756551766525951`17*10^+02    ,0,  0 ,    3.9756551766525951`17*10^+02 };p3={     3.9756551766525951`17*10^+02   ,-5.2461217332126608`17*10^+01   ,  3.0293912899098285`17*10^+02   ,  2.5205960731275806`17*10^+02 };p4={     3.9756551766525951`17*10^+02  ,   5.2461217332126608`17*10^+01  , -3.0293912899098285`17*10^+02   ,-2.5205960731275806`17*10^+02};*)


(* ::Subsubsection::Closed:: *)
(*2*)


(* ::Input:: *)
(*(*SEED 2*)*)
(*p1={-4.8678217388064405`17*10^+02,0,0,-4.8678217388064405`17*10^+02};*)
(*p2={-4.8678217388064405`17*10^+02,0,0, 4.8678217388064405`17*10^+02};*)
(*p3={ 4.8678217388064405`17*10^+02, 1.9365322696179936`17*10^+02 ,1.1431607376733305`17*10^+02,-4.3172577844468481`17*10^+02};*)
(*p4={ 4.8678217388064405`17*10^+02,-1.9365322696179936`17*10^+02,-1.1431607376733305`17*10^+02, 4.3172577844468481`17*10^+02};*)


(* ::Subsubsection:: *)
(*angle*)


(* ::Input:: *)
(*(*SEED 0*)*)
(*p1={-.5,0,0,-.5};*)
(*p2={-.5,0,0, +.5};*)
(*p3={ .5, 0 ,-.5Sin[\[Theta]],-.5 Cos[\[Theta]]};*)
(*p4={ .5,0,.5Sin[\[Theta]], .5Cos[\[Theta]]};*)


(* ::Subsection:: *)
(*Evaluate*)


(* ::Input:: *)
(*ClearAll[B1Lsub,ExpSub]*)
(*ExpSub[x_]:=I \[Pi]^(d/2)/(I \[Pi]^2)(AnalyticBox[s,t]-2  AnalyticTri[t]/(s)-2 AnalyticTri[s]/t)/.d->4-2\[Epsilon]/.{s->msq[p1+p2],t->msq[p2+p3],u->msq[p1+p3]}/.\[Theta]->x//Simplify;*)
(*B1Lsub[\[Theta]_]:=B1Lsub[\[Theta]]=SeriesCoefficient[ExpSub[\[Theta]],{\[Epsilon],0,0}];*)


(* ::Input:: *)
(*Limit[B1Lsub[.5],\[Delta]->0,Direction->"FromAbove"]*)


(* ::Input:: *)
(*Table[{\[Pi]/i,Limit[B1Lsub[\[Pi]/i],\[Delta]->0,Direction->"FromAbove"]},{i,1,10}]//N*)
