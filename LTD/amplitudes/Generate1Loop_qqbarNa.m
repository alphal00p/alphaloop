(* ::Package:: *)

PhotonNumber = 2;


(* ::Subsection::Closed:: *)
(*Functions*)


(* ::Subsubsection::Closed:: *)
(*SPs*)


ClearAll[SP,SQ]
SetAttributes[SP,Orderless]
SP[p1_+p2_,q_]:=SP[p1,q]+SP[p2,q]
SP[p_]:=SP[p,p]
SP[p[i_],p[j_]]:=sij[i,j]/2(*-sij[i,i]/2-sij[j,j]/2*)
SP[p[i_],p[i_] ]:=0
SP[0,_]:=0

sij["\"s"<>ToString[i]<>ToString[j]<>"\""]


(* ::Subsubsection::Closed:: *)
(*Left Loop Reduce*)


(* ::Input:: *)
(*(*Take the current loop and start to shrink it from the right*)*)


LeftLoopReduce[loop_List]:=Module[{
size=Length[Cases[loop,{___,-1,a___,-1,___}:>a,Infinity]]+1,
rule},
rule={factor_,den_,{head___,-1,k[n1_],j[n2_],k[n3_],body___,-1,tail___}}:>{factor/SQ[q[n3]],den /SQ[k[n1]],{head,j[n2],q[n3],-1,k[n3],body,-1,tail}};
Which[
size==2 && !FreeQ[Position[loop[[-1]],-1],1|Length[loop[[-1]]]],{},
size==2 ,{loop},
True,Join[{loop},LeftLoopReduce[loop/.rule]]
]
]


(* ::Subsubsection::Closed:: *)
(*Right Loop Reduce*)


(* ::Input:: *)
(*(*Take the current loop and start to shrink it from the left*)*)


RightLoopReduce[loop_List]:=Module[{
size=Length[Cases[loop,{___,-1,a___,-1,___}:>a,Infinity]]+1,
rule},
rule={factor_,den_,{head___,-1,body___,k[n1_],j[n2_],k[n3_],-1,tail___}}:>{factor/SQ[q[n1]],den /SQ[k[n3]],{head,-1,body,k[n1],-1,q[n1],j[n2],tail}};
Which[
size==2 && !FreeQ[Position[loop[[-1]],-1],1|Length[loop[[-1]]]],{},
size==2 ,{loop},
True,Join[{loop},RightLoopReduce[loop/.rule]]
]
]


(* ::Subsubsection::Closed:: *)
(*Get IR*)


GetIR[maxloop_]:=Module[{rule},
rule={factor_,den_,{-1,k[n1_],head___,k[n2_],tail___,k[n3_],-1}}:>{factor/SQ[q[n2]],den /SQ[k[n2]],{-1,k[n1],head,q[n2],tail,k[n3],-1}};
maxloop//.rule
]


GetIR[maxloop]


(* ::Subsection:: *)
(*Generate qqbarNa*)


constEval={CasimirF->4/3,qU->2/3,qD->-1/3,gS->1.2177157847767195,gEW->0.3079537672443688};


(* ::Subsubsection::Closed:: *)
(*diag list*)


(*{factor, loop_chain}*)
Module[{n=2+PhotonNumber},
maxloop={1,Times@@Array[SQ[k[#]]&,n,0],{-1,Transpose[{Array[k,n-2],Array[j,n-2]}],k[n-1],-1}//Flatten};
vectors = Join[Array[k,n,0],Array[q,n-3,2],Array[j,n-2]];
diaglist=Flatten[RightLoopReduce[#]&/@LeftLoopReduce[maxloop],1];
];
diaglist[[All,1]]=diaglist[[All,1]]*(-I)*gEW^PhotonNumber qD^PhotonNumber gS^2 CasimirF;


diaglist//MatrixForm


(* ::Subsubsection::Closed:: *)
(*Kinematics*)


vectors


evalkinematics = Join[
Table[q[n-1] -> Sum[ToExpression["p"<>ToString[m]],{m,2,n}],{n,2+PhotonNumber}]
];
incomingmom=Table[ToExpression["p"<>ToString[m]],{m,1,2}]
outgoingmom=Table[ToExpression["p"<>ToString[m]],{m,3,2+PhotonNumber}]
allOnShell=Table[sp[#,#]&@ToExpression["p"<>ToString[m]]->0,{m,1,2+PhotonNumber}]

startFromp1=ToExpression["p"<>ToString[2+PhotonNumber]]->-Sum[ToExpression["p"<>ToString[m]],{m,1,1+PhotonNumber}];
photonOutgoing=Table[p->-p,{p,outgoingmom}]
evalkinematics//TableForm
evalkinematicsBIS=evalkinematics/.startFromp1/.photonOutgoing;
evalkinematicsBIS//TableForm


(* ::Subsubsection::Closed:: *)
(*IR*)


(* Get IR CT*)
diaglistFORMir=diaglist[[1]]/.k[n_?(MatchQ[#,0|1|1+PhotonNumber]&)]:>kk[n]/.k[n_]:> q[n]/.kk->k;
diaglistFORMir[[1]]=-diaglistFORMir[[1]]/diaglistFORMir[[2]]/.SQ[_k]:>1;
diaglistFORMir[[2]]=diaglistFORMir[[2]]/.SQ[_q]:>1;
diaglistFORMir


(* ::Subsubsection::Closed:: *)
(*UV*)


(*LO correction*)
diaglistFORMuvLO=Select[Join[diaglist,{diaglistFORMir}],Length[#[[2]]]<4&];
diaglistFORMuvLO[[All,2]]=diaglistFORMuvLO[[All,2]]/.SQ[k[n_]]:>SQ[k[0],muv];
diaglistFORMuvLO[[All,3]]=diaglistFORMuvLO[[All,3]]/.k[n_]:>k[0];
diaglistFORMuvLO[[All,1]]=-diaglistFORMuvLO[[All,1]];



(*NLO corrections needed for the bubbles*)
diaglistFORMuvNLO=Select[Join[diaglist,{diaglistFORMir}],Length[#[[2]]]<3&];
diaglistFORMuvNLO[[All,3]]=diaglistFORMuvNLO[[All,3]]/.{head___,-1,k[n_],-1,tail___}:>{head,-1,k[0],q[n],k[0],-1,tail};
diaglistFORMuvNLO[[All,2]]=diaglistFORMuvNLO[[All,2]]*SQ[k[0],muv]/.SQ[k[n_]]:>SQ[k[0],muv]/.k[n_]:>k[0];
diaglistFORMuvNLO[[All,1]]=(-1)^2*diaglistFORMuvNLO[[All,1]];



diaglistFORMuv=Join[diaglistFORMuvLO,diaglistFORMuvNLO]


(* ::Subsubsection::Closed:: *)
(*Formatting Diagrams*)


(* list of diagrams with CTs *)
diaglistFORM=Join[diaglist,{diaglistFORMir}, diaglistFORMuv];


(* Format Prefactor *)
diaglistFORM[[All,1]]=diaglistFORM[[All,1]]/.{SQ[q[i_]]:>SP[q[i],q[i]]};
diaglistFORM[[All,1]]=diaglistFORM[[All,1]]/.evalkinematics/.SP[a_,b_]:>sp[a,b];


(* Format Denominators *)
diaglistFORM[[All,2]]=diaglistFORM[[All,2]]/.{SQ[k[i_]]:>Den[i,0],SQ[k[0],m_]:>Den[0,m]}(*//.Den[a__]Den[b__]\[RuleDelayed]Den@@Sort[{a,b}]*);


(* Format Numerators *)
diaglistFORM[[All,3]]=diaglistFORM[[All,3]]/.{a_?(MatchQ[#,k|q|j]&)[n_]:>lVec[a[n]],n_Integer:>muL[-n]};
diaglistFORM[[All,3]]=diaglistFORM[[All,3]]/.k[n_]:>k+ q[n]/.evalkinematics;
diaglistFORM[[All,3]]=Map[Apply[gam,#]&,diaglistFORM[[All,3]]]//.gam[head___,lVec[j[n_]],tail___]:>cpol[n,muL[n+1]]gam[head,muL[n+1],tail];
diaglistFORM[[All,3]]=diaglistFORM[[All,3]]/.gam[core__]:>vbarSpinor[1,indS[1]]uSpinor[1,indS[2]] gam[indS[1],core,indS[2]];
diaglistFORM[[All,3]]//TableForm;


(* Remap to Armin's conventions *)
diaglistFORM= diaglistFORM/.startFromp1/.allOnShell/.photonOutgoing/.k->k1;
diaglistFORM[[All,3]]=diaglistFORM[[All,1]]*diaglistFORM[[All,3]];
diaglistFORM = diaglistFORM[[All,2;;]];
diaglistFORM [[All,1]]=Flatten/@(List/@ diaglistFORM[[All,1]]/.Times->List);
diaglistFORM//TableForm


(* ::Subsubsection:: *)
(*To JSON*)


denmap={
Power[Den[n_,m_],p_]:>Module[{name}, {
	"loop_signature"->{1},
	"outgoing_signature"-> Coefficient[q[n]/.evalkinematics/.startFromp1/.photonOutgoing,outgoingmom[[;;-2]],1],
	"incoming_signature"-> Coefficient[q[n]/.evalkinematics/.startFromp1/.photonOutgoing,incomingmom,1],
	"mass" -> ToString[m],
	"power"->p,
	"name"->name}],
Den[n_,m_]:>Module[{name}, {
	"loop_signature"->{1},
	"outgoing_signature"-> Coefficient[q[n]/.evalkinematics/.startFromp1/.photonOutgoing,outgoingmom[[;;-2]],1],
	"incoming_signature"-> Coefficient[q[n]/.evalkinematics/.startFromp1/.photonOutgoing,incomingmom,1],
	"mass" -> m,
	"power"->1,
	"name"->name}]};



(*To JSON *)
diaglistJSON = diaglistFORM;
diaglistJSON[[All,1]]="propagators"->#&/@diaglistJSON[[All,1]]/.denmap;
diaglistJSON[[All,1]]=diaglistJSON[[All,1]]/.Table[#[[i]]->"l"<>ToString[i],{i,Length[#]}]&@Union[Cases[diaglistJSON[[All,1]],Rule["name",a_]:>a,Infinity]];

diaglistJSON[[All,2]]="analytic_num"->StringReplace[{"["->"(","]"->")","I"->"i_", " "->""}][ToString[#,InputForm]]&/@diaglistJSON[[All,2]];

diaglistJSON = {"diagram_list"->diaglistJSON};

Export[NotebookDirectory[]<>"/dd"<>ToString[PhotonNumber]<>"A/qqbar"<>ToString[PhotonNumber]<>"a.json", StringReplace["\\"->""][ExportString[diaglistJSON,"JSON"]],"Text"]


(*To JSON per DIAG*)
diaglistJSON = diaglistFORM;
diaglistJSON[[All,1]]="propagators"->#&/@diaglistJSON[[All,1]]/.denmap;
diaglistJSON[[All,1]]=diaglistJSON[[All,1]]/.Table[#[[i]]->"l"<>ToString[i],{i,Length[#]}]&@Union[Cases[diaglistJSON[[All,1]],Rule["name",a_]:>a,Infinity]];

diaglistJSON[[All,2]]="analytic_num"->StringReplace[{"["->"(","]"->")","I"->"i_", " "->""}][ToString[#,InputForm]]&/@diaglistJSON[[All,2]];

Table[
Export[NotebookDirectory[]<>"/dd"<>ToString[PhotonNumber]<>"A/qqbar"<>ToString[PhotonNumber]<>"a_"<>ToString[diagi]<>".json", StringReplace["\\"->""][ExportString[{"diagram_list"->{diaglistJSON[[diagi]]}},"JSON"]],"Text"],
{diagi,Length[diaglistJSON]}]
