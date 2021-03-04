(* ::Package:: *)

BeginPackage["cLTD`"];


ResourceFunction["RandomString"];

cLTD::usage="\
Compute the correpsonding cLTD expression using FORM.
"<>Style["Input Format",Bold]~ToString~StandardForm<>":
The scalar propagator of the expression must be expressed as: 
	-> prop[4momentum,mass] (e.g. prop[k1+p1,0])
Outside of prop[] the presence of a momentum not contracted 
into a dot product will be considered as being its energy component 

"<>Style["Options",Bold]~ToString~StandardForm<>": {\"FORMpath\"\[Rule]\"form\",\"WorkingDirectory\"\[Rule]Directory[]}

"<>Style["Example",Bold]~ToString~StandardForm<>":
 expr = k1 prop[k1-p1,0]prop[k1-p2,0] - k1.p1 p1.p2 p1 prop[k1-p3,m]prop[k1-p4,m] prop[k1-p5,m]
=> cLTD[expr,\"/home/andrea/BitBucket/alphaloop/\",\"FORMpath\"\[Rule]\"form\"]

WARNING: when summing multiple diagrams they must contain the same 
         number of loops, if this is not possible then call cLTD[] 
         for each of them individually.

"<>Style["Output",Bold]~ToString~StandardForm<>":
  { cLTD expression, Energies definition} 
  The expression contains the den[] and cLTDnorm[] function which are to be evaluated as 
		den[a_] :> 1/a, cLTDnorm[a_] :> 1/a
";

Print[cLTD::usage];

toLTDprop=Module[{SP},
	SetAttributes[SP,Orderless];
	SP[p1_+p2_,p3_]:=SP[p1,p3]+SP[p2,p3];
	SP[-p1_,p2_]:=-SP[p1,p2];
	prop[mom_,m_]:>LTDprop[mom,Sqrt[SP[mom,mom]-m^2]/.SP->Dot]];
(*default loop momenta*)
{k0,k1,k2,k3,k4};
cLTDnorm;


Begin["cLTDPrivate`"];


FORMfile[filename_,expr_,nLoops_,alPATH_]:=Module[{
file="#--\n"},
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf.frm # partial_fractioning_vars\n";
file = file <> "Auto S sp,q;\n";
file = file <> "#define LOOPS \""<>ToString[nLoops]<>"\"\n";
file = file <> "#define topoID \""<>filename<>"\"\n";
file = file <> "L F = "<>expr<>";\n";
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf.frm # partial_fractioning\n.sort\n";
file = file <> "L firstF = firstterm_(F);\n.sort\n#if ( termsin(firstF) != 0 )\n";
file = file <> "#redefine oldextrasymbols \"`extrasymbols_'\"
ExtraSymbols, underscore, den;
hide firstF;
argtoextrasymbol den;
id den(y?) = y;
.sort:den;\n#redefine newextrasymbols \"`extrasymbols_'\"\noff statistics;\n";
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf_numerator.frm # pf_num\non statistics;\n"; 
file = file <> "#if `oldextrasymbols' != `newextrasymbols'\n";
file = file <> "multiply replace_(<den{`oldextrasymbols'+1}_,den(extrasymbol_({`oldextrasymbols'+1}))>\\
                  ,...,<den`extrasymbols_'_,den(extrasymbol_(`extrasymbols_'))>);\n";
file = file <> "#endif\n";
file = file <> "#endif\n.sort\nFormat mathematica;\n#write<"<>filename<>".out> \"%E\" F;\n.end";
file
]


Complement[{k1,k2,k3},{c1,k1}]


Options[cLTD]={"loopmom"->{k0,k1,k2,k3},"FORMpath"->"form","WorkingDirectory"->Directory[]};
cLTD[expression_,alPATH_,OptionsPattern[]]:=Module[{expr=If[Head[expression]===Plus,List@@expression, {expression}]/.toLTDprop,
energies,i=0,loop0subs,spsubs,FORMinput, 
filename=ResourceFunction["RandomString"][10],
runfilename,cLTDfilename,
result, return,cleanKs},
(*Sanitisation*)
cleanKs=Table[k->ToExpression["nonLoopk"<>ToString[i++]],{k, Complement[Table[ToExpression["k"<>ToString[n]],{n,0,9}],OptionValue["loopmom"]]}];i=0;

(*map FORM symbols to input values*)
energies = Table[ToExpression["E"<>ToString[i++]]->e,{e,Last/@Union[Cases[expr,LTDprop[__],Infinity]]}];i=0;
loop0subs = Table[l->ToExpression["LTDk"<>ToString[i++]],{l,Select[OptionValue["loopmom"],!FreeQ[expr,#]&]}];i=0;
If[Length[loop0subs] == 0, Print["No loop momenta"]; Abort[]];
expr = expr/.ReplaceAll[energies,Rule[a_,b_]:>Rule[b,a]]/.cleanKs;

(*Check for scalar products in the numerator*)
spsubs = Table[x->ToExpression["sp"<>ToString[i++]],{x,Union[Cases[expr,_Dot,Infinity]]}];i=0;
expr = expr/.spsubs/.loop0subs;

(*Create strings to send to FORM*)
FORMinput=StringReplace[{"LTD"->"","[":>"(","]"->")"}][ToString[Plus@@expr,InputForm]];
(*Print[FORMinput];*)

(*Call FORM*)
runfilename=OptionValue["WorkingDirectory"]<>"/cLTD"<>filename<>".frm";
cLTDfilename=OptionValue["WorkingDirectory"]<>"/"<>filename<>".out";
(*Print[runfilename];*)
Export[runfilename,FORMfile[filename,FORMinput,Length[loop0subs],alPATH],"Text"];
return = RunProcess[{OptionValue["FORMpath"],runfilename},ProcessDirectory->OptionValue["WorkingDirectory"]];
If[return["ExitCode"] != 0, Print[return]];
result=ToExpression[StringReplace[" "|"\\"|"\n"->""][Import[cLTDfilename,"Text"]]];
DeleteFile[{runfilename,cLTDfilename}];
result = result/.ReplaceAll[spsubs,Rule[a_,b_]:>Rule[b,a]]
				/.ReplaceAll[cleanKs,Rule[a_,b_]:>Rule[b,a]]
				/.ToExpression["norm"]->ToExpression["cLTDnorm"];
{Collect[result,_cLTDnorm], energies}
]


End[];
EndPackage[];


(*cLTD[ prop[k0-p1,0]prop[k1-p2,0],"/home/andrea/BitBucket/alphaloop/","loopmom"\[Rule]{k0}]*)


(*cLTD[ c1 c1.k1 prop[c1-p1,0]prop[c1-p2,0],"/home/andrea/BitBucket/alphaloop/","loopmom"\[Rule]{c1}]*)


(*cLTD[k1 prop[k1-p1,0]prop[k1-p2,0] - k1.p1 p1.p2 p1 prop[k1-p3,m]prop[k1-p4,m] prop[k1-p5,m],"/home/andrea/BitBucket/alphaloop/"]*)
