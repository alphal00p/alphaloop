(* ::Package:: *)

BeginPackage["cLTD`"]
cLTD::usage="\
Compute the correpsonding cLTD expression using FORM.
"<>Style["Input Format",Bold]~ToString~StandardForm<>":
The scalar propagator of the expression must be expressed as: 
	-> prop[4momentum,mass] (e.g. prop[k1+p1,0])
Outside of prop[] the presence of a loop momentum will be considered as being 
its energy component 

"<>Style["Options",Bold]~ToString~StandardForm<>": {\"FORMpath\"\[Rule]\"form\",\"WorkingDirectory\"\[Rule]Directory[]}

"<>Style["Example",Bold]~ToString~StandardForm<>"
 expr = k1 prop[k1-p1,0]prop[k1-p2,0] - prop[k1-p3,m]prop[k1-p4,m]
=> cLTD[expr,\"/home/andrea/BitBucket/alphaloop/\",\"FORMpath\"\[Rule]\"form\"]

WARNING: when summing multiple diagrams they must contain the same number of loops,
		if this is not possible then call cLTD[] for each of them individually.
";
Options[cLTD]={"loopmom"->{k0,k1,k2,k3},"FORMpath"->"form","WorkingDirectory"->Directory[]};
toLTDprop=Module[{SP},
	SetAttributes[SP,Orderless];
	SP[p1_+p2_,p3_]:=SP[p1,p3]+SP[p2,p3];
	SP[-p1_,p2_]:=-SP[p1,p2];
	prop[mom_,m_]:>LTDprop[mom,Sqrt[SP[mom,mom]-m^2]/.SP->Dot]];
Begin["Private`"]


FORMfile[filename_,expr_,nLoops_,alPATH_]:=Module[{
file="#--\n"},
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf.frm # partial_fractioning_vars\n";
file = file <> "#define LOOPS \""<>ToString[nLoops]<>"\"\n";
file = file <> "#define topoID \""<>filename<>"\"\n";
file = file <> "L F = "<>expr<>";\n";
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf.frm # partial_fractioning\n.sort\n";
file = file <> "#redefine oldextrasymbols \"`extrasymbols_'\"
ExtraSymbols, underscore, den;
argtoextrasymbol den;
id den(y?) = y;
.sort:den;\noff statistics;\n";
file = file <> "#include "<>ToString[alPATH]<>"alpha_loop/pf_numerator.frm # pf_num\non statistics;\n"; 
file = file <> "multiply replace_(<den{`oldextrasymbols'+1}_,den(extrasymbol_({`oldextrasymbols'+1}))>\\
                  ,...,<den`extrasymbols_'_,den(extrasymbol_(`extrasymbols_'))>);";
file = file <> ".sort\nFormat mathematica;\n#write<"<>filename<>".out> \"%E\" F;\n.end";
file
]


cLTD[expression_,alPATH_,OptionsPattern[]]:=Module[{expr=If[Head[expression]===Plus,List@@expression, {expression}]/.toLTDprop,
energies,i=0,loopsubs,FORMinput, 
filename=ResourceFunction["RandomString"][10],
runfilename,cLTDfilename,
result},
(*map FORM symbols to input values*)
energies = Table[ToExpression["E"<>ToString[i++]]->e,{e,Last/@Union[Cases[expr,LTDprop[__],Infinity]]}];i=0;
loopsubs = Table[ToExpression["k"<>ToString[i++]]->l,{l,Select[OptionValue["loopmom"],!FreeQ[expr,#]&]}];
(*Create strings to send to FORM*)
expr = expr/.ReplaceAll[energies,Rule[a_,b_]:>Rule[b,a]]/.ReplaceAll[loopsubs,Rule[a_,b_]:>Rule[b,a]];
FORMinput=StringReplace[{"LTD"->"","[":>"(","]"->")"}][ToString[Plus@@expr,InputForm]];
(*Call FORM*)
runfilename=OptionValue["WorkingDirectory"]<>"/cLTD"<>filename<>".frm";
cLTDfilename=OptionValue["WorkingDirectory"]<>"/"<>filename<>".out";
(*Print[runfilename];*)
Export[runfilename,FORMfile[filename,FORMinput,Length[loopsubs],alPATH],"Text"];
RunProcess[{OptionValue["FORMpath"],runfilename},ProcessDirectory->OptionValue["WorkingDirectory"]];
result=ToExpression[StringReplace[" "|"\\"|"\n"->""][Import[cLTDfilename,"Text"]]];
DeleteFile[{runfilename,cLTDfilename}];
{result, energies}
]


End[]
EndPackage[]


(*cLTD[k1 prop[k1-p1,0]prop[k1-p2,0],"/home/andrea/BitBucket/alphaloop/"]*)
