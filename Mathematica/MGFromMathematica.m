(* ::Package:: *)

Print["
Package for accessing MadGraph Fortran numerators from Mathematica.

Note that a MG5aMC process folder are generated from the MG5aMC alphaLoop plugin, for example with:
./bin/mg5_aMC --mode=alphaloop PLUGIN/alphaloop/examples/epem_jjj.aL

The function: 
MGFM$GetLTDHook[OutputProcessRootPath_,OptionsPattern[
{
DEBUG->False
}

Allows you to generate a hook for that process output, which will automatically be placed in the list:

MGFM$AllHooksStarted

Then you can test that one hook is active with:

MGFM$CheckHookStatus[hook_]

And use them to access information from Rust with the following three API entry points (for now):

MGFM$GetMGNumerator[hook_, ComplexMomenta_, OptionsPattern[DEBUG->False]]
"];


(* ::Text:: *)
(*All `private` variables of this package will start with MGFM$*)


(* ::Text:: *)
(*Set your environment variables below*)


MGFM$AllHooksStarted={};


(* ::Text:: *)
(*Then start the process hook to the Python bindings*)


(* ::Input::Initialization:: *)
MGFM$GetLTDHook[OutputProcessRootPath_,OptionsPattern[
{
DEBUG->False
}
]]:= Block[
{
RunResult,OldEnv,
Arguments={
OutputProcessRootPath<>"/SubProcesses/IO_bindings"
}
},
If[OptionValue[DEBUG],Print[Environment["PATH"]]];
RunResult=If[OptionValue[DEBUG],
RunProcess[Arguments],
AppendTo[MGFM$AllHooksStarted,
StartProcess[Arguments]
];
If[Length[MGFM$AllHooksStarted]>0,
WriteLine[MGFM$AllHooksStarted[[-1]],OutputProcessRootPath<>"/Cards/param_card.dat\n"];
];
MGFM$AllHooksStarted[[-1]]
];
RunResult
]


(* ::Text:: *)
(*We can verify that it is indeed running and waiting to be fed with data*)


(* ::Input::Initialization:: *)
MGFM$CheckHookStatus[Hook_]:=ProcessStatus[Hook]=="Running"


(* ::Text:: *)
(*And it is useful to have an easy way to kill all processes started so far *)


(* ::Input::Initialization:: *)
MGFM$KillAllHooks[]:=Block[{},
For[i=1,i<=Length[MGFM$AllHooksStarted],i++,
If[MGFM$CheckHookStatus[MGFM$AllHooksStarted[[i]]],
KillProcess[MGFM$AllHooksStarted[[i]]];
];
];
MGFM$AllHooksStarted={};
]


(* ::Text:: *)
(*The hook can be debuged by running a command like the one below*)


(* ::Input::Initialization:: *)
MGFM$KillAllHooks[]


MGFM$FormatFloat[float_,OptionsPattern[{NDigits->16}]]:=Block[{},
If[float==0,"0.0D0",
ToString[Format[ToString@NumberForm[float,OptionValue[NDigits],NumberFormat->(Row[{#1,"D",If[#3=="","0",#3]}]&),NumberPadding->{"",""}],OutputForm]]
]
]


(* ::Input::Initialization:: *)
MGFM$SendNumbersToHook[hook_,Numbers_,OptionsPattern[{DEBUG->False}]]:=Block[{StrToSend},
(* Write the real momenta to sys.stdin of the hook *)
StrToSend=StringJoin[
Riffle[
Table[MGFM$FormatFloat[N],{N, Numbers}],
" "]
]<>"\n";
If[OptionValue[DEBUG],Print["Raw intput sent: "<>StrToSend];];
WriteLine[hook,StrToSend];
]


(* ::Input::Initialization:: *)
MGFM$ParseNumbers[StrNumbers_,OptionsPattern[DEBUG->False]]:=Block[{},
If[StringContainsQ[StrNumbers,"ERROR"],
Print["Error received from Rust LTD worker: "<> StrNumbers];
None
,
Table[
Read[
If[ke=="nan",
If[OptionValue[DEBUG],Print["NaN received!"]];
StringToStream["0"],
StringToStream[ke]],
Number],{ke,StringSplit[StrNumbers]}
]
]
]


MGFM$ReadFromHook[hook_]:=Block[{msgFromHook},
(*Print["Waiting for msg..."];*)
msgFromHook=ReadLine[hook];
(*Print[msgFromHook];*)
While[Not[StringStartsQ[msgFromHook," TOMATHEMATICA "]],
msgFromHook=ReadLine[hook];
];
(* Drop the "TOMATHEMATICA " tagging prefix *)
StringDrop[msgFromHook,15]
]


(* ::Text:: *)
(*We can now build  functions using the hook to access the API*)


(* ::Input::Initialization:: *)
MGFM$GetMGNumerator[hook_,ProcID_,DiagIDs_,ComplexMomenta_, OptionsPattern[DEBUG->False]]:=Module[
{RawOutput, iMom,Numerator, ParsedOutput,StrToSend},

(* Write the selected process and diagram IDS sys.stdin of the hook *)
StrToSend=ToString[ProcID]<>" "<>StringJoin[Riffle[
Table[ToString[DiagID],{DiagID, DiagIDs}],
" "]]<>"\n";
If[OptionValue[DEBUG],Print["Raw intput sent: "<>StrToSend];];
WriteLine[hook,StrToSend];

(* Write the real momenta to sys.stdin of the hook *)
For[iMom=1,iMom<=Length[ComplexMomenta],iMom++,
MGFM$SendNumbersToHook[hook,
Join@@Table[{Re[ke],Im[ke]},{ke,ComplexMomenta[[iMom]]}],DEBUG->OptionValue[DEBUG]
];
];

(* Now recover the deformed momenta *)
RawOutput = MGFM$ReadFromHook[hook];
If[OptionValue[DEBUG],Print["Raw output received: "<>RawOutput];];

(* Parse it *)
ParsedOutput =MGFM$ParseNumbers[RawOutput,DEBUG->OptionValue[DEBUG]];
If[ParsedOutput==None,Return[0]];
Numerator=ParsedOutput[[1]]+I ParsedOutput[[2]];

(* Return *)
Numerator
]


(* ::Text:: *)
(*Example usage:*)


(*TT=MGFM$GetLTDHook["/Users/valentin/Documents/MG5/3.0.2.py3/NUMERATORS_epem_a_ddx",DEBUG\[Rule]False]
MGFM$CheckHookStatus[TT]
MGFM$CheckHookStatus[MGFM$AllHooksStarted[[-1]]]
MGFM$GetMGNumerator[MGFM$AllHooksStarted[[-1]],
0,{1,1},
{
{500.0,0.,0.,500.0},
{500.0,0.,0.,-500.0},
{500.0,110.9242844438328,444.8307894881214,-199.5529299308788},
{500.0,-110.9242844438328,-444.8307894881214,199.5529299308788}
},
DEBUG\[Rule]False]
MGFM$KillAllHooks[]*)
