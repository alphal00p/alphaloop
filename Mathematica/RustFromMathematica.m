(* ::Package:: *)

Print["
Package for accessing alphaLoop implementation in Rust.

Variables you may want to overwrite are:

RFM$LTDFolder = <Path>;
RFM$PythonInterpreter = <Path>;
RFM$CUBAPath = <Path>;
RFM$SCSPath=<Path>;
RFM$ECOSPath=<Path>;

The function: 
RFM$GetLTDHook[name_,OptionsPattern[
{
RunMode->'LTD',
HyperparametersPath->'LTD/hyperparameters.yaml',
TopologiesPath->'LTD/topologies.yaml',
AmplitudesPath->'LTD/anplitudes.yaml',
DEBUG->False
}

Allows you to generate a hook, which will automatically be placed in the list:

RFM$AllHooksStarted

Then you can test that one hook is active with:

RFM$CheckHookStatus[hook_]

And use them to access information from Rust with the following three API entry points (for now):

RFM$GetLTDDeformation[hook_, RealMomenta_,OptionsPattern[DEBUG->False]]
RFM$GetCrossSectionDeformation[hook_, CutID_,RealMomenta_,OptionsPattern[DEBUG->False]]
RFM$GetRescaling[hook_, CutID_,RealMomenta_,OptionsPattern[DEBUG->False]]
"];


(* ::Text:: *)
(*All `private` variables of this package will start with RFM$*)


(* ::Text:: *)
(*First set the current working directory to be the notebook directory:*)


SetDirectory[NotebookDirectory[]];


(* ::Input:: *)
(**)


(* ::Text:: *)
(*Set your environment variables below*)


RFM$LTDFolder = NotebookDirectory[]<>"..";
RFM$PythonInterpreter = "/opt/local/bin/python3";
RFM$CUBAPath = "/Users/valentin/Documents/HEP_softs/Cuba-4.2/lib";
RFM$SCSPath="/Users/valentin/Documents/HEP_softs/scs/out";
RFM$ECOSPath="/Users/valentin/Documents/HEP_softs/ecos";
RFM$DYLDPATHS=RFM$CUBAPath<>":"<>RFM$SCSPath<>":"<>RFM$ECOSPath;


RFM$AllHooksStarted={};


(* ::Text:: *)
(*Then start the process hook to the Python bindings*)


(* ::Input::Initialization:: *)
RFM$GetLTDHook[name_,OptionsPattern[
{
RunMode->"LTD",
HyperparametersPath->"LTD/hyperparameters.yaml",
TopologiesPath->"LTD/topologies.yaml",
AmplitudesPath->"LTD/anplitudes.yaml",
DEBUG->False
}
]]:= Block[
{
Arguments={
RFM$PythonInterpreter,
RFM$LTDFolder<>"/"<>"API_accessor.py",
"--name",name,
"--mode",OptionValue[RunMode],
"--hyperparameter-path",OptionValue[HyperparametersPath],
"--topologies-path",OptionValue[TopologiesPath],
"--amplitudes-path",OptionValue[AmplitudesPath]
}
},
If[OptionValue[DEBUG],
RunProcess[Arguments,
ProcessEnvironment-><|
"DYLD_LIBRARY_PATH"->RFM$DYLDPATHS
|>
],
AppendTo[RFM$AllHooksStarted,
StartProcess[Arguments,
ProcessEnvironment-><|
"DYLD_LIBRARY_PATH"->RFM$DYLDPATHS
|>
]
];
RFM$AllHooksStarted[[-1]]
]
]


(* ::Text:: *)
(*We can verify that it is indeed running and waiting to be fed with data*)


(* ::Input::Initialization:: *)
RFM$CheckHookStatus[Hook_]:=ProcessStatus[Hook]=="Running"


(* ::Text:: *)
(*And it is useful to have an easy way to kill all processes started so far *)


(* ::Input::Initialization:: *)
RFM$KillAllHooks[]:=Block[{},
For[i=1,i<=Length[RFM$AllHooksStarted],i++,
If[RFM$CheckHookStatus[RFM$AllHooksStarted[[i]]],
KillProcess[RFM$AllHooksStarted[[i]]];
];
];
RFM$AllHooksStarted={};
]


(* ::Text:: *)
(*The hook can be debuged by running a command like the one below*)


(* ::Input::Initialization:: *)
RFM$KillAllHooks[]


(* ::Input:: *)
(*(*RFM$GetLTDHook["LTD/ee_to_dd_2l_bubble.yaml",RunMode\[Rule]"cross_section",DEBUG\[Rule]True]*)*)


(* ::Input::Initialization:: *)
RFM$SendNumbersToHook[hook_,Prefix_, Numbers_,OptionsPattern[DEBUG->False]]:=Block[{StrToSend},
(* Write the real momenta to sys.stdin of the hook *)
StrToSend=Prefix<>" "<>StringJoin[
Riffle[
Table[ToString[NumberForm[N,Infinity,ExponentFunction->(Null&)]],{N, Numbers}],
" "]
];
If[OptionValue[DEBUG],Print["Raw intput sent: "<>StrToSend];];
WriteLine[hook,StrToSend];
]


(* ::Input::Initialization:: *)
RFM$ParseNumbers[StrNumbers_]:=Block[{},
Table[
Read[
If[ke=="nan",
Print["NaN received!"];
StringToStream["0"],
StringToStream[ke]],
Number],{ke,StringSplit[StrNumbers]}
]
]


(* ::Text:: *)
(*We can now build  functions using the hook to access the API*)


(* ::Input::Initialization:: *)
RFM$GetLTDDeformation[hook_, RealMomenta_,OptionsPattern[DEBUG->False]]:=Module[
{RawOutput, DeformedMomenta, ParsedOutput, Jacobian},

(* Write the real momenta to sys.stdin of the hook *)
RFM$SendNumbersToHook[hook,"get_deformation",
Join@@Table[Table[ke,{ke,k}],{k,RealMomenta}],DEBUG->OptionValue[DEBUG]
];

(* Now recover the deformed momenta *)
RawOutput = ReadLine[hook];
If[OptionValue[DEBUG],Print["Raw output received: "<>RawOutput];];

(* Parse it *)
ParsedOutput =RFM$ParseNumbers[RawOutput];
Jacobian=ParsedOutput[[1]]+I ParsedOutput[[2]];
DeformedMomenta=ArrayReshape[ParsedOutput[[3;;]],{Length[ParsedOutput[[3;;]]]/3,3}];

(* Return *)
<|"Jacobian"->Jacobian,"DeformationVectors"->DeformedMomenta|>
]


(* ::Input::Initialization:: *)
RFM$GetCrossSectionDeformation[hook_, CutID_,RealMomenta_,OptionsPattern[DEBUG->False]]:=Module[
{RawOutput, DeformedMomenta, ParsedOutput, Jacobian},

(* Write the real momenta to sys.stdin of the hook *)
RFM$SendNumbersToHook[hook,"get_deformation",
Prepend[Join@@Table[Table[ke,{ke,k}],{k,RealMomenta}],CutID],DEBUG->OptionValue[DEBUG]
];

(* Now recover the deformed momenta *)
RawOutput = ReadLine[hook];
If[OptionValue[DEBUG],Print["Raw output received: "<>RawOutput];];

(* Parse it *)
ParsedOutput =RFM$ParseNumbers[RawOutput];
DeformedMomenta=ArrayReshape[ParsedOutput,{Length[ParsedOutput]/4,4}];

For[i=1,i<=Length[DeformedMomenta]/2,i++,
For[j=1,j<=4,j++,
DeformedMomenta[[i]][[j]]=DeformedMomenta[[i]][[j]]+I DeformedMomenta[[i+(Length[DeformedMomenta]/2)]][[j]]
];
];
DeformedMomenta=DeformedMomenta[[;;(Length[DeformedMomenta]/2)]];
(* Return *)
<|"DeformedMomenta"->DeformedMomenta|>
]


(* ::Input::Initialization:: *)
RFM$GetRescaling[hook_, CutID_,RealMomenta_,OptionsPattern[DEBUG->False]]:=Module[
{RawOutput, solutions, ParsedOutput, Jacobian, tValues, Jacobians},

(* Write the real momenta to sys.stdin of the hook *)
RFM$SendNumbersToHook[hook,"get_scaling",
Prepend[Join@@Table[Table[ke,{ke,k}],{k,RealMomenta}],CutID],DEBUG->OptionValue[DEBUG]
];

(* Now recover the output *)
RawOutput = ReadLine[hook];
If[OptionValue[DEBUG],Print["Raw output received: "<>RawOutput];];

(* Parse it *)
ParsedOutput =RFM$ParseNumbers[RawOutput];
solutions=ArrayReshape[ParsedOutput,{Length[ParsedOutput]/2,2}];
tValues=Table[sol[[1]],{sol,solutions}];
Jacobians=Table[sol[[2]],{sol,solutions}];

(* Return *)
<|"tSolutions"->tValues,"tJacobians"->Jacobians|>
]
