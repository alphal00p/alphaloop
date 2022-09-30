(* ::Package:: *)

Print["
Package for accessing alphaLoop implementation in Rust.
This V2 version uses the faster approach of StartExternalSession[Python] to talk to alphaLoop.

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
RFM$Parameterize[hook_,LoopIndex_,ECM_,Xs_,OptionsPattern[{DEBUG->False,f128->False}]]
RFM$InvParameterize[hook_,LoopIndex_,ECM_,Momentum_,OptionsPattern[{DEBUG->False,f128->False}]]
RFM$Evaluate[hook_,Momenta_,OptionsPattern[{DEBUG->False,f128->False}]]
RFM$EvaluateCut[hook_,CutID_,scalingFactor_,Momenta_,OptionsPattern[{DEBUG->False,f128->False,Deformation->None,DeformationJacobian->None}]]
RFM$EvaluateIntegrand[hook_,Xs_,OptionsPattern[{DEBUG->False,f128->False}]]

Set your environment variables as follows:

RFM$CUBAPath = \"<PATHTO>/Cuba-4.2/lib\";
RFM$SCSPath=\"<PATHTO>/scs/out\";
RFM$ECOSPath=\"<PATHTO>/ecos\";
RFM$DYLDPATHS=RFM$CUBAPath<>\":\"<>RFM$SCSPath<>\":\"<>RFM$ECOSPath;

and then call:
RFM$SetEnvironment[];

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


RFM$LoadYaml[name_]:=Block[
{WexPath="Mathematica/wex",Arguments},
If[Not[FileExistsQ[name]],
Print["Could not find target specified yaml file at path: "<>name];
Return[];
];
If[Not[FileExistsQ[RFM$LTDFolder<>"/"<>WexPath]],
Print["Could not find the converter binary 'wex' at path: "<>RFM$LTDFolder<>"/"<>WexPath];
Print["Please first install it with 'cargo install wxf-converter' or following instructions on: https://github.com/GalAster/wolfram-exchange-cli"];
Return[];
];
Arguments={
RFM$LTDFolder<>"/"<>WexPath,
"-f","yaml","-t",name
};
RunProcess[Arguments];
If[Not[FileExistsQ[name<>".m"]],
Print["Wex apparently failed the processing of the specified yaml file since output file "<>name<>".m"<>" cannot be found."];
Return[];
];
Get[name<>".m"]
];


(* ::Text:: *)
(*Facility to set proper environment variables from RFM$DYLDPATHS*)


RFM$SetEnvironment[]:=Block[{},
SetEnvironment["DYLD_LIBRARY_PATH"->If[GetEnvironment["DYLD_LIBRARY_PATH"][[2]]===None,"",GetEnvironment["DYLD_LIBRARY_PATH"][[2]]<>":"]<>RFM$DYLDPATHS];
SetEnvironment["LD_LIBRARY_PATH"->If[GetEnvironment["LD_LIBRARY_PATH"][[2]]===None,"",GetEnvironment["LD_LIBRARY_PATH"][[2]]<>":"]<>RFM$DYLDPATHS];
SetEnvironment["LIBRARY_PATH"->If[GetEnvironment["LIBRARY_PATH"][[2]]===None,"",GetEnvironment["LIBRARY_PATH"][[2]]<>":"]<>RFM$DYLDPATHS];
]


(* ::Input::Initialization:: *)
RFM$GetLTDHook[name_,OptionsPattern[
{
RunMode->"LTD",
HyperparametersPath->"LTD/hyperparameters.yaml",
TopologiesPath->"LTD/topologies.yaml",
AmplitudesPath->"LTD/amplitudes.yaml",
MGNumeratorPath->"N/A",
DEBUG->False
}
]]:= Block[
{
Arguments={
RFM$LTDFolder,
name,
OptionValue[RunMode],
OptionValue[HyperparametersPath],
OptionValue[TopologiesPath],
OptionValue[AmplitudesPath]
},
SessionProlog="import os;os.environ['DYLD_LIBRARY_PATH']='"<>
RFM$DYLDPATHS<>":'+os.environ.get('DYLD_LIBRARY_PATH','');"
<>"os.environ['LD_LIBRARY_PATH']='"<>RFM$DYLDPATHS<>":'+os.environ.get('LD_LIBRARY_PATH','');"
<>"os.environ['LIBRARY_PATH']='"<>RFM$DYLDPATHS<>":'+os.environ.get('LIBRARY_PATH','');"},
SessionProlog=If[OptionValue[MGNumeratorPath]!="N/A",SessionProlog<>"os.environ['MG_NUMERATOR_PATH']='"<>
(OptionValue[MGNumeratorPath]<>"/")<>"';",SessionProlog];
If[OptionValue[DEBUG],
Print["Arguments="<>ToString[Arguments,InputForm]];
Print["Prolog="<>SessionProlog];
];
AppendTo[RFM$AllHooksStarted,StartExternalSession[<|
"System"->"Python","Version"->"3","SessionProlog"->SessionProlog
|>]
];
If[OptionValue[DEBUG],Print["External session started."];];
ExternalEvaluate[RFM$AllHooksStarted[[-1]],File[RFM$LTDFolder<>"/"<>"API_accessor.py"]];
If[OptionValue[DEBUG],Print["API accessor evaluated."];];
ExternalEvaluate[RFM$AllHooksStarted[[-1]],<|"Command"->"API_initialise","Arguments"->Arguments|>];
If[OptionValue[DEBUG],Print["API accessor evaluated."];];
RFM$AllHooksStarted[[-1]]
]


(* ::Text:: *)
(*We can verify that it is indeed running and waiting to be fed with data*)


(* ::Input::Initialization:: *)
RFM$CheckHookStatus[Hook_]:=ExternalEvaluate[Hook,<|"Command"->"API_is_alive","Arguments"->{}|>]


(* ::Text:: *)
(*And it is useful to have an easy way to kill all processes started so far *)


(* ::Input::Initialization:: *)
RFM$KillAllHooks[]:=Block[{},
For[i=1,i<=Length[RFM$AllHooksStarted],i++,
Quiet[DeleteObject[RFM$AllHooksStarted[[i]]]];
];
RFM$AllHooksStarted={};
]


(* ::Text:: *)
(*The hook can be debugged by running a command like the one below*)


(* ::Input::Initialization:: *)
RFM$KillAllHooks[]


(* ::Input:: *)
(*(*RFM$GetLTDHook["LTD/ee_to_dd_2l_bubble.yaml",RunMode\[Rule]"cross_section",DEBUG\[Rule]True]*)*)


(* ::Text:: *)
(*We can now build  functions using the hook to access the API*)


(* ::Input:::GetLTDDeformation:: *)
RFM$GetLTDDeformation[hook_, RealMomenta_,OptionsPattern[DEBUG->False]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_get_deformation","Arguments"->{False,RealMomenta,-1,-1}|>];
If[Not[MemberQ[Keys[Res],"jac"]],Res,
(* Return *)
<|"Jacobian"->Res["jac"],"DeformationVectors"->Res["kappas"]|>
]
]


(* ::Input::GetCrossSectionDeformation:: *)
RFM$GetCrossSectionDeformation[hook_, CutID_, scalingFactor_, RealMomenta_,OptionsPattern[{DEBUG->False}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_get_deformation","Arguments"->{False,RealMomenta,CutID,scalingFactor}|>];
If[Not[MemberQ[Keys[Res],"deformed_momenta"]],Res,
(* Return *)
<|"DeformedMomenta"->Res["deformed_momenta"],"DeformationJacobian"->Res["deformation_jacobian"]|>
]
]


(* ::Input::GetRescaling:: *)
RFM$GetRescaling[hook_, CutID_,RealMomenta_,OptionsPattern[DEBUG->False]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_get_scaling","Arguments"->{False,CutID,RealMomenta}|>];
(* Return *)
If[Not[MemberQ[Keys[Res],"solutions"]],Res,
<|"tSolutions"->Table[s[[1]],{s,Res["solutions"]}],"tJacobians"->Table[s[[2]],{s,Res["solutions"]}]|>
]
]


(* ::Input::Parameterize:: *)
RFM$Parameterize[hook_,LoopIndex_,ECM_,Xs_,OptionsPattern[{DEBUG->False,f128->False}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_parameterize","Arguments"->{OptionValue[f128],LoopIndex,ECM,Xs}|>];
If[Not[MemberQ[Keys[Res],"momentum"]],Res,
(* Return *)
<|"Momentum"->Res["momentum"],"Jacobian"->Res["jacobian"]|>
]
]


(* ::Input::InvParameterize:: *)
RFM$InvParameterize[hook_,LoopIndex_,ECM_,Momentum_,OptionsPattern[{DEBUG->False,f128->False}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_inv_parameterize","Arguments"->{OptionValue[f128],LoopIndex,ECM,Momentum}|>];
If[Not[MemberQ[Keys[Res],"xs"]],Res,
(* Return *)
<|"Xs"->Res["xs"],"Jacobian"->Res["jacobian"]|>
]
]


(* ::Input::Evaluate:: *)
RFM$Evaluate[hook_,Momenta_,OptionsPattern[{DEBUG->False,f128->False}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_evaluate","Arguments"->{OptionValue[f128],Momenta}|>];
(* Return *)
If[Not[MemberQ[Keys[Res],"res"]],Res,
Res["res"]
]
]


(* ::Input::EvaluateCut:: *)
RFM$EvaluateCut[hook_,CutID_,scalingFactor_,Momenta_,OptionsPattern[{DEBUG->False,f128->False,Deformation->None,DeformationJacobian->None}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_evaluate_cut","Arguments"->{
    OptionValue[f128],CutID,scalingFactor,Momenta,
    If[Not[OptionValue[Deformation]===None],OptionValue[Deformation],-1],
    If[Not[OptionValue[DeformationJacobian]===None],{Re[OptionValue[DeformationJacobian]],Im[OptionValue[DeformationJacobian]]},-1]}|>];
If[Not[MemberQ[Keys[Res],"res"]],Res,
(* Return *)
Res["res"]
]
]


(* ::Input::EvaluateIntegrand:: *)
RFM$EvaluateIntegrand[hook_,Xs_,OptionsPattern[{DEBUG->False,f128->False}]]:=Module[
{Res},

Res=ExternalEvaluate[hook,<|"Command"->"API_evaluate_integrand","Arguments"->{OptionValue[f128],Xs}|>];
If[Not[MemberQ[Keys[Res],"res"]],Res,
(* Return *)
Res["res"]
]
]
