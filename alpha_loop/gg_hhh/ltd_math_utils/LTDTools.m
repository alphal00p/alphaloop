(* ::Package:: *)

Print["The documented functions in this package are: \n ?constructCuts \n ?createSuperGraph \n ?exportAmplitude \n ?extractTensCoeff \n ?findIsomorphicGraphs \n ?getCutStructure \n ?getLoopLines \n ?getSymCoeff \n ?getSymCoeffSP \n ?importGraphs \n ?plotGraph \n ?processNumerator \n ?writeMinimalJSON \n ?writeLTDSqrtJSON \n ?translateToFeynCalc
 ----------------------------------------- 
 Needs the package FeynCalc which can installed with Import[\"https://raw.githubusercontent.com/FeynCalc/feyncalc/master/install.m\"]; InstallFeynCalc[]
 Needs the package IGraphM which can be downloaded from https://github.com/szhorvat/IGraphM. !!! \n Run: Get[\"https://raw.githubusercontent.com/szhorvat/IGraphM/master/IGInstaller.m\"] for installation. "
];
Needs["IGraphM`"];
Needs["FeynCalc`"];
constructCuts::usage="Finds all cutkovsky cuts of a graph.
	Input: 
		graph (see output from qgraf \"orientedGraph.sty\")
	Optional Input:
		NumberFinalStates (default:All): number of cuts
		DisplayCuts (default: True): plots graphs with high-lighted cuts
		DetectSelfEnergy (default: True): append entry to graph, which contains the information on the self energy insertions
	Output: all possible cut graphs
"
findIsomorphicGraphs::usage="Finds isomorphisms between graphs.
	Input: 
		List of graphs (see output from qgraf \"orientedGraph.sty\")
	Optional Input:
		edgeTagging (default: {}): List of Keyword(s) of association(s) in the graphs on which to tag by (e.g. \"partilcleType\")
	Output: List where every entry has the form
	{
		{Position of graph to be mapped, Positon of graph it maps to},
		{
		{isomorphism, permutation of edges, orientation changes: +1 = no orientation change, -1: orientation change}
		}
	}
"
plotGraph::usage="plots graph
		Input: 
			graph (see output from qgraf \"orientedGraph.sty\")
		Optional Input:
			edgeLabels(default: {}): List of Keyword(s) of association(s) in the graphs on which to tag by (e.g. \"partilcleType\")
			plotSize (default: Scaled[0.5])
"
importGraphs::usage="Import QGraf graphs generated with orientedGraphs.sty
			Input: 
				\"PathToQGrafOutput\"
			Optional Input:
				sumIsoGraphs\[Rule]True/False (* Default: True *): Align momenta of isomorphic scalar topologies and sum numerators
"
getLoopLines::usage="Appends loop-lines to graph.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
			include free edges (Boolean), Default: False, If true, edges of propagators which are free of loop-momenta are appended to the loop-lines and cut-structure
		Output:
			graphs with appended association for loop-lines
"
getCutStructure::usage="Appends cut-strucuture to graph.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
				- if graphs have no attribute \"loopLines\", it will be generated
				- if graphs have attribute \"contourClosure\" ->{\"Above\",...,\"Below\",..} this particular closure will be used. Otherwise
				  default closure from above.
		Output:
			graphs with appended association for cut-structure
"
writeMinimalJSON::usage="Exports JSON-File for LTD.
		Input: 
			graphs (List or single graphs) (see: importGraphs[output from qgraf \"orientedGraph.sty\"])
				- if graphs have no attribute \"cutStructure\", it will be generated
				- if graphs have no attribute \"name\", name will be genrated as PROCESSNAME_DIAG_DIAGNUMBER
			numAssociation (type: association): association of masses and external momenta to numerical values
		Optional Input:
			Options: {processName\[Rule] \"NameOfProcess\"(Default: \"\"), exportDirectory\[Rule] \"PathToExportDir/\" (Default: \"./\"), exportDiagramwise (Default: False),writeNumerator (Default: False), additionalKeys-{}}
		Output:
			None"
extractTensCoeff::usage=" Extracts from a given numerator the loop-momentum tensor coefficients and returns them ordered by rank.
	INPUT: 
		graph (or list of graphs): see output from importGraphs[], needs attribute  \"numerator\" .
	OPTIONAL INPUT: 
		consistencyCheckLevel (Default=1): INTEGER: 2: Run inbetween checks (for debugging);
											  1: Perform Dirac-/Lorentz-algebra and compare final result vs input expression
											  0: No-checks
	OUTPUT: graph with appended entry \"analyticTensorCoeff\" (see example)
"
getSymCoeff::usage=" Gives the symmetric sum over tensor coefficients with explicitly expanded numeric Lorentz indices. The tensor coefficients will be of the form k^0*k^1 T^{01}: All effects of the metric are already included. 
	INPUT: 
		graph (or list of graphs) with attribute \"numerator\" , \"momentumMap\", \"loopMomenta\" and/or \"analyticTensorCoeff\". 
    OPTIONAL INPUT
        outputFormat (default: \"long\"): 
             - \"short\": only non-vanishing coefficients are produced with appended index information.
             - \"long\":  all coefficients are produced in canonical order. No index information is appended.
	OUTPUT: graph with appended entry \"symmetrizedExpandedTensorCoeff\"
"
processNumerator::usage="Processes the numerator of a graph.
  INPUT:
   - graph (or list of graphs)
   - absolute path to Feynamn rules
  OPTIONAL INPUT:
   {additionalRules -> {1->1} (*used for simplification*), spinChainSimplify\[Rule]True,dimensions\[Rule]4-2 eps (* or 4 *),spinSum->False (True) ,extractTensors->False, symCoefficients\[Rule]False,coeffFormat->\"long\" (* for symmetric coefficients *)}
  OUTPUT:  Graph with processed numerator
"
createSuperGraph::usage="Creates a supegraph for a list of diagrams by loop-line merging.
  INPUT:
   - graph (or list of graphs)
  OUTPUT:  <|supergraph-> supergraph, diagrams->diagrams with mapping to supergraph|>
"

exportAmplitude::usage="Writes supergraph and amplitude .json-files for rust.
  INPUT: 
	- process: either all graphs of the process OR output from create supergraph. If expression is not preprocessed, everything will be derived on the fly.
		- comments: 			
		    - Default values are used if process entries: \"sets\" (all diagrams are put into the UV- and IR- sets) and vectors \"vectors=\" do not exist.
			- Default values are used for diagram entries: \"loopSignature\", \"positions\", \"ct\" (counterterm, False), \"factor\"
	- numeric_replacements: association with numeric replacemets of all external momenta and couplings
  OPTIONAL INPUT: See Options[exportAmplitude]
"
translateToFeynCalc::usage=" Translates old input format for algebra to FeynCalc input format*)
	INPUT:
		- expression in old format
	OPTIONAL INPUT:
		- dimensions\[Rule] 4-2eps (Default) / 4 : dimension of algebra
"
getSymCoeffSP::usage="Construct list of symmetric coefficients for given graph, if the numerator has only scalar products. The tensor coefficients will be of the form k^0*k^1 T^{01}: All effects of the metric are already included.
	INPUT:
		- graph (or list of graphs) with attribute \"numerator\" , \"momentumMap\", \"loopMomenta\" and/or \"analyticTensorCoeff\". 
	OPTIONAL INPUT:
		- numericReplacement\[Rule]{1\[Rule]1}: List with replacements for vectors/constants
"
writeLTDSqrtJSON::usage="Writes json-file for rust backend for squared-ltd graphs with Cutkosky-cuts.
  	INPUT:
		- graph with cut-info as obtained from construct-cuts
		- numeric replacements: Association for numeric values for all variables in the process (including external momenta)
	  OPTIONAL INPUT:
		- process name (string, default: new_process): the output will be saved as new_process.json
"


(* translates old input format for algebra to FeynCalc Input *)
Protect[dimensions]
ClearAll[translateToFeynCalc]
Options[translateToFeynCalc]={dimensions->4-2eps};
translateToFeynCalc[expr_,opts:OptionsPattern[]]:=
	Block[{vector,gamma,dim=OptionValue[dimensions],g,delta,T,f,spinorU,spinorV,spinorUbar,spinorVbar},
		delta[sunF[a_],sunF[b_]]:=SUNFDelta[sunF[a],sunF[b]];
		delta[sunA[a_],sunA[b_]]:=SUNDelta[sunA[a],sunA[b]];
		T[sunA[a_],sunF[i_],sunF[j_]]:=SUNTF[{sunA[a]},sunF[i],sunF[j]];
		f[sunA[a_],sunA[b_],sunA[c_]]:=SUNF[sunA[a],sunA[b],sunA[c]];
		delta[a_,b_]:=DiracIndexDelta[DiracIndex[a],DiracIndex[b]]/;Head@a=!=sunF&&Head@a=!=sunA;
		If[dim===(4-2eps),
			vector[mom_,ind_]:=FVD[mom,ind];
			gamma[s1_,k__,s2_]:=DCHN[GAD[k],s1,s2];
			g[a_,b_]:=MTD[a,b];
			spinorU[{mom_,mass_},ind_]:=DCHN[ind,SpinorUD[mom,mass]];
			spinorV[{mom_,mass_},ind_]:=DCHN[ind,SpinorVD[mom,mass]];
			spinorUbar[{mom_,mass_},ind_]:=DCHN[SpinorUBarD[mom,mass],ind];
			spinorVbar[{mom_,mass_},ind_]:=DCHN[SpinorVBarD[mom,mass],ind];
		,
			vector[mom_,ind_]:=FV[mom,ind];
			gamma[s1_,k__,s2_]:=DCHN[GA[k],s1,s2];
			g[a_,b_]:=MT[a,b];
			spinorU[{mom_,mass_},ind_]:=DCHN[ind,SpinorU[mom,mass]];
			spinorV[{mom_,mass_},ind_]:=DCHN[ind,SpinorV[mom,mass]];
			spinorUbar[{mom_,mass_},ind_]:=DCHN[SpinorUBar[mom,mass],ind];
			spinorVbar[{mom_,mass_},ind_]:=DCHN[SpinorVBar[mom,mass],ind];
		];
		expr //FCI
]


Options[importGraphs]={sumIsoGraphs->True};
importGraphs[file_:String,opts:OptionsPattern[]]:=Block[{graphs, loopMom,isos,giso,kNew,relabelLoopMom,graphToBeMapped,graphToMapTo,momToMapTo,momToBeMapped,solMom},

graphs=ToExpression["{"<>StringDrop[StringReplace[Import[file,"Text"],{"**NewDiagram**"->",","\n"->""}]<>"}",1]];
loopMom=Complement[Variables@(Union@Flatten@graphs[[All,"momentumMap"]]/.f_[x_]:>x),Union@(Cases[graphs[[All,"momentumMap"]],f_[x_]:>x,Infinity] )]
;
graphs=Append[#,"loopMomenta"->loopMom]&/@graphs;
(* only relevant for forward scattering *)
graphs=DeleteCases[graphs,x_/;!FreeQ[x[["momentumMap"]],0],{1}];
If[OptionValue[sumIsoGraphs]==True,
isos=findIsomorphicGraphs[graphs,edgeTagging->{"particleType"}];
giso=Length@isos;
(* filter isomorphisms *)

Do[
isos=Join[isos[[1;;countIso]],DeleteCases[isos[[countIso+1;;]],x_/;!FreeQ[x[[1]],isos[[countIso,1,1]] ],{1}]];
If[countIso==Length@isos-1,
Break[]
]
,{countIso,giso}];
(* build minimal graphs *)
relabelLoopMom=Array[kNew,Length@loopMom];
Do[

graphToBeMapped=graphs[[iso[[1,1]]]];
graphToMapTo=graphs[[iso[[1,2]]]];
momToMapTo=graphToMapTo[["momentumMap"]];
momToBeMapped=Permute[graphToBeMapped[["momentumMap"]],First[Sort[iso[[2,All,"permutation"]]]]]/. Thread[Rule[loopMom,relabelLoopMom]];
(* we consider the isomorphism as a pure shift in the loop-momenta which mappes the scalar subtopologies *)
solMom=(Solve[(momToMapTo^2-momToBeMapped^2)==0,relabelLoopMom]);
If[Length@solMom==0,
Print["Error in alignment of momenta. Please run with: ``sumIsoGraphs\[Rule]False``"];
Abort[],
solMom=solMom[[1]]
];
graphs[[iso[[1,2]],"numerator"]]=graphs[[iso[[1,2]],"numerator"]]+(graphToBeMapped[["numerator"]]//. Thread[Rule[loopMom,relabelLoopMom]]//. solMom);
graphs[[iso[[1,1]]]]=Append[graphs[[iso[[1,1]]]],"delete"->{True,iso[[1]]}];
,{iso,isos}];
(* delete mapped graphs *)
graphs=DeleteCases[graphs,x_/;KeyExistsQ[x,"delete"],{1}];
(* last check for aligned scalar topologies: check that every summand in every numerator has the same scalar props.... maybe a bit overcautious? *)
Do[
If[Length@((Sort/@(Cases[#,_prop,Infinity]&/@(Flatten@Apply[List,(Flatten@{gg[["numerator"]]}),{1}])/. prop[f_[x_,y_]]:>x^2//FullSimplify))//Union)>1,
Print["Error in alignment of momenta. Please run with: ``sumIsoGraphs\[Rule]False``"];
Abort[]
]
,{gg,graphs}]
];
graphs
]



Protect[NumberFinalStates,DisplayCuts,DetectSelfEnergy];
Options[constructCuts] ={ NumberFinalStates->All,DisplayCuts->True,DetectSelfEnergy-> True};
constructCuts[graph_,opts:OptionsPattern[]]:=Block[{finalRes,finalStateNum,detectSE,displayCuts,particles,momenta, trees,edges,cutEdges,forests,external,vertices,tmpDiag,cutTag,loopNum,cutGraphs,seInfo,props,cutProps,cutPos,prop,seEdges,probEdge,selfEnergy,subGEdges,sePostion,sePos,findAllSpanningTrees,externalMom},
SetAttributes[UndirectedEdge,Orderless];
(* explorative function to construct all spanning trees *)
findAllSpanningTrees[edges_,external_,loopNumber_]:=
	Block[{internalEdges,spanningTree,partitions,openEdges,addEdge},
	(* add edge connected to tree *)
		addEdge[openEdge_,tree_]:=Block[{shortEdges,addTree},
				shortEdges=List/@Cases[openEdge,UndirectedEdge[a_,b_]/;(MemberQ[VertexList@(Graph@tree),a]&&!MemberQ[VertexList@(Graph@tree),b])||(!MemberQ[VertexList@(Graph@tree),a]&&MemberQ[VertexList@(Graph@tree),b]),Infinity];
				If[Length@shortEdges>0,
					addTree=Join[tree,#]&/@shortEdges,
					addTree={tree}]
				];
		internalEdges=Complement[edges,external];
		spanningTree={{internalEdges[[1]]}};
		(* add edges until all trees are found *)
		Do[
			openEdges=Complement[internalEdges,#]&/@spanningTree;
			spanningTree=Flatten[Table[addEdge[openEdges[[treeCount]],spanningTree[[treeCount]]],{treeCount,Length@openEdges}],1];
		,{kk,Length@internalEdges-loopNumber}];
		Join[external,#]&/@((Union@(Sort/@spanningTree)))
	];

(* detects selfenery insertions *)
selfEnergy[edgesSubG_,problematicEdge_,extEdges_]:=
	Block[{graphs,result},
		graphs=Complement[edgesSubG,problematicEdge];(* two-component graph*)
		graphs=EdgeList/@(WeaklyConnectedGraphComponents@graphs);
		(* the component not having externals is the SE *)
		If[ContainsAny[graphs[[1]],extEdges],
			result=graphs[[2]],
			result=graphs[[1]]
		];
		result
	];

(* initialize data *)
If[!MemberQ[Keys@graph,"momentumMap"],
	Print[" \"momentumMap\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"edges"],
	Print[" \"edges\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"particleType"],
	Print[" \"particleType\" missing"];
	Abort[]
];
If[!MemberQ[Keys@graph,"loopMomenta"],
	Print[" \"loopMomenta\" missing. Import with function importGraphs[]"];
	Abort[]
];

finalStateNum=OptionValue[NumberFinalStates];
displayCuts=OptionValue[DisplayCuts];
detectSE=OptionValue[DetectSelfEnergy];

edges=graph[["edges"]]/. Rule[a_,b_]:>UndirectedEdge[a,b] /. DirectedEdge[a_,b_]:>UndirectedEdge[a,b]/. TwoWayRule[a_,b_]:>UndirectedEdge[a,b];
momenta=graph[["momentumMap"]];
particles=graph[["particleType"]];

externalMom=Cases[momenta,x_/;FreeQ[x,Alternatives@@(Join[graph[["loopMomenta"]],graph[["loopMomenta"]]])]]//Union;

external=Extract[edges,Position[graph[["momentumMap"]],Alternatives@@externalMom]];

vertices=VertexList@edges;
loopNum=Length@(graph[["loopMomenta"]]);
If[Length@(WeaklyConnectedGraphComponents@Complement[edges,external])>1,
Return[{Append[Append[graph,"cutInfo"->"factorizedLoopDiagram"],"selfEnergieInfo"->"factorizedLoopDiagram"]}];
Break[]
];

trees=findAllSpanningTrees[edges,external,loopNum];
(* deleting edge from spanning tree gives forest *)
forests=DeleteDuplicates@(Sort/@(Join[external,#]&/@(Flatten[Subsets[Complement[#,external],{Length@#-1-Length@external}]&/@(trees),1])));

(* a two-forest must have all vertices of the original graph: if it has not, the cut is improper \[Rule] delete all of them *)
forests=DeleteCases[forests,x_/;Length@(VertexList@x)!=Length@vertices,{1}];

(* all edges not in the 2-forest are cut *)
cutEdges={Flatten@(EdgeList/@WeaklyConnectedGraphComponents[#,in[x_]]),Flatten@(EdgeList/@WeaklyConnectedGraphComponents[#,out[x_]]), Complement[edges,#]}&/@forests;

(* Filtering *)
(* 1.) A component having ''in'' as well as ''out'' states is not allowed! *)
cutEdges=DeleteCases[cutEdges,x_/;!FreeQ[x[[1]],out[a_]],{1}];
cutEdges=DeleteCases[cutEdges,x_/;!FreeQ[x[[2]],in[a_]],{1}];

(* 2.) A cut edge is not allowed to be complete in the right or left diagram. If it is, uncut it *)
cutEdges=
	Table[
		tmpDiag=cutDiag;
		Do[
			If[(ContainsAll[VertexList@cutDiag[[1]],VertexList@({edge})]),
			tmpDiag[[1]]=Join[tmpDiag[[1]],{edge}];
			tmpDiag[[-1]]=Complement[tmpDiag[[-1]],{edge}];
		];
		If[(ContainsAll[VertexList@cutDiag[[2]],VertexList@({edge})]),
			tmpDiag[[2]]=Join[tmpDiag[[2]],{edge}];
			tmpDiag[[-1]]=Complement[tmpDiag[[-1]],{edge}];
		];
		,{edge,cutDiag[[-1]]}];
		Sort/@tmpDiag
	,{cutDiag,cutEdges}];
cutEdges=DeleteDuplicates@cutEdges;

(* Add assuciation cut-graph *)
cutTag=Table[edges /. Thread@Rule[cutEdges[[cutCount,1]],"left"]/.Thread@Rule[cutEdges[[cutCount,2]],"right"] /. Thread@Rule[cutEdges[[cutCount,3]],"cut"],{cutCount,Length@cutEdges}];

If[NumberQ[finalStateNum],
	cutTag=DeleteCases[cutTag,x_/;Count[x,"cut"]!=finalStateNum,{1}]
];

If[displayCuts==True,
	Do[
		Print@HighlightGraph[edges,Extract[edges,Position[cutTag[[cc]],"cut",Heads->False]]]
	,{cc,Length@cutTag}]
];
cutGraphs=Table[Append[graph,"cutInfo"->cutTag[[ct]]],{ct,Length@cutTag}];

If[detectSE==True,
	seInfo=Table[
		cutPos=Position[gg[["cutInfo"]],"cut"];
		props=prop[#^2//Expand //FullSimplify]&/@gg[["momentumMap"]];
		cutProps=Extract[props,cutPos];
		(* info on edges which will blow up *)
		probEdge=(Extract[#,Complement[Position[props,Alternatives@@cutProps],cutPos]]&/@{gg[["edges"]],gg[["particleType"]],gg[["momentumMap"]],gg[["cutInfo"]]});
		(* get the momenta of the selfenergy insertions *)
		seEdges=Table[
			subGEdges=Extract[gg[["edges"]],Position[gg[["cutInfo"]],probEdge[[-1,countInsertions]]]];
			selfEnergy[subGEdges,{probEdge[[1,countInsertions]]},external]
		,{countInsertions,Length@probEdge[[-1]]}];
		If[Length@(Flatten@probEdge)==0,
			probEdge={}
		];
		<|"problematicEdges"->probEdge ,
		"seDiagrams"->Table[
			sePos=Position[edges,Alternatives@@se];
			<|"EdgesSE"->Extract[graph[["edges"]],sePos],
			"particleTypeSE"->Extract[gg[["particleType"]],sePos],
			"momentumMapSE"->Extract[gg[["momentumMap"]],sePos],
			"cutInfoSE"-> Extract[gg[["cutInfo"]],sePos]
		|>
	,{se,seEdges}]|>
	,{gg,(cutGraphs/. Rule[a_,b_]:>UndirectedEdge[a,b] /. DirectedEdge[a_,b_]:>UndirectedEdge[a,b]/. TwoWayRule[a_,b_]:>UndirectedEdge[a,b])}];
	];
	finalRes=Table[Evaluate/@<|"diagram_splitting"->cutGraphs[[infoCount,"cutInfo"]],"self_energy"->seInfo[[infoCount]]|>,{infoCount,Length@cutGraphs}];
	finalRes=Append[graph,Evaluate/@<|"cutDiagInfo"->finalRes|>];
	finalRes=constructLeftRightDiagrams[finalRes]
]


(* represent graph by: {{edge1: (Rule[a,b] for directed edge, TwoWayRule[a,b] for undirected edge),  ,tag1:List},...,{edgeN,tagN:List}} *)
Protect[edgeTagging];
Options[findIsomorphicGraphs] ={ edgeTagging->{}};
findIsomorphicGraphs[graphs_:List,tag:OptionsPattern[]]:=Block[{edges,myTaggedGraphs,tagKeys,tagging,taggingToNum,$uniqueTag, graphsEdgeColored,graphsEdgeColoredSimple,isIsomorph,$myAssoc,isomorphisms,shifts,g1,g2,shiftsAndPermutations,orientation,mappedGraph,graphID1,graphID2,permutation,mapping},
	SetAttributes[TwoWayRule,Orderless];
	edges=graphs[[All,"edges"]];
	If[!FreeQ[edges,_Missing],
		Print["graph is missing Key :\"edges\""];
		Abort[]		
		];
	tagKeys = OptionValue[edgeTagging];

	(* taggings become colors on edge-colored graphs, multiple taggings will be translated to single prime which becomes the color of the edge *)
	If[Length@tagKeys!=0,
		tagging=Values/@graphs[[All,tagKeys]];
		If[Length@(Union@(Length/@tagging))!=1||Length@(Union@Map[Length,tagging,{2}])!=1,
			Print["Not every edge has a tagging for choosen edgeTagging: "<>ToString[tagKeys]];
			Abort[]		
		];
		If[!FreeQ[tagging,_Missing],
			Print["Not every edge has a tagging for choosen edgeTagging: "<>ToString[tagKeys]];
			Abort[]		
		];
		tagging=Transpose/@tagging,		
		(* if there is not tagging, introduce random tagging *)
		tagging=(edges/. Rule[a_,b_]:>{$uniqueTag} /. TwoWayRule[a_,b_]:> {$uniqueTag}/. DirectedEdge[a_,b_]:> {$uniqueTag}  /. UndirectedEdge[a_,b_]:> {$uniqueTag});
		];
	
(* start prime range late, so that translation of multi-edge graphs is still unique *)
	taggingToNum=Thread@(Rule[Variables@tagging,Prime[Range[100,Length@(Variables@tagging)+99]]]);

(* translate to undirected graphs and apply tagging *)
	myTaggedGraphs=MapThread[Riffle[{#1},{#2}]&,{edges,tagging},2];

	graphsEdgeColored=myTaggedGraphs /. taggingToNum /. {DirectedEdge[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x} /. {UndirectedEdge[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x}/. {Rule[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x} /. {TwoWayRule[a_,b_],x_List}:>{TwoWayRule[a,b],Times@@x};

(* get isomorphism for graphs which are isomorphic; IGVF2 needed for multi-edge graphs: 
	There is no algorithm working for multi-edge graphs directly. Therefore translate into edge-colored simple graphs! The color is the multiplicity of the edge e.g.: g1Multiedge={{1\[TwoWayRule]4,5},1\[TwoWayRule]4,5}} --> g1ColouredSimple={1\[TwoWayRule]4\[Rule]2*5}*)
	graphsEdgeColoredSimple=(Counts/@(Sort/@graphsEdgeColored)) /.Association:>$myAssoc/.Rule[List[TwoWayRule[a_,b_],color1_],color2_]:> Rule[TwoWayRule[a,b],color1*color2] /. TwoWayRule->TwoWayRule /. $myAssoc:>Association;

(* thats a bit of laziness, since over computation of half the entries.... *)
	isIsomorph=Table[IGVF2IsomorphicQ[{Graph[VertexList[Keys@graph1],Keys[graph1]],"EdgeColors"->graph1},{Graph[VertexList[Keys@graph2],Keys[graph2]],"EdgeColors"->graph2}],{graph1,graphsEdgeColoredSimple},{graph2,graphsEdgeColoredSimple}]//LowerTriangularize;
	isomorphisms=(Position[isIsomorph,True]/. {i_,i_}-> Nothing);

(* compute isomorphism for isomorphic graphs *)
	shifts=ConstantArray[0,Length[isomorphisms]];
	Do[
		g1=graphsEdgeColoredSimple[[isomorphisms[[i,1]]]];
		g2=graphsEdgeColoredSimple[[isomorphisms[[i,2]]]];
		shifts[[i]]={{isomorphisms[[i,1]],isomorphisms[[i,2]]},IGVF2FindIsomorphisms[{Graph[VertexList[Keys@g1],Keys[g1]],"EdgeColors"->g1},{Graph[VertexList[Keys@g2],Keys[g2]],"EdgeColors"->g2}]};
	,{i,1,Length[isomorphisms]}];
(* find edge permutations *)
	shiftsAndPermutations=Table[
		Join[{shifts[[countShifts,1]]},
		{Table[
		
			{shifts[[countShifts,2,isoCount]],
			If[Length@(Complement[graphsEdgeColored[[shifts[[countShifts,1,1]]]]/.shifts[[countShifts,2,isoCount]],graphsEdgeColored[[shifts[[countShifts,1,2]]]]])==0,
			FindPermutation[graphsEdgeColored[[shifts[[countShifts,1,1]]]]/.shifts[[countShifts,2,isoCount]],graphsEdgeColored[[shifts[[countShifts,1,2]]]]],
			"errorPerm"
			]}
		,{isoCount,Length@shifts[[countShifts,2]]}]}]
	,{countShifts,Length@shifts}];
	
	shiftsAndPermutations=shiftsAndPermutations /. {x_,"errorPerm"}:>Nothing;
	
(* find if edge direction changed: No=1,Yes=-1*)
	Do[
		graphID1=shiftsAndPermutations[[shiftCount,1,1]];
		graphID2=shiftsAndPermutations[[shiftCount,1,2]];
		Do[
			mapping=shiftsAndPermutations[[shiftCount,2,isoCount,1]];
			permutation=shiftsAndPermutations[[shiftCount,2,isoCount,-1]];
			mappedGraph=Permute[myTaggedGraphs[[graphID1]],permutation]/.mapping /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b});
			orientation=(Keys[mappedGraph[[All,1]]]-Values[mappedGraph[[All,1]]])/(Keys[(myTaggedGraphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))]-Values[(myTaggedGraphs[[graphID2,All,1]] /. TwoWayRule[a_,b_]:>Rule@@(Sort@{a,b}))])//Simplify;
			shiftsAndPermutations[[shiftCount,2,isoCount]]=Join[shiftsAndPermutations[[shiftCount,2,isoCount]],{orientation}];
		,{isoCount,Length@shiftsAndPermutations[[shiftCount,2]]}];
	,{shiftCount,Length@shiftsAndPermutations}];
	shiftsAndPermutations=Table[{shiftsAndPermutations[[countIso,1]],Association/@Table[Thread@Rule[{"vertexMap","permutation","orientationChange"},iso],{iso,shiftsAndPermutations[[countIso,2]]}] },{countIso,Length@shiftsAndPermutations}]
];


(* GraphComputation`GraphPlotLegacy is a hack, since mathematica updated GraphPlot *)
Protect[edgeLabels,plotSize]
Options[plotGraph] ={ edgeLabels->{},plotSize->Scaled[0.5]};
plotGraph[graph_,opt:OptionsPattern[]]:=Block[{mygraph,imSize=OptionValue[plotSize]},
	Do[
	mygraph=Transpose@(Values@(gg[[Join[{"edges"},Flatten@{OptionValue[edgeLabels]}]]]))/. TwoWayRule->Rule /. UndirectedEdge->Rule /. DirectedEdge->Rule /. {Rule[a_,b_],c__}:>{Rule[a,b],Flatten@{c}} /. {Rule[a_,b_]}:>{Rule[a,b],{Rule[a,b]}};
	Print[
	If[TrueQ[$VersionNumber >= 12.0], GraphComputation`GraphPlotLegacy, GraphPlot][mygraph,EdgeRenderingFunction->({{RGBColor[0.22,0.34,0.63],Arrowheads[0.015],Arrow[#]},If[#3=!=None,Text[ToString[#3],Mean[#1],Background->White],{}]}&),MultiedgeStyle->.3,SelfLoopStyle->All,ImageSize->imSize]];
	,{gg,Flatten@{graph}}]
];



getLoopLines[graphs_,includeFreeEdges_:False]:=Block[{eMom,lMom,lLines,lProp,orientation,res,mass,pos,eProp},
res=Table[
If[!KeyExistsQ[graph,"loopMomenta"],
	Print["Error: Graphs have no entry \"loopMomenta\". Import with importGraph[]"];
	Abort[]
];
eMom=Complement[Variables@(Union@Flatten@graph[["momentumMap"]]),graph[["loopMomenta"]]];
lLines=CoefficientArrays[((graph[["momentumMap"]] /. Thread[Rule[eMom,0]]/. 0:>Nothing)^2 //FullSimplify)/.(x_)^2:>x// DeleteDuplicates,graph[["loopMomenta"]]][[2]]//Normal;
lProp=Table[
	orientation=(ReplaceAll[#,x_/;!NumberQ[x]:>0]&/@(FullSimplify@((graph[["momentumMap"]]/. Thread[Rule[eMom,0]])/(Dot[graph[["loopMomenta"]],ll])) ));
	pos=Position[Abs@orientation,1,Heads->False];
	<|"signature"->ll,
	"propagators"->Table[
	<|"mass"->mass[graph[["particleType",pp]]],
	  "edge_name"-> If[KeyExistsQ[graph,"edge_names"],
	  graph[["edge_names",pp]],
	  "edge"<>ToString@(pp)
	  ],
	  "power"->If[KeyExistsQ[graph,"powers"],
	  graph[["powers",pp]],
	  1
	  ],
	  "q"->(graph[["momentumMap",pp]]*orientation[[pp]] /. Thread[Rule[graph[["loopMomenta"]],0]] /. 0->ConstantArray[0,4])
	  |>
	  ,{pp,Flatten@pos}]	  
	|>	 
,{ll,lLines}];
graph=Append[graph,Association@{"loopLines"->lProp}];
If[includeFreeEdges===True,
	graph=getCutStructure[graph];
	eProp=Association@{Association@("signature"->ConstantArray[0,Length@(graph[["loopMomenta"]])]),
	Association@("propagators"->
	Table[
	pos=Position[graph[["momentumMap"]],freeEdge,Heads->False];
	Association@{
	"mass"->mass[Extract[graph[["particleType"]],pos][[1]]]
	,
	 "edge_name"->"edge"<>ToString@((Flatten@(pos))[[1]])
	 ,
	 "power"->If[KeyExistsQ[graph,"powers"],
	  Extract[graph[["powers"]],pos][[1]],
	  1
	  ],
	  "q"->Extract[graph[["momentumMap"]],pos][[1]]	  
	}
	,{freeEdge,Cases[graph[["momentumMap"]],x_/;FreeQ[x,Alternatives@@(Flatten@{graph[["loopMomenta"]],in,out})],1]}])};
	If[Length@(Values@eProp[["propagators"]])>0,
	graph[["loopLines"]]=Append[graph[["loopLines"]],eProp];
	graph[["cutStructure"]]=Append[#,0]&/@graph[["cutStructure"]];
	graph
	,
	graph
	]
	,
	graph
]

,{graph,Flatten@{graphs}}];
If[Length@res==1,
res=res[[1]]
];
res
]




$fileDirectory=If[$Input==="",FileNameJoin[{NotebookDirectory[],""}],FileNameJoin[{DirectoryName@$InputFileName,""}]];
getCutStructure[graphs_]:=Block[{session,result,pyTest,myGraphs,ctStruc,ll},

If[!FileExistsQ[$fileDirectory<>"/compute_cut_structure.py"],
	Print["No cut-structure script found in: \""<>$fileDirectory<>"/compute_cut_structure.py"<>"\""];
	Abort[]
];
session=StartExternalSession[<|"System"->"Python","Version"->"3","ReturnType"->"String"|>];
If[!StringMatchQ[session["Version"],"3"~~___],
	Print["Error: No external python3 evaluator registered."];
	Print["see: \" https://reference.wolfram.com/language/workflow/ConfigurePythonForExternalEvaluate.html \""];
	Abort[]
];
pyTest=ExternalEvaluate[session,"import sys"];

If[pyTest=!="Null",
	Print["Error in ExternalEvaluate[session,\"import sys\"]"];
	Print["possible cause: Usage of python 3.8"];
	Print["possible fix: \n Change \" yourPathTo/WolframClientForPython/wolframclient/utils/externalevaluate.py \""];
	Print["66c66
		<         exec(compile(ast.Module(expressions, []), '', 'exec'), current)
		---
		>         exec(compile(ast.Module(expressions), '', 'exec'), current)
	"];
	Abort[]
];
ExternalEvaluate[session,"sys.path.append('"<>$fileDirectory<>"')"];
ExternalEvaluate[session,"import compute_cut_structure as cct"];

If[ContainsAny[KeyExistsQ[#,"loopLines"]&/@(Flatten@{graphs}),{False}],
	myGraphs=getLoopLines[graphs],
	myGraphs=graphs
];
result=Table[
	ll=gg[["loopLines"]][[All,"signature"]];
	ExternalEvaluate[session,"loop_line_signatures = "<>ExportString[ll,"PythonExpression"]];
	ExternalEvaluate[session,"cut_structure_generator = cct.CutStructureGenerator(loop_line_signatures)"];
	If[KeyExistsQ[gg,"contourClosure"],
		ExternalEvaluate[session,"contour_closure = "<>StringReplace[ExportString[(gg[["contourClosure"]]/. "Above"->"cct.CLOSE_ABOVE"/. "Below"->"cct.CLOSE_BELOW"),"PythonExpression"],"\""->"" ]],
		ExternalEvaluate[session,"contour_closure = "<>StringReplace[ExportString[ConstantArray["cct.CLOSE_ABOVE",Length@ll[[1]]],"PythonExpression"],"\""->"" ]];
	];
	ctStruc=ImportString[ExternalEvaluate[session,"cut_structure_generator(contour_closure)"],"PythonExpression"];
	Append[gg,"cutStructure"->ctStruc]
,{gg,Flatten@{myGraphs}}];
Clear[session];
If[Length@result==1,
	result=result[[1]]
];
result
]



Protect[processName,exportDirectory,exportDiagramwise,writeNumerator,additionalKeys,export];
Options[writeMinimalJSON]={processName->"",exportDirectory->"./",exportDiagramwise->False,writeNumerator->False,additionalKeys->{},export->True}
writeMinimalJSON[graphs_,numAssociation_Association,opts:OptionsPattern[]]:=Block[{mygraphs=Flatten@{graphs},inMom,outMom,extAsso,cutAsso,llAsso,lnAsso,nameAsso,diagCount=0,procName,exportDir,pLong,fullAsso,diagramwise,processAsso,expAsso,evaluateNumerator,numeratorAsso},
(* small wrapper for plugging in numerical values in the numerator *)

(* small wrapper for numerical evaluation of numerators *)	 	 
evaluateNumerator[numerator_,numRepl_]:=Block[{numericNum=numerator,Pair},
  numericNum=numericNum//FCI;
  Pair[Momentum[x_List,d___],Momentum[y_List,d___]]:=x[[1]]*y[[1]]-x[[2;;]].y[[2;;]];
  (* vectors are assumed to be contravariant *)
   numericNum=(numerator//. numRepl //.x_List[y_Integer]:>x[[y+1]]);
  If[Length@(Variables@numericNum)!=0,
  Print["Error: The numerator coefficients: "<>ToString[Variables@numericNum,InputForm]<>" have no numeric value!"];
  Abort[];  
  ];
  (* short vs long export format *)
  If[!ListQ[numericNum[[1]]],
    numericNum=((numericNum)/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]),
    numericNum[[All,2]]=(numericNum[[All,2]]/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]);    
  ];
  numericNum
];


diagramwise=OptionValue[exportDiagramwise];
procName=OptionValue[processName];
exportDir=OptionValue[exportDirectory];
If[ContainsAny[KeyExistsQ[#,"loopLines"]&/@(Flatten@{graphs}),{False}]||ContainsAny[KeyExistsQ[#,"cutStructure"]&/@(Flatten@{graphs}),{False}],
	mygraphs=Flatten@{getCutStructure[graphs]},
	mygraphs=Flatten@{graphs}
];
(* loop over graphs *)
processAsso=Table[
	diagCount+=1;
	(* check input *)
	If[KeyExistsQ[mygraph,"momentumMap"],
		inMom=Cases[mygraph[["momentumMap"]],_in,Infinity]//Union;
		outMom=Cases[mygraph[["momentumMap"]],_out,Infinity]//Union
		,
		If[KeyExistsQ[mygraph,"externalMomenta"],
			inMom=Cases[mygraph[["externalMomenta"]],_in,Infinity]//Union;
			outMom=Cases[mygraph[["externalMomenta"]],_out,Infinity]//Union
			,
			Print["Can not determine external momenta."];
			Abort[]	
		]	
	];
	(* check if stuff is numeric *)
	If[MemberQ[NumericQ/@(Flatten@(inMom /. in[x_]:>x /. numAssociation)),False],
		Print["Error: association is missing numeric value for incoming momenta"];
		Abort[];
	];
	If[MemberQ[NumericQ/@(Flatten@(outMom/. out[x_]:>-x /. numAssociation)),False],
		Print["Error: association is missing numeric value for outgoing momenta"];
		Abort[];
	];
(* check momentum conservation *)
	If[Abs@(Total@((Total@inMom-Total@outMom) /. f_[x_]:>x /. numAssociation))/Abs@(Total@((Total@inMom+Total@outMom+10^-16) /. f_[x_]:>x /. numAssociation))>10^-10,
		Print["Error: rel. error of momentum conservation is worse than 10^(-10)"];
		Abort[]
	];
(* check mass replacements *)
	If[MemberQ[NumericQ/@(Flatten@(mygraph[["loopLines",All,"propagators",All,"mass"]]/.numAssociation)),False],
		Print["No numerical assignment for :"<>ToString[(Flatten@(mygraph[["loopLines",All,"propagators",All,"mass"]]/.numAssociation))/;x_/;NumericQ[x]:>Nothing//Union]];
	];
(* replacement for yaml *)
	extAsso=<|"external_kinematics"->Join[inMom,outMom]/. in[x_]:>x /. out[x_]:>-x/.numAssociation|>;
	cutAsso=<|"ltd_cut_structure"->mygraph[["cutStructure"]]|>;
	(* all the BS is because mathematica exports 0.e-31, which is not readable by yaml *)
	llAsso=
		Map[Evaluate,<|"loop_lines"->Table[
			<|<|"end_node"->69|>,
			<|"propagators"->KeyMap[Replace[{"mass"->"m_squared","edge_name"->"name"}]]/@(l[["propagators"]]/.mass[x_]:>mass[x]^2/.numAssociation /. {x_,y__}/;NumericQ[x]:>ImportString[ExportString[{x,y},"Real64"],"Real64"])|>,
			l[[{"signature"}]],
			<|"start_node"->70|>
			|>,{l,(mygraph[["loopLines"]])}]|>,5];
	lnAsso=Map[Evaluate,<|"n_loops"->Length@mygraph[["loopLines",1,"signature"]]|>];
	If[KeyExistsQ[mygraph,"name"],
		nameAsso=<|"name"->mygraph[["name"]]|>,
		nameAsso=<|"name"->procName<>"diag_"<>ToString[diagCount]|>
	];
	If[KeyExistsQ[mygraph,"maximum_ratio_expansion_threshold"],
		expAsso=<|"maximum_ratio_expansion_threshold"->mygraph[["maximum_ratio_expansion_threshold"]]|>,
		expAsso=<|"maximum_ratio_expansion_threshold"->-1|>
	];
	fullAsso=extAsso;
	If[OptionValue[writeNumerator],
	    If[!KeyExistsQ[mygraph,"symmetrizedExpandedTensorCoeff"],
	         mygraph=getSymCoeff[mygraph]
	      ];
		numeratorAsso=evaluateNumerator[mygraph[["symmetrizedExpandedTensorCoeff"]],numAssociation];
		If[ListQ[numeratorAsso[[-1]]],
		numeratorAsso=Association@("numerator_tensor_coefficients_sparse"->evaluateNumerator[numeratorAsso,numAssociation]),
		numeratorAsso=Association@("numerator_tensor_coefficients"->evaluateNumerator[numeratorAsso,numAssociation])		
		]		
		,
		numeratorAsso=<|"numerator_tensor_coefficients"->{{1,0}}|>
	];
	
	Do[fullAsso=Append[fullAsso,asso],{asso,{llAsso,cutAsso,expAsso,lnAsso,numeratorAsso,nameAsso}}];
	If[Length@OptionValue[additionalKeys]>0,
	  
	  If[Union@(Flatten@{KeyExistsQ[mygraph,#]&/@OptionValue[additionalKeys]})==={True},
	  	  fullAsso=Prepend[fullAsso,mygraph[[OptionValue[additionalKeys]]]],
	  	  Print["Warning: Could not find all additional keys in diagram "<>ToString@diagCount<>". They will be ommitted on export."]	  
	  ]
	];
	fullAsso
	,{mygraph,mygraphs}];
	If[OptionValue[export]===True,
	If[DirectoryQ[exportDir],
		Export[exportDir<>"allDiags"<>procName<>".json",processAsso,"JSON"],
		Print["Couldn't find exportDirectory. Export to standard location."];
		Export["./allDiags"<>procName<>".json",processAsso,"JSON"]
	];
	];
	processAsso
]


Protect[consistencyCheckLevel]
Options[extractTensCoeff]={consistencyCheckLevel->1}
extractTensCoeff[graphs_List,opts:OptionsPattern[]]:=extractTensCoeff[#,opts]&/@graphs;
extractTensCoeff[graph_,opts:OptionsPattern[]]:=Block[
	{amp,loopMom,extVect,maxRank,epsK1K1,sp,mySP,tensDecomp,indexSpace,allStruc,allStrucReplRule,allStrucReplRuleVals,res,dummyIndex=Unique[alpha],dummyIndex2=Unique[beta],cleanUpRules
	,tensIndexSet,tensIndices,tensFinal,tensTMP
	,resTMP,replRuleTensExtraction,tens,reorderTensors,repl1,repl2,repl3,tensFinalCheck,finalRes,lMomInd,rewriteRule,ampMod,consistencyCheck=OptionValue[consistencyCheckLevel],$randomDummy,dim,vector,g,
	myCheckII,$vv},
SetAttributes[sp,Orderless];
loopMom=graph[["loopMomenta"]];
amp=$randomDummy*graph[["numerator"]];
ampMod=amp//FCI//ScalarProductExpand//MomentumExpand//FCI;
(* expect dimensionality to be consistent *)
dim=Flatten@FCGetDimensions[ampMod];

If[dim==={D},
	vector[a_,b_]:=FCI@FVD[a,b];
	dummyIndex[x_]:=LorentzIndex[dummyIndex2[x],D]/;NumericQ@x;
	g[x_,y_]:=FCI@MTD[x,y];
	,
	vector[a_,b_]:=FCI@FV[a,b];
	dummyIndex[x_]:=FCI@LorentzIndex[dummyIndex2[x]]/;NumericQ@x;
	g[x_,y_]:=FCI@MT[x,y]
];

extVect=Complement[(graph[["momentumMap"]]/. out[x_]:>x /. in[x_]:>x),loopMom];
tensIndexSet=Table[ToExpression@("lInd"<>ToString@lCount),{lCount,Length@loopMom}];


If[amp==0,
	Return[Append[graph,"analyticTensorCoeff"->{{ConstantArray[0,Length@loopMom],0}}]]
];
(* determine max tensorrank *)
maxRank=Max@@(Exponent[(ampMod //. f_[x__]/;!FreeQ[{x},Alternatives@@loopMom]&&StringMatchQ[ToString@f,Alternatives["Pair","DiracChain","DiracGamma"]]:>epsK1K1^Length@(Cases[{x},Alternatives@@loopMom,Infinity]) Unique[f]),epsK1K1]);

tensDecomp=(ampMod//Expand)/. Pair[a_,b_]^pow_/;(!FreeQ[a,Alternatives@@loopMom]||!FreeQ[b,Alternatives@@loopMom]):>mySP[{a,b},pow]; (* scalarproducts can be to higher powers *)
(* find all structures *)
allStruc=Join[
	Cases[tensDecomp,f_[x__]/;!FreeQ[{x},Alternatives@@loopMom]&&StringMatchQ[ToString@f,Alternatives["Pair","DiracChain","DiracGamma"]]:>{f[x],Length@(Cases[{x},Alternatives@@loopMom,Infinity])},Infinity]//Union,
	Cases[tensDecomp, mySP[x_,pow_]:>{mySP[x,pow],Length@(Cases[{x},Alternatives@@loopMom,Infinity])*pow} ,Infinity]//Union
	];
(* generate consecutive indexing *)
(* allStruc[[1,-1]] goes to  [startindex(i),endIndex(i)] where startindex(i+1)=endindex(i)+1 *)
If[Length@allStruc>0,
  allStruc[[1,-1]]={2,allStruc[[1,-1]]+1};
Do[
 allStruc[[countStruc,-1]]={allStruc[[countStruc-1,-1,-1]]+1,allStruc[[countStruc-1,-1,-1]]+allStruc[[countStruc,-1]]}
	,{countStruc,2,Length@allStruc}];
	
	
(* define unique indices *)
allStrucReplRuleVals=Table[
	res={struc[[1]],{}};
	Do[
		res={ReplacePart[#,FirstPosition[#,Momentum[x_,y___]/;MemberQ[loopMom,x]]:>dummyIndex[index]]&@(res[[1]] ),Join[res[[2]],{vector[Extract[#,FirstPosition[#,x_/;MemberQ[loopMom,x]]]&@(res[[1]] ),dummyIndex[index]]}]}
	,{index,Range[struc[[2,1]],struc[[2,-1]]]}
	];
	Times@@(Flatten@res)
	,{struc,(allStruc/. mySP[x_,pow_]:>ConstantArray[sp@@x,pow])}];
	

cleanUpRules=Dispatch@{sp[x_,y_]/;(FreeQ[x,Alternatives@@Join[extVect,loopMom]]&&FreeQ[y,Alternatives@@Join[extVect,loopMom]]):>g[x,y],
						 sp[x_,y_]/;(!FreeQ[x,Alternatives@@Join[extVect,loopMom]]&&FreeQ[y,Alternatives@@Join[extVect,loopMom]]):>vector[x,y],
						 sp[x_,y_]/;(FreeQ[x,Alternatives@@Join[extVect,loopMom]]&&!FreeQ[y,Alternatives@@Join[extVect,loopMom]]):>vector[y,x]
						 };	
allStrucReplRuleVals=allStrucReplRuleVals //.cleanUpRules ;
allStrucReplRule=Thread@(Rule[allStruc[[All,1]],allStrucReplRuleVals])
,
allStrucReplRule=1->1
];
tensDecomp=(tensDecomp//. allStrucReplRule);

(* consistency check --> new expression matches input: Extraction of loop-momenta worked *)
If[consistencyCheck>1,
	Print["check I started"];
	If[Simplify@(ampMod-DiracSimplify@(Contract[tensDecomp]))=!=0
		,
		Print["Error in tensor extraction: Check Input"]; Abort[]
		,
		PrintTemporary["check I passed"]
	 ]
];

tensDecomp=Table[SeriesCoefficient[(tensDecomp /. vector[x_,ind_]/;MemberQ[loopMom,x]:>epsK1K1 vector[x,ind]),{epsK1K1,0,rank}],{rank,0,maxRank}];
(* consistency check --> new expression matches input: rank extraction worked *)
If[consistencyCheck>1,
	PrintTemporary["check II started"];
	myCheckII=Simplify@(ampMod-Total@(DiracSimplify@(Contract[tensDecomp])));
	If[Simplify@(myCheckII)=!=0
		,
		Print["Error in tensor extraction: Check Input"]; Abort[]
		,
		PrintTemporary["check II passed"]
	 ]
];
(* redefine mathematica default ordering in times *)
reorderTensors[myexpr_, pats_List] :=
  Module[{h, rls},
    rls = MapIndexed[x : # :> h[#2, Replace[x, rls, -1]] &, pats];
    HoldForm @@ {myexpr /. rls} //. h[_, x_] :> x
  ];

tensFinal=Table[
	resTMP=tensDecomp[[rankCount+1]];
	Do[
		tensIndices=Table[tensIndexSet[[loopCount]][ind],{ind,1,rankCount}];
		(*works since we defined unique indices*)
		tensTMP=Cases[tensDecomp[[rankCount+1]],vector[loopMom[[loopCount]],ind_]:>ind,Infinity]//Union;
		
		replRuleTensExtraction=Dispatch@{
		tens[x_,y_]tens[z_,v_]:>tens[Join[x,z],Join[y,v]]		 	
		};
		repl1=f_[x___]/;(MemberQ[{Pair,DiracChain,DiracGamma},Head@(f[x])]&&!FreeQ[Flatten@{x},Alternatives@@tensTMP]):>tens[{f[x]}, Cases[{x},Alternatives@@tensTMP,Infinity]];
		repl2=tens[{start___,tens[a_,x_],rest___} ,y_]:> tens[Flatten@Join[{start},{a},{rest}], Join[x,y]];
		repl3 = tens[x_,y_]/;!FreeQ[y,Alternatives@@tensTMP] :> ReplaceAll[tens[x, y],Thread@(Rule[Cases[y,Alternatives@@tensTMP,Infinity],Flatten@{tensIndices[[(*rankCount+1,*)1;;Length@Cases[y, Alternatives@@tensTMP,Infinity]]]} ] )];
		resTMP=(ReleaseHold@(ReplaceRepeated[reorderTensors[Expand[(# //.vector[loopMom[[loopCount]],ind_]-> 1 /. repl1//. repl2)],{tens}],replRuleTensExtraction])&@((resTMP)))/.repl3;  
	,{loopCount,Length@loopMom}];
	
	resTMP
,{rankCount,0,maxRank}];

(* order entries by loop-momenta tensors *)
tensFinal=Flatten[Table[
	Table[
		{tensCoeff,(tensFinal[[rankCount+1]] //.tens[x_,ind_]/;(Count[ind,#]&/@(Blank/@(tensIndexSet[[;;Length@loopMom]]))!=tensCoeff):>0 ) }
		,{tensCoeff,Sort[(Flatten[#,1]&@(Permutations/@(PadRight[#,Length@loopMom]&/@IntegerPartitions[rankCount,Length@loopMom])))]//Reverse}]	
,{rankCount,0,maxRank}],1];
tensFinal=tensFinal //. tens[x_,y_]:>Times@@x;

tensFinal=tensFinal //. {loopList_List,{x_}}/;(Length@loopList==Length@loopMom && FreeQ[loopList,Alternatives@@extVect]):>{loopList,x};
 tensFinal=DeleteCases[tensFinal,{x_List,0},1];
 (* consistency check --> final expression matches input *)

(* final check*)
If[consistencyCheck>=1,
	PrintTemporary["Final check started"];
	(* replace integers by loop-momenta *)
	tensFinalCheck=(tensFinal//. {indX_List,x_}/;AllTrue[indX,IntegerQ]:>{loopMom^2 loopMom^indX,x});
	tensFinalCheck=tensFinalCheck /.mom_^pow_/;MemberQ[loopMom,mom]:>(Times@@(Table[Evaluate@$vv[mom,Extract[tensIndexSet,Position[loopMom, mom]][[1]][cc]],{cc,-2,pow-2}]));
	tensFinalCheck=tensFinalCheck//.$vv[ll_,_[num_]]/;(num <= 0):>1 /.$vv:>vector;
	
	(* restore input structure *)
	tensFinalCheck=DiracSimplify[(Contract[(tensFinalCheck //. {loopList_List,x_}/;(Length@loopList==Length@loopMom && FreeQ[loopList,Alternatives@@extVect]):>(Times@@(loopList)*x))]),Expanding->True];
	
	tensFinalCheck=(Total@tensFinalCheck);
	
	tensFinalCheck=FullSimplify@(FullSimplify@((ampMod-tensFinalCheck)//ExpandScalarProduct));
	
	If[tensFinalCheck=!=0
	,
	Print["Error in final check: likely a bug :("];Abort[],
	PrintTemporary["Final check passed"]
	]
	];
	
	(* replace lIndX[Y] by lMomInd[x][y]... x: loop-momentum number, y: lorentz-index: mu[y]...  ugly but was only fix which came to mind so \.af\_(\:30c4)_/\.af  *)
	rewriteRule=Thread[RuleDelayed[Compose[#[[1]],#[[2]]]&/@(Transpose@{tensIndexSet,ConstantArray[Pattern[x,Blank[]],Length@tensIndexSet]}),Evaluate[Compose[#[[1]],#[[2]]]&/@(Transpose@{Array[lMomInd,Length@tensIndexSet],ConstantArray[x,Length@tensIndexSet]})]]];
	
	tensFinal=tensFinal/.rewriteRule /. $randomDummy->1;
	finalRes=Append[graph,"analyticTensorCoeff"->tensFinal]
]



Protect[outputFormat];

Options[getSymCoeff] ={ outputFormat->"long", consistencyCheckLevel->1 };

getSymCoeff[graphs_,opts:OptionsPattern[]]:=Block[
{tensDecompNum,resTmp,res,indSet,g,rank,replV,replG,format,indShift,replRest,fullResult,resTmpTmp,dim,vector},
    format=OptionValue[outputFormat];	
	g[int1_Integer,int2_Integer]:=If[int1!=int2,0, If[int1==0,1,-1]];
	fullResult=Table[
    If[KeyExistsQ[graph,"analyticTensorCoeff"],
      tensDecompNum=graph[["analyticTensorCoeff"]]
      ,
      tensDecompNum=extractTensCoeff[graph,consistencyCheckLevel->OptionValue[consistencyCheckLevel]];
      tensDecompNum=tensDecompNum[["analyticTensorCoeff"]]
    ];
    
    dim=Flatten@FCGetDimensions[tensDecompNum[[1]]];
    
    If[format==="long",
      rank=tensDecompNum[[All,1]];
      (* generate complete rank upto max rank *)
      rank=Flatten[Table[Sort@Flatten[Permutations/@IntegerPartitions[cR,{Length@rank[[-1]]},Range[0,cR]],1],{cR,0,Total@rank[[-1]]}],1];
       ,
      (* only non-vanishing coefficients*)
      rank=tensDecompNum[[All,1]]
    ];
    

	If[dim==={D},
		vector[a_,b_]:=FCI@FVD[a,b];
		g[x_,y_]:=FCI@MTD[x,y];
	,
		vector[a_,b_]:=FCI@FV[a,b];
		g[x_,y_]:=FCI@MT[x,y]
];
   res=
   Flatten[
   Table[
    (* translate given rank to all possible (symmetrized sets) of lorentz momenta *)

    (* thats basically a cartesian product *)
	indSet=Replace[List[x__]:>Flatten[Outer[List,x,1],Length@rank[[rankCount]]-1]]/@(Distribute[(rank[[rankCount]]/. x_/;NumericQ[x]:>(Permutations/@(Union@(Sort/@(Tuples[Range[0,3],{x}]))))),List]);
(*	    Print["-------"];
	Print[rank[[rankCount]]];
	Print[indSet];
	Print["-------"];*)
	If[MemberQ[tensDecompNum[[All,1]],rank[[rankCount]]],
		resTmp=Extract[tensDecompNum[[All,2]],Position[tensDecompNum[[All,1]],rank[[rankCount]]]];
		resTmp=resTmp[[1]]
		,
		resTmp=0	
	];
	indShift=Table[ConstantArray[(llCount-1)*4,(rank[[rankCount,llCount]])],{llCount,Length@rank[[rankCount]]}];
	(* loop over all symmetrized sets *)
	Table[
(*	    Print["new"];*)
	    resTmpTmp=
	    Total@Table[
(*	     Print[ind];
	     Pause[.2];*)
	     replV=vector[x_,lMomInd[loopMomNum_][loopIndex_]]:>If[ind[[loopMomNum,loopIndex]]==0,x[ind[[loopMomNum,loopIndex]]],-x[ind[[loopMomNum,loopIndex]]]];
		 replG=g[lMomInd[loopMomNum_][loopIndex_],lMomInd[loopMomNum2_][loopIndex2_]]:>g[ind[[loopMomNum,loopIndex]],ind[[loopMomNum2,loopIndex2]]];
		 replRest=lMomInd[loopMomNum_][loopIndex_]:>covartiantLI[ind[[loopMomNum,loopIndex]]];
	     resTmp //. replV //. replG //. replRest
	    ,{ind,myInd}];
	    {Flatten@(myInd[[1]]+indShift),resTmpTmp}
	,{myInd,indSet}]
  ,{rankCount,Length@rank}]
  ,1];
 res=SortBy[res,#[[2]]&];
(*res=res //. gamma[x__,y_,z__]/;(!Head[y]===lVec&&!NumericQ[y])\[RuleDelayed]Sum[(gamma[x,y,z]/. y\[Rule]lInd),{lInd,0,3}];*)

    If[format==="long",
      res=res[[All,2]]
      ,
      (* only non-vanishing coefficients*)
      res=DeleteCases[res,{0,x_List},1]
    ];
    
Append[graph,"symmetrizedExpandedTensorCoeff"->res]
,{graph,Flatten@{graphs}}];

If[Length@fullResult==1,
fullResult[[1]],
fullResult
]
]



Options[processNumerator]={additionalRules -> {1->1}, spinChainSimplify->True,dimensions->4-2 eps,extractTensors->False, symCoefficients->False,coeffFormat->"long",spinSum->False, consistencyCheckLevel->1};
processNumerator[graphs_,pathToFeynmanRules_,opts:OptionsPattern[]]:=Block[{myGraphs=Flatten@{graphs},loopMom,$rules,evaluateNumerator},
	
	evaluateNumerator[num_,path_]:=Block[{v,prop,in,out},
	If[FileExistsQ[path],
		Get[path]
		,
		Print["Error: Could not find Feynman rules in: \""<>path<>"\""];
		Abort[]
		];
		Evaluate@num (*//. scalarProp[x_,y_]/;!FreeQ[x,Alternatives@@loopMom]:>1 //. scalarProp[x_,y_]/;FreeQ[x,Alternatives@@loopMom]:>1 /(sp[x,x]-y^2)*)
		];
		
		$rules=Thread@Rule[translateToFeynCalc[Keys@OptionValue[additionalRules],dimensions->OptionValue[dimensions]],translateToFeynCalc[Values@OptionValue[additionalRules],dimensions->OptionValue[dimensions]]];
		
	loopMom=myGraphs[[1,"loopMomenta"]];
	myGraphs=Map[Evaluate,#]&/@myGraphs;
	
	myGraphs=Append[#,<|"numerator_unevaluated"->Zeta[3]|>]&/@myGraphs;
	myGraphs[[All,"numerator_unevaluated"]]=myGraphs[[All,"numerator"]];
	
	
	myGraphs[[All,"numerator"]]=evaluateNumerator[myGraphs[[All,"numerator"]],pathToFeynmanRules];
	
	(* relevant for uv-cts only *)
	If[KeyExistsQ[myGraphs[[1]],"tree_level_numerator"],
		myGraphs[[All,"tree_level_numerator"]]=evaluateNumerator[myGraphs[[All,"tree_level_numerator"]],pathToFeynmanRules];
	];
	myGraphs[[All,"numerator"]]=translateToFeynCalc[#,dimensions->OptionValue[dimensions]]&/@myGraphs[[All,"numerator"]];
	myGraphs[[All,"numerator"]]=SUNFSimplify[SUNSimplify[myGraphs[[All,"numerator"]],SUNTrace->True,Factoring->True,Explicit->True]];
	myGraphs[[All,"numerator"]]=Contract[myGraphs[[All,"numerator"]]]//. $rules;
	If[OptionValue[spinChainSimplify]==True,
	 myGraphs[[All,"numerator"]]=DiracSimplify[myGraphs[[All,"numerator"]]]//.$rules
	];
	If[OptionValue[spinSum]==True,
	 myGraphs[[All,"numerator"]]=DiracSimplify[FermionSpinSum[myGraphs[[All,"numerator"]]]]//.$rules
	];
	
	
	
	myGraphs=DeleteCases[myGraphs,x_/;x[["numerator"]]==0,1];
	If[OptionValue[extractTensors]==True,
		myGraphs=extractTensCoeff[myGraphs,consistencyCheckLevel->OptionValue[consistencyCheckLevel]]
	];
	If[OptionValue[symCoefficients]==True,
		myGraphs=getSymCoeff[myGraphs,outputFormat->OptionValue[coeffFormat],consistencyCheckLevel->OptionValue[consistencyCheckLevel]]
	];
	If[Length@myGraphs==1,
	myGraphs[[1]]
	,
	myGraphs
	]
]


createSuperGraph[allGraphs_]:=Block[{myGraphs=Flatten@{allGraphs},sigs,superGraph,diagramInfo,llDiag,llSG,propPos,pows,extMom,powsAll,propsPosAll},
	If[Union[Flatten@(KeyExistsQ[#]&/@myGraphs)]!={True},
		myGraphs=getLoopLines[myGraphs];
	];
	(* signatures *)
	sigs=Flatten[myGraphs[[All,"loopLines",All,"signature"]],1]//Union;
	(* create supergraph *)
	superGraph=<|"loopLines"->Table[Evaluate/@<|"signature"-> sig ,<|"propagators"->Union@((Flatten@(Cases[myGraphs[[All,"loopLines"]],x_/;MatchQ[x,KeyValuePattern["signature"->sig]]:>KeyDrop[x[["propagators"]],{"edge_name","power"}],Infinity])))|>|>,{sig,sigs}]|>;
	diagramInfo=
		Table[
		powsAll={};
		propsPosAll={};
			Table[
				llDiag=Flatten@Cases[graph[["loopLines"]],x_/;MatchQ[x,KeyValuePattern["signature"->sigs[[sigCount]]]]:>KeyDrop[x[["propagators"]],"edge_name"],Infinity];
				pows=llDiag[[All,"power"]];
				llDiag=KeyDrop[llDiag,"power"];
				llSG=superGraph[["loopLines",sigCount,"propagators"]];
				propPos=Position[llSG,#]&/@llDiag;
				propPos=(Flatten[propPos,1]-1) /. {x_}/;IntegerQ[x]:>{sigCount-1,x};
				powsAll=Join[powsAll,pows];
				propsPosAll=Join[propsPosAll,propPos]
			,{sigCount,Length@sigs}];
			<|"denomintatorsInSuperGraph"->propsPosAll,"powers"->powsAll|>
		,{graph,myGraphs}];
	myGraphs=KeyDrop[#,"loopLines"]&/@myGraphs;
	myGraphs=Table[ AppendTo[myGraphs[[gg]],diagramInfo[[gg]]],{gg,Length@myGraphs}];
	extMom=Cases[allGraphs[[1,"momentumMap"]],Alternatives@@{_in,_out},Infinity];
	myGraphs=<|"supergraph"->Append[superGraph,<|"externalMomenta"->extMom|>],"diagrams"->myGraphs|>
]



Protect[processName,exportDirectory,exportDiagramwise,writeNumerator,processName];
Options[exportAmplitude]={processName->"new_process",exportDirectory->"./",exportDiagramwise->False,writeNumerator->False};
exportAmplitude[amp_,numerics_,opts:OptionsPattern[]]:=Block[{res,myAmp,sgName="supergraph_"<>OptionValue[processName],myOptions},
 If[KeyExistsQ[amp,"supergraph"],
    myAmp=amp,
    myAmp=createSuperGraph[amp] 
 ];
 myOptions=Thread@Rule[Keys@Options[exportAmplitude],OptionValue[Keys@Options[exportAmplitude]]];
 writeMinimalJSON[Append[myAmp[["supergraph"]],<|"name"->sgName|>],numerics,myOptions]; 
 writeAmplitudeJSON[myAmp,numerics,Join[Normal@KeyDrop[myOptions,{exportDiagramwise,writeNumerator}],{supergraphName->sgName}]];
]


Options[writeAmplitudeJSON]={exportDirectory->"./",processName->"new_Process",supergraphName->"supergraph_new_process",amplitudeType->"DiHiggsTopologyLO"};
writeAmplitudeJSON[amp_,numAssociation_Association,opts:OptionsPattern[]]:=Block[
{mygraphs,factAsso,chainAsso,ctAsso,posAsso,loopSigAsso,powAsso,nameAsso,diagCount=0,procName,denomAsso,exportDir,pLong,fullAsso,diagramwise,processAsso,expAsso,evaluateNumerator,numeratorAsso,vectorAsso,setAsso,ampTypeAsso,loopNumAsso,topoAsso,inMom,outMom,extAsso},

(* small wrapper for numerical evaluation of numerators *)	 	 
evaluateNumerator[numerator_,numRepl_]:=Block[{numericNum=numerator,Pair},
  numericNum=numericNum//FCI;
  Pair[Momentum[x_List,d___],Momentum[y_List,d___]]:=x[[1]]*y[[1]]-x[[2;;]].y[[2;;]];
  (* vectors are assumed to be contravariant *)
   numericNum=(numerator//. numRepl //.x_List[y_Integer]:>x[[y+1]]);
  If[Length@(Variables@numericNum)!=0,
  Print["Error: The numerator coefficients: "<>ToString[Variables@numericNum,InputForm]<>" have no numeric value!"];
  Abort[];  
  ];
  (* short vs long export format *)
  If[!ListQ[numericNum[[1]]],
    numericNum=((numericNum)/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]),
    numericNum[[All,1]]=(numericNum[[All,1]]/.x_/;NumericQ[x]:>ImportString[ExportString[(ReIm@x),"Real64"],"Real64"]);    
  ];
  numericNum
];


(* Create default entries *)
procName=OptionValue[processName];
exportDir=OptionValue[exportDirectory];
If[KeyExistsQ[amp,"diagrams"],
	mygraphs=Flatten@{amp[["diagrams"]]},
	Print["amplitude has no entry \"diagrams\". Can not write amplitude file."];
	Abort[];
];
If[KeyExistsQ[amp,"vectors"],
	vectorAsso=amp[[{"vectors"}]],
	vectorAsso=Evaluate/@(<|"vectors"->{{{},ConstantArray[0,4]}}|>)
	 ];
ampTypeAsso=<|"amp_type"->OptionValue["amplitudeType"]|>;
loopNumAsso=<|"n_loops"->Length@mygraphs[[1,"loopMomenta"]]|>;
topoAsso=<|"topology"->OptionValue[supergraphName]|>;

	If[KeyExistsQ[mygraphs[[1]],"momentumMap"],
		inMom=Cases[mygraphs[[1]][["momentumMap"]],_in,Infinity]//Union;
		outMom=Cases[mygraphs[[1]][["momentumMap"]],_out,Infinity]//Union
		,
		If[KeyExistsQ[mygraphs[[1]],"externalMomenta"],
			inMom=Cases[mygraphs[[1]][["externalMomenta"]],_in,Infinity]//Union;
			outMom=Cases[mygraphs[[1]][["externalMomenta"]],_out,Infinity]//Union
			,
			Print["Can not determine external momenta."];
			Abort[]	
		]	
	];
extAsso=<|"ps"->Join[inMom,outMom]/. in[x_]:>x /. out[x_]:>-x/.numAssociation|>;
(* loop over graphs *)
processAsso=Table[
    fullAsso=<||>;
	diagCount+=1;
	If[KeyExistsQ[mygraph,"name"],
		nameAsso=<|"name"->mygraph[["name"]]|>,
		nameAsso=<|"name"->procName<>"diag_"<>ToString[diagCount]|>
	];
	If[!KeyExistsQ[mygraph,"symmetrizedExpandedTensorCoeff"],
	         mygraph=getSymCoeff[mygraph]
	 ];
	numeratorAsso=<|"tensor_coefficients_split"->evaluateNumerator[mygraph[["symmetrizedExpandedTensorCoeff"]],numAssociation]|>;
	numeratorAsso=Evaluate/@numeratorAsso;
    denomAsso=<|"denominators"->mygraph[["denomintatorsInSuperGraph"]]|>;
    powAsso=<|"pows"->mygraph[["powers"]]|>;     	
    If[KeyExistsQ[mygraph,"loopSignature"],
      loopSigAsso=<|"loop_signature"->mygraph[["loopSignature"]]|>
      ,
      loopSigAsso=<|"loop_signature"->0|> 
    ];
    If[KeyExistsQ[mygraph,"positions"],
      posAsso=<|"positions"->mygraph[["positions"]]|>
      ,
      posAsso=<|"positions"->{}|> 
    ];
    If[KeyExistsQ[mygraph,"ct"],
      ctAsso=<|"ct"->mygraph[["ct"]]|>
      ,
      ctAsso=<|"ct"->False|> 
    ];
    If[KeyExistsQ[mygraph,"chain"],
      chainAsso=<|"chain"->mygraph[["chain"]]|>
      ,
      chainAsso=<|"chain"->{}|> 
    ];
    If[KeyExistsQ[mygraph,"factor"],
      factAsso=<|"factor"->mygraph[["factor"]]|>
      ,
      factAsso=<|"factor"->{1,0}|> 
    ];
	Do[fullAsso=Append[fullAsso,asso],{asso,{chainAsso,ctAsso,posAsso,loopSigAsso,denomAsso,powAsso,numeratorAsso,nameAsso,factAsso}}];
	fullAsso
	,{mygraph,mygraphs}];
	If[KeyExistsQ[amp,"sets"],
	  setAsso=amp[[{"sets"}]],
	  setAsso=Evaluate@<|"sets"->{processAsso[[All,"name"]],processAsso[[All,"name"]]}|>
	];
	processAsso={Append[<|ampTypeAsso,topoAsso,vectorAsso,setAsso,extAsso,loopNumAsso,"diagrams"->processAsso,"pols_type"->{},"mu_r_sq"->100.|>,Evaluate/@<|"name"->"amp_"<>procName|>]};
	If[DirectoryQ[exportDir],
		Export[exportDir<>"amplitude_"<>procName<>".json",processAsso,"JSON"];
		Print["Amplitude written to: "<>exportDir<>"amplitude_"<>procName<>".json"];,
		Print["Couldn't find exportDirectory. Export to standard location."];
		Export["./amplitude_"<>procName<>".json",processAsso,"JSON"];
		Print["Amplitude written to: "<>"./amplitude_"<>procName<>".json"];
	];
	
]



Protect[numericReplacement];
Options[getSymCoeffSP]={numericReplacement->{1->1}};
getSymCoeffSP[graph_,opts:OptionsPattern[]]:=Block[{myGraphs=Flatten@{graph},res,num,sp,dim,vector,Pair,loopMom,loopComponents,formatRule,indexMapping,formatRuleLight,verificationRules1,verificationRules2,loopMomMapping,coeffAsso,nonExpandedMomenta,nonExpandedComponents},
(* The output of getSymCoefficients is k^0 k^1*T_{(0,1)}. Here the output will be k^0 k^1 T^{(0,1)} *)

      
    res=Table[
    num=(mygraph[["numerator"]]//FCI//ExpandScalarProduct);
    dim=Flatten@FCGetDimensions[num];
    If[dim==={D},
          Pair[Momentum[x_,D],Momentum[y_,D]]:=sp[x,y]
		,
          Pair[Momentum[x_],Momentum[y_]]:=sp[x,y]
	];

	sp[x_List,y_List]:=x[[1]]*y[[1]]-x[[2]]*y[[2]]-x[[3]]*y[[3]]-x[[4]]*y[[4]];
	loopMom=mygraph[["loopMomenta"]];
	loopComponents=Map[Map[#,Range[0,3],{1}]&,loopMom];
	loopMomMapping=Thread@Rule[loopMom,loopComponents];

	num=(num//FCI//ExpandScalarProduct)//. Dispatch@OptionValue[numericReplacement]/. Dispatch@loopMomMapping;
	nonExpandedMomenta=Replace[#,x_/;ListQ[x]:>Nothing]&/@Flatten[(Cases[num,sp[x_,y_]:>{x,y},Infinity]//Union),1];
	nonExpandedComponents=Map[Map[#,Range[0,3],{1}]&,nonExpandedMomenta];
	num=num/.Thread@Rule[nonExpandedMomenta,nonExpandedComponents];
	
	loopComponents=Flatten@Map[Map[#,Range[0,3],{1}]&,loopMom];
	indexMapping=Thread@Rule[loopComponents,Flatten@(loopComponents/.x_[i_]:>i+4*(Position[loopMom,x]-1))];
	formatRule={Rule[a_List,b_]:>{Flatten@((Flatten@(Replace[#,(x_)^(y_/;y>1):>Flatten@{ConstantArray[x,y]}]&/@(Replace[#,1->Nothing]&/@(((loopComponents))^a))))),b}};
	formatRuleLight={Rule[a_List,b_]:> {a,b}};
	
	num=((CoefficientRules[num,loopComponents])/. formatRule)/.indexMapping;
	coeffAsso=SortBy[num,#[[1]]&];
	Append[mygraph,<|"symmetrizedExpandedTensorCoeff"->coeffAsso|>]
	,{mygraph,myGraphs}];
	If[Length@myGraphs==1,
	 res[[1]],
	 res
	]
]



constructCutInfo[graph_]:=Block[{mygraph=graph,cutInfo,cut,loopMom,extMom,leftEdgePos,rightEdgePos,cutEdgePos,cutInfoKeys,cutInfoValues,mass,basisMapping},
loopMom=mygraph[["loopMomenta"]];
extMom=Cases[mygraph[["momentumMap"]],in[x_]:>x,Infinity]//Union;
(* construct cut basis to LTD basis *)
basisMapping[cutMom_(*n-1*)]:=Block[{cutMap,replaceLoop,ltdLoopMom,basisTrafo},
cutMap=Thread@Rule[( Array[cMom,Length@cutMom]),cutMom]/. cMom[x_]:>ToExpression["cMom"~~ToString@x];
replaceLoop=(Intersection[Variables@cutMom,loopMom])[[1;;Length@cutMom]];
replaceLoop=(Solve[(cutMap/.Rule->Equal),replaceLoop])[[1]];
ltdLoopMom=Complement[loopMom,Keys@replaceLoop];
basisTrafo=Inverse@(Normal@CoefficientArrays[Flatten@{cutMom,ltdLoopMom},loopMom][[2]]);
Evaluate/@<|"cb_to_lmb"->Flatten@basisTrafo|>
];

(* cut specific information for yaml *)
Do[
cut=mygraph[["cutDiagInfo",countCut,"diagram_splitting"]];

leftEdgePos=Position[cut,"left"];
rightEdgePos=Position[cut,"right"];
cutEdgePos=Position[cut,"cut"];
cutInfoKeys={"level","m_squared","sign","name","signature"};
cutInfoValues=Table[
{
(*cut level*)
0,
(* mass *)
mass[Extract[mygraph[["particleType"]],{cutEdgePos[[cc]]}][[1]]]^2
,
 (* sign *) 
If[
MemberQ[Flatten@({Keys[#],Values[#]}&@( Extract[mygraph[["edges"]],leftEdgePos])),(Flatten@({Keys[#]}&@( Extract[mygraph[["edges"]],{cutEdgePos[[cc]]}])))[[1]] ],
1,
-1
]
,
 (* name *)
"edge"<>ToString@(cutEdgePos[[cc]][[1]] )
,
(* signature *)
 {
 Flatten@Table[(Extract[mygraph[["momentumMap"]],{cutEdgePos[[cc]]}])//. ll->1 //. Thread@Rule[Join[loopMom,extMom],ConstantArray[0,Length@loopMom+Length@extMom]],{ll,loopMom}]
,
 Join[Flatten@Table[(Extract[mygraph[["momentumMap"]],{cutEdgePos[[cc]]}])//.ee->1 //. Thread@Rule[Join[loopMom,extMom],ConstantArray[0,Length@loopMom+Length@extMom]],{ee,extMom}],ConstantArray[0,Length@extMom]]
}
}
,{cc,Length@cutEdgePos}];
cutInfo=Association/@(Thread/@Thread[Rule[ConstantArray[cutInfoKeys,Length@cutEdgePos],cutInfoValues]]);

mygraph[["cutDiagInfo",countCut]]=Append[mygraph[["cutDiagInfo",countCut]],<|"cuts"->cutInfo|>];
mygraph[["cutDiagInfo",countCut]]=Append[mygraph[["cutDiagInfo",countCut]],basisMapping[Extract[mygraph[["momentumMap"]],cutEdgePos][[1;;-2]] ]];
,{countCut,Length@mygraph[["cutDiagInfo"]]}];
mygraph
]



constructLeftRightDiagrams[graph_]:=Block[{mygraph=graph,leftDiag,rightDiag,signs,leftEdgePos,rightEdgePos,cutEdgePos,cut,cutEdge,compParaShift,cutSol,constructLeft,constructRight,loopLines,props,loopEdges,loopProps,freeProps,loopSig,getLTDSqLL,
leftName,rightName},
(* small wraper to get LTD-loop lines plus ltd-cut-structure *)
getLTDSqLL[subdiag_,parashift_]:=Block[{subgraph=subdiag,eMom,lMom,lLines,lProp,orientation,res,mass,pos},
		If[!KeyExistsQ[subgraph,"loopMomenta"],
			Print["Error: Graphs have no entry \"loopMomenta\". Import with importGraph[]"];
		Abort[]
		];
		eMom=Complement[Variables@(Union@Flatten@subgraph[["momentumMap"]]),subgraph[["loopMomenta"]]];

	lLines=CoefficientArrays[((subgraph[["momentumMap"]] /. Thread[Rule[eMom,0]]/. 0:>Nothing)^2 //FullSimplify)/.(x_)^2:>x// DeleteDuplicates,subgraph[["loopMomenta"]]][[2]]//Normal;
	

		lProp=
	Table[
			orientation=(ReplaceAll[#,x_/;!NumberQ[x]:>0]&/@(FullSimplify@((subgraph[["momentumMap"]]/. Thread[Rule[eMom,0]])/(Dot[subgraph[["loopMomenta"]],ll])) ));
	  
              pos=Position[Abs@orientation,1,Heads->False];
	      Association@{
	"end_node"->"99",
			"propagators"->Association/@Table[
	{
		                "mass"->mass[subgraph[["particleType",pp]]]^2,
      "name"-> If[
		KeyExistsQ[subgraph,"name"],
		subgraph[["name",pp]],
		"prop"<>ToString@pp
	  ],
     "parametric_shift"->orientation[[pp]]*parashift[[pp]] ,
      "power"->If[KeyExistsQ[subgraph,"power"], subgraph[["power",pp]],  1  ],
     "q"->(ConstantArray[0,4])}
	   ,{pp,Flatten@pos}],
"signature"->ll
,"start_node"->"99"
}
	,{ll,lLines}];

	subgraph=AppendTo[subgraph,<|"loopLines"->lProp|>];

getCutStructure[subgraph]
          ];
(* -------------------------------------------------------------------------------------------------------------- *)
(* function to compute looplines: cutMom momenta of n-1 cut edges *)
loopLines[diag_,loopMom_,cutMom_]:=Block[{mydiag=diag,cutMap,cMom,replaceLoop,extMom,masses,mass,fullProp,freeEdges,loopGraph,ltdLoopMom,parashift,loopMomMap},
extMom=Complement[Variables@(mydiag[["momentumMap"]]/. in[x_]:>x /. out[x_]:>x),loopMom];
cutMap=Thread@Rule[( Array[cMom,Length@cutMom]),cutMom]/. cMom[x_]:>ToExpression["cMom"~~ToString@x];
replaceLoop=(Intersection[Variables@cutMom,loopMom])[[1;;Length@cutMom]];
replaceLoop=(Solve[(cutMap/.Rule->Equal),replaceLoop])[[1]];
mydiag[["momentumMap"]]=mydiag[["momentumMap"]]/. replaceLoop;

(* edges free of any ltd-loop-momentum *)
freeEdges=Position[(mydiag[["momentumMap"]]),x_/;(FreeQ[Variables@x,Alternatives@@(loopMom)]&&Head[x]=!=in && Head[x]=!=out),1,Heads->False];
loopEdges=Position[(mydiag[["momentumMap"]]),x_/;!FreeQ[Variables@x,Alternatives@@(loopMom)],1,Heads->False];
parashift={Join[Coefficient[#,Keys@cutMap],{0}],Join[Coefficient[#,extMom],ConstantArray[0,Length@extMom]]}&/@(diag[["momentumMap"]](*/. in[x_]\[RuleDelayed]Nothing /. out[x_]\[RuleDelayed]Nothing*));
loopMomMap={loopMom,ConstantArray[0,2*Length@extMom]};

If[Length@(Flatten@loopEdges)>0,
loopGraph=mydiag;
ltdLoopMom=Complement[loopMom,Keys@replaceLoop];
loopGraph=AppendTo[loopGraph,Association@("loopMomenta"->ltdLoopMom)];
loopGraph=getLTDSqLL[loopGraph,parashift];
loopMomMap=loopMomMap/. Thread[Rule[ltdLoopMom,ConstantArray[1,Length@ltdLoopMom]]]
];

 loopMomMap=loopMomMap/. Thread[Rule[loopMom,ConstantArray[0,Length@loopMom]]];
If[Union@(Flatten@loopMomMap)==={0},
loopMomMap={}
];
masses=mass[#]^2&/@(mydiag[["particleType"]]);


If[(Length@(Flatten@loopEdges)>0)&&(Length@(Flatten@freeEdges)>0),
freeProps=Table[Association@(Association/@{("m_squared"->Extract[masses,cE]),("name"->Extract[mydiag[["name"]],cE]),("parametric_shift"->Extract[parashift,cE]),
("power"->1),"q"->ConstantArray[0,4]}),{cE,freeEdges}];
mydiag=Append[mydiag,Association@("loopMomenta"->ltdLoopMom)];
mydiag=Association@("graphInfo"->mydiag);
mydiag=AppendTo[mydiag,Association@("cut_momentum_map"->replaceLoop)];
mydiag=AppendTo[mydiag,Association@("loop_lines"->Join[{Association@{"end_node"->99,"propagators"->freeProps,"signature"->ConstantArray[0,Length@ltdLoopMom],"start_node"->99}},loopGraph[["loopLines"]]])];
mydiag=AppendTo[mydiag,Association@("ltd_cut_structure"->Flatten/@(Transpose@{ConstantArray[0,Length@loopGraph[["cutStructure"]]],loopGraph[["cutStructure"]]})) ];
mydiag=AppendTo[mydiag,Association@("loop_momentum_map"->loopMomMap)]
,
If[Length@(Flatten@loopEdges)>0,
(* only loops. not relevant yet *)
mydiag=Append[mydiag,Association@("loopMomenta"->ltdLoopMom)];
mydiag=Association@("graphInfo"->mydiag);
mydiag=AppendTo[mydiag,Association@("cut_momentum_map"->replaceLoop)];
mydiag=AppendTo[mydiag,Association@("loop_lines"->loopGraph[["loopLines"]])];
mydiag=AppendTo[mydiag,Association@("ltd_cut_structure"->loopGraph[["cutStructure"]]) ];
mydiag=AppendTo[mydiag,Association@("loop_momentum_map"->loopMomMap)]
,
(* only free edges *)
freeProps=Table[Association@(Association/@{("m_squared"->Extract[masses,cE]),("name"->Extract[mydiag[["name"]],cE]),("parametric_shift"->Extract[parashift,cE]),
("power"->1),"q"->ConstantArray[0,4]}),{cE,freeEdges}];
mydiag=Append[mydiag,Association@("loopMomenta"->{})];
mydiag=Association@("graphInfo"->mydiag);
mydiag=AppendTo[mydiag,Association@("cut_momentum_map"->replaceLoop)];
mydiag=AppendTo[mydiag,Association@("loop_lines"->{Association@{"end_node"->99,"propagators"->freeProps,"signature"->{},"start_node"->99}})];
mydiag=AppendTo[mydiag,Association@("ltd_cut_structure"->{{0}}) ];
mydiag=AppendTo[mydiag,Association@("loop_momentum_map"->loopMomMap)]
]
]


];
(* -------------------------------------------------------------------------------------------------------------- *)
(* -------------------------------------------------------------------------------------------------------------- *)
(* small wrappers to counstruct left- and right diagrams *)
constructLeft[asso_]:=Block[{},Association@( asso->Join[Extract[mygraph[[asso]],leftEdgePos],
Table[
If[asso==="edges",
If[
signs[[cP]]==-1,
	Extract[mygraph[[asso]],cutEdgePos[[cP]]] /. Rule[a_,b_]:>Rule[a,out[cutEdge[cP]]]
,
	Extract[mygraph[[asso]] ,cutEdgePos[[cP]]]/. Rule[a_,b_]:>Rule[in[cutEdge[cP]],b ]
]
,
If[
signs[[cP]]==-1,
	out@Extract[mygraph[[asso]],cutEdgePos[[cP]]]
,
	in@Extract[mygraph[[asso]] ,cutEdgePos[[cP]]]
]
]
,{cP,Length@cutEdgePos}]
])
];
(* -------------------------------------------------------------------------------------------------------------- *)
constructRight[asso_]:=Block[{},Association@(asso->Join[Extract[mygraph[[asso]],rightEdgePos],
Table[
If[asso==="edges",
If[
signs[[cP]]==-1,
	Extract[mygraph[[asso]],cutEdgePos[[cP]]] /. Rule[a_,b_]:>Rule[in[cutEdge[cP]],b]
,
	Extract[mygraph[[asso]] ,cutEdgePos[[cP]]]/. Rule[a_,b_]:>Rule[a,out[cutEdge[cP]]]
]
,
If[
signs[[cP]]==-1,
	in@Extract[mygraph[[asso]],cutEdgePos[[cP]]]
,
	out@Extract[mygraph[[asso]] ,cutEdgePos[[cP]]]
]
]
,{cP,Length@cutEdgePos}]
])
];
(* -------------------------------------------------------------------------------------------------------------- *)
(* -------------------------------------------------------------------------------------------------------------- *)

(* ------------------ INPUT Checks-------------------------------- *)
(* check if cutInfo exists *)
If[KeyExistsQ[mygraph[["cutDiagInfo",1]],"cuts"],
mygraph=mygraph,
mygraph=constructCutInfo[mygraph];
];



(* ------------------ Definitions Cut Specific -------------------------------- *)
(* cut specific data *)
Do[
cut=mygraph[["cutDiagInfo",cutCount,"diagram_splitting"]];
signs=mygraph[["cutDiagInfo",cutCount,"cuts",All,"sign"]];

cutEdgePos=Position[cut,"cut"];

leftEdgePos=Position[cut,"left"];
rightEdgePos=Position[cut,"right"];

leftDiag=Association@(constructLeft/@{"edges","momentumMap","particleType"});
leftName="edge"<>ToString[#]&/@(Flatten@({leftEdgePos,cutEdgePos}));
rightName="edge"<>ToString[#]&/@(Flatten@({rightEdgePos,cutEdgePos}));
leftDiag=AppendTo[leftDiag,Association@("name"->leftName)];
rightDiag=Association@(constructRight/@{"edges","momentumMap","particleType"});
rightDiag=AppendTo[rightDiag,Association@("name"->rightName)];

leftDiag=Association@{
Association@("external_kinematics"->{}(*ConstantArray[{0,0,0,0},Length@(graph[["edges"]])]*)),
loopLines[leftDiag,mygraph[["loopMomenta"]],Extract[mygraph[["momentumMap"]],cutEdgePos][[1;;-2]] ],
Association@("name"->"left_diag_cut_"<>StringReplace[ToString@(Flatten@cutEdgePos),{"{"->"e",","->"e","}"->""," "->""}])
};
rightDiag=Association@{
Association@("external_kinematics"->{}(*ConstantArray[{0,0,0,0},Length@(graph[["edges"]])]*)),
loopLines[rightDiag,mygraph[["loopMomenta"]],Extract[mygraph[["momentumMap"]],cutEdgePos][[1;;-2]] ],
Association@("name"->"right_diag_cut_"<>StringReplace[ToString@(Flatten@cutEdgePos),{"{"->"e",","->"e","}"->""," "->""}])
};

mygraph[["cutDiagInfo",cutCount]]=AppendTo[mygraph[["cutDiagInfo",cutCount]],Association@("left_diagram"->leftDiag)];
mygraph[["cutDiagInfo",cutCount]]=AppendTo[mygraph[["cutDiagInfo",cutCount]],Association@("right_diagram"->rightDiag)];
,{cutCount,1,Length@mygraph["cutDiagInfo"]}];
mygraph
]


writeLTDSqrtJSON[graph_,numAssociation_,procName_:"new_process"]:=Block[{mygraph=graph,uvCTAsso,ccCutAsso,cutAsso,fullAsso={},uvBuild,evaluateNumerator,superG,ecm},
(* -------------------------------------------------- *)
(* small wrapper for numerical evaluation of numerators *)	 	 
evaluateNumerator[numerator_,numRepl_]:=Block[{numericNum=numerator,Pair},
  numericNum=numericNum//FCI;
  Pair[Momentum[x_List,d___],Momentum[y_List,d___]]:=x[[1]]*y[[1]]-x[[2;;]].y[[2;;]];
  (* vectors are assumed to be contravariant *)
   numericNum=(numerator//. numRepl //.x_List[y_Integer]:>x[[y+1]]);
  If[Length@(Variables@numericNum)!=0,
  Print["Error: The numerator coefficients: "<>ToString[Variables@numericNum,InputForm]<>" have no numeric value!"];
  Abort[];  
  ];
  (* short vs long export format *)
  If[!ListQ[numericNum[[1]]],
           Print["Error: numerator should be in sparse format... Recomputing numerator"];
numericNum=False;    
  ];
  numericNum
];
(* -------------------------------------------------- *)

(* small wrapper for extracting uv-limits *)
uvBuild[cInfoGraphs_]:=Block[{preface,diagAsso,tensorAsso,symCoeff},
preface={
cInfoGraphs[[{"cb_to_lmb"}]],
Association@("conjugate_deformation"->{False,True})
};
diagAsso=Flatten/@Association@("diagrams"->
{
Association@(Association/@{
cInfoGraphs[["left_diagram",{"external_kinematics"}]],
Map[Evaluate,#,-1]&/@(ReplaceAll[#,numAssociation]&/@cInfoGraphs[["left_diagram",{"loop_lines"}]]),
cInfoGraphs[["left_diagram",{"ltd_cut_structure"}]],
cInfoGraphs[["left_diagram",{"loop_momentum_map"}]],
"maximum_ratio_expansion_threshold"->-1,
"n_loops"->Length@cInfoGraphs[["left_diagram","graphInfo","loopMomenta"]],
cInfoGraphs[["left_diagram",{"name"}]],
"numerator_tensor_coefficients"->{{0,0}}
})
,
Association@(Association/@{
cInfoGraphs[["right_diagram",{"external_kinematics"}]],
Map[Evaluate,#,-1]&/@(ReplaceAll[#,numAssociation]&/@cInfoGraphs[["right_diagram",{"loop_lines"}]]),
cInfoGraphs[["right_diagram",{"ltd_cut_structure"}]],
cInfoGraphs[["right_diagram",{"loop_momentum_map"}]],
"maximum_ratio_expansion_threshold"->-1,
"n_loops"->Length@cInfoGraphs[["right_diagram","graphInfo","loopMomenta"]],
cInfoGraphs[["right_diagram",{"name"}]],
"numerator_tensor_coefficients"->ConstantArray[{0,0},1(*15*)]
})
});
If[KeyExistsQ[cInfoGraphs,"numerator"],
If[KeyExistsQ[cInfoGraphs,"symmetrizedExpandedTensorCoeff"],
symCoeff=evaluateNumerator[cInfoGraphs[["symmetrizedExpandedTensorCoeff"]],numAssociation];
If[symCoeff===False,
symCoeff=(getSymCoeffSP[cInfoGraphs,numericReplacement->numAssociation])[["symmetrizedExpandedTensorCoeff"]]
]
,
symCoeff=(getSymCoeffSP[cInfoGraphs,numericReplacement->numAssociation])[["symmetrizedExpandedTensorCoeff"]]
]
,
If[KeyExistsQ[mygraph,"symmetrizedExpandedTensorCoeff"],
symCoeff=evaluateNumerator[mygraph[["symmetrizedExpandedTensorCoeff"]],numAssociation];
If[symCoeff===False,
symCoeff=(getSymCoeffSP[mygraph,numericReplacement->numAssociation])[["symmetrizedExpandedTensorCoeff"]]
]
,
symCoeff=(getSymCoeffSP[mygraph,numericReplacement->numAssociation])[["symmetrizedExpandedTensorCoeff"]]
],
symCoeff=(getSymCoeffSP[mygraph,numericReplacement->numAssociation])[["symmetrizedExpandedTensorCoeff"]];
];
symCoeff[[All,2]]=((ReIm@symCoeff[[All,2]])/. {x_,y__}/;NumericQ[x]:>ImportString[ExportString[{x,y},"Real64"],"Real64"]);
symCoeff=Association@("numerator_tensor_coefficients_sparse"->symCoeff);
Association@("uv_limits"->{Association@(Flatten@({preface,diagAsso,symCoeff,Association@("symmetry_factor"->1)}))})
];
(* -------------------------------------------------- *)
(* Write JSON *)
(* -------------------------------------------------- *)
If[!KeyExistsQ[mygraph,"symmetrizedExpandedTensorCoeff"],
mygraph=getSymCoeffSP[mygraph,numericReplacement->numAssociation]
];
(* Loop over all cutkosky cuts *)
fullAsso=Association@("cutkosky_cuts"->Flatten@Table[
Association@{Map[Evaluate,#,-1]&/@(ReplaceAll[#,numAssociation]&/@cutInfo[[{"cuts"}]]),Association@("n_bubbles"->0),uvBuild[cutInfo]}
,{cutInfo,mygraph[["cutDiagInfo"]]}
]);
superG=getLoopLines[mygraph,True];
superG=Association@("topo"->writeMinimalJSON[superG,numAssociation,{export->False,processName->procName}][[1]]);
ecm=((Total@(Cases[mygraph[["momentumMap"]],in[x_]:>x]//Union))/.numAssociation)/. x_List:>x[[1]]*x[[1]]-x[[2;;]].x[[2;;]];
Export[procName<>".json",
Association@{Association@("MG_numerator"-><||>),
fullAsso,
Association@("external_momenta"->superG[["topo","external_kinematics"]]),
Association@("e_cm_squared"->ecm),
Association@("n_incoming_momenta"->Count[mygraph[["particleType"]],_in]),
Association@("n_loops"->Length[mygraph[["loopMomenta"]]]),
Association@("name"->superG[["topo","name"]]),
Evaluate/@Association@("overall_numerator"->If[KeyExistsQ[mygraph,"overall_numerator"],mygraph[["overall_numerator"]],1]),
Association@("numerator_in_loop_momentum_basis"->True),
superG}
,"JSON"]
]

