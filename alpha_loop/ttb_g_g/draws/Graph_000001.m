GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=Labeled[GraphClass[{
Labeled[Style[CreateEdge["1","4",2],Blue,Thickness[0.002000]],"\!\(\*OverscriptBox[\(t\), \(_\)]\)|q1|p2"],
Labeled[Style[CreateEdge["2","3",4],Blue,Thickness[0.002000]],"t|q2|p1"],
Labeled[Style[CreateEdge["7","5",8],Blue,Thickness[0.002000]],"\!\(\*OverscriptBox[\(t\), \(_\)]\)|q3|p2"],
Labeled[Style[CreateEdge["7","6",10],Blue,Thickness[0.002000]],"t|q4|p1"],
Labeled[Style[CreateEdge["3","4",5],Blue,Thickness[0.002000]],"t|pq1|p1-k1"],
Labeled[Style[CreateEdge["4","7",1],Red,Thickness[0.002000]],"g|pq2|p2+p1-k1"],
Labeled[Style[CreateEdge["3","7",3],Red,Thickness[0.005000]],"g|pq3|k1"]
},
EdgeShapeFunction -> {
CreateEdge["1","4",2]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["2","3",4]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","5",8]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","6",10]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["3","4",5]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","7",1]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["3","7",3]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000]
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> {"PackingLayout"->"ClosestPacking"},
ImageSize -> {660.000000, 510.000000}
],"MG: SG_QG1 | FORM: #1"];
Export["Graph_000001.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];