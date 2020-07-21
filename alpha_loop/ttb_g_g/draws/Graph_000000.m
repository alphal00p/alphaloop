GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=Labeled[GraphClass[{
Labeled[Style[CreateEdge["2","3",3],Blue,Thickness[0.002000]],"t|q1|p1"],
Labeled[Style[CreateEdge["1","3",1],Blue,Thickness[0.002000]],"\!\(\*OverscriptBox[\(t\), \(_\)]\)|q2|p2"],
Labeled[Style[CreateEdge["7","5",8],Blue,Thickness[0.002000]],"\!\(\*OverscriptBox[\(t\), \(_\)]\)|q3|p2"],
Labeled[Style[CreateEdge["7","6",10],Blue,Thickness[0.002000]],"t|q4|p1"],
Labeled[Style[CreateEdge["4","3",5],Red,Thickness[0.002000]],"g|pq1|-p1-p2"],
Labeled[Style[CreateEdge["4","7",2],Red,Thickness[0.005000]],"g|pq2|k1"],
Labeled[Style[CreateEdge["4","7",0],Red,Thickness[0.002000]],"g|pq3|p2+p1-k1"]
},
EdgeShapeFunction -> {
CreateEdge["2","3",3]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["1","3",1]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","5",8]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["7","6",10]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","3",5]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","7",2]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["4","7",0]->GraphElementData["Arrow", "ArrowSize" -> 0.015000]
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> {"PackingLayout"->"ClosestPacking"},
ImageSize -> {660.000000, 510.000000}
],"MG: SG_QG0 | FORM: #0"];
Export["Graph_000000.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];