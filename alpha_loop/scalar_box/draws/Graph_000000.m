GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=Labeled[GraphClass[{
Labeled[Style[CreateEdge["2","3",3],Red,Thickness[0.002000]],"g|q1|p1"],
Labeled[Style[CreateEdge["1","4",1],Red,Thickness[0.002000]],"g|q2|p2"],
Labeled[Style[CreateEdge["9","8",12],Red,Thickness[0.002000]],"g|q3|p1"],
Labeled[Style[CreateEdge["10","7",10],Red,Thickness[0.002000]],"g|q4|p2"],
Labeled[Style[CreateEdge["4","3",5],Blue,Thickness[0.005000]],"t|pq1|k1"],
Labeled[Style[CreateEdge["3","5",6],Blue,Thickness[0.002000]],"t|pq2|k1+p1"],
Labeled[Style[CreateEdge["6","4",7],Blue,Thickness[0.002000]],"t|pq3|k1-p2"],
Labeled[Style[CreateEdge["5","6",8],Blue,Thickness[0.002000]],"t|pq4|k1+p1-k3"],
Labeled[Style[CreateEdge["9","10",14],Blue,Thickness[0.005000]],"t|pq5|k2"],
Labeled[Style[CreateEdge["11","9",15],Blue,Thickness[0.002000]],"t|pq6|k2+p1"],
Labeled[Style[CreateEdge["10","12",16],Blue,Thickness[0.002000]],"t|pq7|k2-p2"],
Labeled[Style[CreateEdge["12","11",17],Blue,Thickness[0.002000]],"t|pq8|k2+p1-k3"],
Labeled[Style[CreateEdge["5","11",2],Green,Thickness[0.005000]],"h|pq9|k3"],
Labeled[Style[CreateEdge["6","12",0],Green,Thickness[0.002000]],"h|pq10|p2+p1-k3"]
},
EdgeShapeFunction -> {
CreateEdge["2","3",3]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["1","4",1]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["9","8",12]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","7",10]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["4","3",5]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["3","5",6]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["6","4",7]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","6",8]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["9","10",14]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["11","9",15]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["10","12",16]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["12","11",17]->GraphElementData["Arrow", "ArrowSize" -> 0.015000],
CreateEdge["5","11",2]->GraphElementData["HalfFilledDoubleArrow", "ArrowSize" -> 0.025000],
CreateEdge["6","12",0]->GraphElementData["Arrow", "ArrowSize" -> 0.015000]
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> 10, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> {"PackingLayout"->"ClosestPacking"},
ImageSize -> {660.000000, 510.000000}
],"MG: SG_QG0 | FORM: #0"];
Export["Graph_000000.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {825.000000, 637.500000}];