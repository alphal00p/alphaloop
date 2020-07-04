graphs=[]


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(4+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(2+4,1+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(1+4,3+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1-p2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1-p3",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1-p1","p1","-k1"),
    "indices": (4+4+1,-1+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-1+4,1+4)
    },
 2+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1","p2","-k1-p2"),
    "indices": (2+4+1,-3+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-3+4,3+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1-p1+p3","-p3","-k1+p1"),
    "indices": (8+4+1,-2+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-2+4,2+4)
    },
 4+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1+p2","-p4","-k1+p1-p3"),
    "indices": (6+4+1,-4+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-4+4,4+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(4+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(1+4,2+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(3+4,1+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1+p2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(2+4,4+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1+p3",
    "indices": (8+4+1,7+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1","p1","k1-p1"),
    "indices": (2+4+1,-1+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-1+4,2+4)
    },
 2+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1-p2","p2","k1"),
    "indices": (6+4+1,-3+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-3+4,1+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p1","-p3","k1-p1+p3"),
    "indices": (4+4+1,-2+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-2+4,4+4)
    },
 4+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p1-p3","-p4","k1+p2"),
    "indices": (8+4+1,-4+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-4+4,3+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(4+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(2+4,1+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(1+4,3+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1-p2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1-p4",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1-p1","p1","-k1"),
    "indices": (4+4+1,-1+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-1+4,1+4)
    },
 2+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1","p2","-k1-p2"),
    "indices": (2+4+1,-3+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-3+4,3+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1-p1+p4","-p4","-k1+p1"),
    "indices": (8+4+1,-4+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-4+4,2+4)
    },
 4+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1+p2","-p3","-k1+p1-p4"),
    "indices": (6+4+1,-2+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-2+4,4+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(4+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(1+4,2+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(3+4,1+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1+p2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(2+4,4+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1+p4",
    "indices": (8+4+1,7+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1","p1","k1-p1"),
    "indices": (2+4+1,-1+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-1+4,2+4)
    },
 2+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1-p2","p2","k1"),
    "indices": (6+4+1,-3+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-3+4,1+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p1","-p4","k1-p1+p4"),
    "indices": (4+4+1,-4+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-4+4,4+4)
    },
 4+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p1-p4","-p3","k1+p2"),
    "indices": (8+4+1,-2+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-2+4,3+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,4+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(2+4,1+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(1+4,3+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p3",
    "indices": (6+4+1,5+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1-p4",
    "indices": (8+4+1,7+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1-p1","p1","-k1"),
    "indices": (4+4+1,-1+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-1+4,1+4)
    },
 2+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1","-p3","-k1+p3"),
    "indices": (2+4+1,-2+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-2+4,3+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("k1-p1+p4","-p4","-k1+p1"),
    "indices": (8+4+1,-4+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-4+4,2+4)
    },
 4+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1-p3","p2","-k1+p1-p4"),
    "indices": (6+4+1,-3+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-3+4,4+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


graphs.append(
{
"edges":{
(-1+4):{
    "name":"q1",
    "PDG": 21,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"q2",
    "PDG": 21,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,4+4)
 },
(-2+4):{
    "name": "q3",
    "PDG": 25,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },
(-4+4):{
    "name": "q4",
    "PDG": 25,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(1+4):{
    "name":"p"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(1+4,2+4)
 },
(2+4):{
    "name":"p"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(3+4,1+4)
 },
(3+4):{
    "name":"p"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p3",
    "indices": (6+4+1,5+4+1,),
    "vertices":(2+4,4+4)
 },
(4+4):{
    "name":"p"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1-p1+p4",
    "indices": (8+4+1,7+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (21,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (21,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (25,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (25,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1","p1","k1-p1"),
    "indices": (2+4+1,-1+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,-1+4,2+4)
    },
 2+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p3","-p3","k1"),
    "indices": (6+4+1,-2+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-2+4,1+4)
    },
 3+4:{
    "PDGs": (-6,25,6),
    "momenta": ("-k1+p1","-p4","k1-p1+p4"),
    "indices": (4+4+1,-4+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-4+4,4+4)
    },
 4+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1+p1-p4","p2","k1-p3"),
    "indices": (8+4+1,-3+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,-3+4,3+4)
    }
},
"overall_factor": "-1",
"in_momenta":("p1","p2"),
"out_momenta":("p3","p4")
}
)


