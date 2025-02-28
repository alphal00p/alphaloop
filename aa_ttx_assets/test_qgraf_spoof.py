graphs=[]


graphs.append(
{
"edges":{
(-1+4):{
    "name":"p1",
    "PDG": 22,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
(-3+4):{
    "name":"p2",
    "PDG": 22,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },
(-2+4):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p1",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },
(-4+4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p2",
    "indices": (-4+4+1,),
    "vertices":(4+4,-4+4+1)
 },
(1+4):{
    "name":"q"+str(1),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (2+4+1,1+4+1,),
    "vertices":(5+4,1+4)
 },
(2+4):{
    "name":"q"+str(2),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (4+4+1,3+4+1,),
    "vertices":(1+4,6+4)
 },
(3+4):{
    "name":"q"+str(3),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k2",
    "indices": (6+4+1,5+4+1,),
    "vertices":(2+4,5+4)
 },
(4+4):{
    "name":"q"+str(4),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k2-p2",
    "indices": (8+4+1,7+4+1,),
    "vertices":(7+4,2+4)
 },
(5+4):{
    "name":"q"+str(5),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k3",
    "indices": (10+4+1,9+4+1,),
    "vertices":(6+4,3+4)
 },
(6+4):{
    "name":"q"+str(6),
    "PDG": 6,
    "type": "virtual",
    "momentum": "-k3-p1",
    "indices": (12+4+1,11+4+1,),
    "vertices":(3+4,8+4)
 },
(7+4):{
    "name":"q"+str(7),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1+k2-k3-p1-p2",
    "indices": (14+4+1,13+4+1,),
    "vertices":(4+4,7+4)
 },
(8+4):{
    "name":"q"+str(8),
    "PDG": 6,
    "type": "virtual",
    "momentum": "k1+k2-k3-p1",
    "indices": (16+4+1,15+4+1,),
    "vertices":(8+4,4+4)
 },
(9+4):{
    "name":"q"+str(9),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-k1-k2",
    "indices": (18+4+1,17+4+1,),
    "vertices":(8+4,5+4)
 },
(10+4):{
    "name":"q"+str(10),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-k3-p1",
    "indices": (20+4+1,19+4+1,),
    "vertices":(7+4,6+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (22,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1,
    "edge_ids": (-1+4,)
 },
 -3+4+1:{
    "PDGs": (22,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1,
    "edge_ids": (-3+4,)
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p1"),
    "indices": (-2+4+1,),
    "vertex_id": -2,
    "edge_ids": (-2+4,)
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p2"),
    "indices": (-4+4+1,),
    "vertex_id": -2,
    "edge_ids": (-4+4,)
 },
 1+4:{
    "PDGs": (-6,22,6),
    "momenta": ("k1-p1","p1","-k1"),
    "indices": (4+4+1,-1+4+1,1+4+1),
    "vertex_id": 0,
    "edge_ids": (2+4,-1+4,1+4)
    },
 2+4:{
    "PDGs": (-6,22,6),
    "momenta": ("-k2","p2","k2-p2"),
    "indices": (6+4+1,-3+4+1,7+4+1),
    "vertex_id": 0,
    "edge_ids": (3+4,-3+4,4+4)
    },
 3+4:{
    "PDGs": (-6,22,6),
    "momenta": ("k3+p1","-p1","-k3"),
    "indices": (12+4+1,-2+4+1,9+4+1),
    "vertex_id": 0,
    "edge_ids": (6+4,-2+4,5+4)
    },
 4+4:{
    "PDGs": (-6,22,6),
    "momenta": ("-k1-k2+k3+p1+p2","-p2","k1+k2-k3-p1"),
    "indices": (14+4+1,-4+4+1,15+4+1),
    "vertex_id": 0,
    "edge_ids": (7+4,-4+4,8+4)
    },
 5+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k1","-k1-k2","k2"),
    "indices": (2+4+1,17+4+1,5+4+1),
    "vertex_id": 0,
    "edge_ids": (1+4,9+4,3+4)
    },
 6+4:{
    "PDGs": (-6,21,6),
    "momenta": ("k3","k1-k3-p1","-k1+p1"),
    "indices": (10+4+1,19+4+1,3+4+1),
    "vertex_id": 0,
    "edge_ids": (5+4,10+4,2+4)
    },
 7+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k2+p2","-k1+k3+p1","k1+k2-k3-p1-p2"),
    "indices": (8+4+1,20+4+1,13+4+1),
    "vertex_id": 0,
    "edge_ids": (4+4,10+4,7+4)
    },
 8+4:{
    "PDGs": (-6,21,6),
    "momenta": ("-k1-k2+k3+p1","k1+k2","-k3-p1"),
    "indices": (16+4+1,18+4+1,11+4+1),
    "vertex_id": 0,
    "edge_ids": (8+4,9+4,6+4)
    }
},
"overall_factor": "25"
}
)
