graphs=[]


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (1+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(1+4,-2+4+1)
 },

    (2+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(2+4,-4+4+1)
 },
(1+4,3+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "p1-p3",
    "indices": (1+4+1,2+4+1,),
    "vertices":(1+4,3+4)
 },
(4+4,2+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-p2+p4",
    "indices": (3+4+1,4+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4,3+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p1+p3",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,3+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-p1+p3","-p3","p1"),
    "indices": (2+4+1,-2+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("p2","-p4","-p2+p4"),
    "indices": (-3+4+1,-4+4+1,3+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("-k1","k1-p1+p3","p1-p3"),
    "indices": (8+4+1,5+4+1,1+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2-p4","-k1+p1-p3","k1"),
    "indices": (4+4+1,6+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (2+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },

    (1+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(1+4,-4+4+1)
 },
(1+4,3+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "p1-p4",
    "indices": (1+4+1,2+4+1,),
    "vertices":(1+4,3+4)
 },
(4+4,2+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-p2+p3",
    "indices": (3+4+1,4+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4,3+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p1+p4",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,3+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-p1+p4","-p4","p1"),
    "indices": (2+4+1,-4+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("p2","-p3","-p2+p3"),
    "indices": (-3+4+1,-2+4+1,3+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("-k1","k1-p1+p4","p1-p4"),
    "indices": (8+4+1,5+4+1,1+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2-p3","-k1+p1-p4","k1"),
    "indices": (4+4+1,6+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,1+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },

    (2+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },

    (3+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(4+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (1+4+1,2+4+1,),
    "vertices":(4+4,1+4)
 },
(3+4,2+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(3+4,2+4)
 },
(2+4,4+4,3+4):{
    "name":"q"+str(3),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1-p3",
    "indices": (5+4+1,6+4+1,),
    "vertices":(2+4,4+4)
 },
(4+4,3+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p4",
    "indices": (7+4+1,8+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","-p1-p2","p1"),
    "indices": (-3+4+1,1+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1+p3","-p3","-k1"),
    "indices": (6+4+1,-2+4+1,3+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1","-p4","-k1+p4"),
    "indices": (4+4+1,-4+4+1,7+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("k1-p4","p1+p2","-k1-p3"),
    "indices": (8+4+1,2+4+1,5+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "-1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,1+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },

    (2+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(2+4,-2+4+1)
 },

    (3+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(4+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-p1-p2",
    "indices": (1+4+1,2+4+1,),
    "vertices":(4+4,1+4)
 },
(2+4,3+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(2+4,3+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1+p3",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1-p4",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","-p1-p2","p1"),
    "indices": (-3+4+1,1+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-k1","-p3","k1+p3"),
    "indices": (4+4+1,-2+4+1,5+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-k1+p4","-p4","k1"),
    "indices": (8+4+1,-4+4+1,3+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("-k1-p3","p1+p2","k1-p4"),
    "indices": (6+4+1,2+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "-1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (1+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(1+4,-2+4+1)
 },

    (3+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(1+4,4+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "p1-p3",
    "indices": (1+4+1,2+4+1,),
    "vertices":(1+4,4+4)
 },
(3+4,2+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(3+4,2+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p2",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4,3+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p4",
    "indices": (7+4+1,8+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-p1+p3","-p3","p1"),
    "indices": (2+4+1,-2+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","k1-p2","-k1"),
    "indices": (-3+4+1,5+4+1,3+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1","-p4","-k1+p4"),
    "indices": (4+4+1,-4+4+1,7+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("k1-p4","-k1+p2","p1-p3"),
    "indices": (8+4+1,6+4+1,1+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (3+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },

    (1+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(1+4,-4+4+1)
 },
(1+4,4+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "p1-p4",
    "indices": (1+4+1,2+4+1,),
    "vertices":(1+4,4+4)
 },
(3+4,2+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(3+4,2+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p2",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(4+4,3+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p3",
    "indices": (7+4+1,8+4+1,),
    "vertices":(4+4,3+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-p1+p4","-p4","p1"),
    "indices": (2+4+1,-4+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","k1-p2","-k1"),
    "indices": (-3+4+1,5+4+1,3+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1","-p3","-k1+p3"),
    "indices": (4+4+1,-2+4+1,7+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("k1-p3","-k1+p2","p1-p4"),
    "indices": (8+4+1,6+4+1,1+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,2+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,2+4)
 },
 (-3+4+1,1+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },

    (1+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(1+4,-2+4+1)
 },

    (3+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(4+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-p2+p3",
    "indices": (1+4+1,2+4+1,),
    "vertices":(4+4,1+4)
 },
(2+4,3+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(2+4,3+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p1",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1-p4",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("p2","-p3","-p2+p3"),
    "indices": (-3+4+1,-2+4+1,1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("-k1","k1-p1","p1"),
    "indices": (4+4+1,5+4+1,-1+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-k1+p4","-p4","k1"),
    "indices": (8+4+1,-4+4+1,3+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2-p3","-k1+p1","k1-p4"),
    "indices": (2+4+1,6+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,2+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,2+4)
 },
 (-3+4+1,1+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,1+4)
 },

    (3+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },

    (1+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(1+4,-4+4+1)
 },
(4+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-p2+p4",
    "indices": (1+4+1,2+4+1,),
    "vertices":(4+4,1+4)
 },
(2+4,3+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(2+4,3+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": 21,
    "type": "virtual",
    "momentum": "k1-p1",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "k1-p3",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("p2","-p4","-p2+p4"),
    "indices": (-3+4+1,-4+4+1,1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("-k1","k1-p1","p1"),
    "indices": (4+4+1,5+4+1,-1+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("-k1+p3","-p3","k1"),
    "indices": (8+4+1,-2+4+1,3+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2-p4","-k1+p1","k1-p3"),
    "indices": (2+4+1,6+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (3+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(3+4,-2+4+1)
 },

    (4+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(4+4,-4+4+1)
 },
(2+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (1+4+1,2+4+1,),
    "vertices":(2+4,1+4)
 },
(1+4,3+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(1+4,3+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1-p2",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p1-p3",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("k1-p1","-k1","p1"),
    "indices": (4+4+1,1+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","k1","-k1-p2"),
    "indices": (-3+4+1,2+4+1,5+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1-p1+p3","-p3","-k1+p1"),
    "indices": (8+4+1,-2+4+1,3+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1+p2","-p4","-k1+p1-p3"),
    "indices": (6+4+1,-4+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


graphs.append(
{
"edges":{
 (-1+4+1,1+4,1):{
    "name":"p1",
    "PDG": +2,
    "type": "in",
    "momentum": "p1",
    "indices": (-1+4+1,),
    "vertices": (-1+4+1,1+4)
 },
 (-3+4+1,2+4,2):{
    "name":"p2",
    "PDG": -2,
    "type": "in",
    "momentum": "p2",
    "indices": (-3+4+1,),
    "vertices": (-3+4+1,2+4)
 },

    (4+4,-2+4+1,3):{
    "name": "p3",
    "PDG": 22,
    "type": "out",
    "momentum": "p3",
    "indices": (-2+4+1,),
    "vertices":(4+4,-2+4+1)
 },

    (3+4,-4+4+1,4):{
    "name": "p4",
    "PDG": 22,
    "type": "out",
    "momentum": "p4",
    "indices": (-4+4+1,),
    "vertices":(3+4,-4+4+1)
 },
(2+4,1+4,1+4):{
    "name":"q"+str(1),
    "PDG": 21,
    "type": "virtual",
    "momentum": "-k1",
    "indices": (1+4+1,2+4+1,),
    "vertices":(2+4,1+4)
 },
(1+4,3+4,2+4):{
    "name":"q"+str(2),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p1",
    "indices": (3+4+1,4+4+1,),
    "vertices":(1+4,3+4)
 },
(4+4,2+4,3+4):{
    "name":"q"+str(3),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1-p2",
    "indices": (5+4+1,6+4+1,),
    "vertices":(4+4,2+4)
 },
(3+4,4+4,4+4):{
    "name":"q"+str(4),
    "PDG": +2,
    "type": "virtual",
    "momentum": "-k1+p1-p4",
    "indices": (7+4+1,8+4+1,),
    "vertices":(3+4,4+4)
 },

},
"nodes": {
 -1+4+1:{
    "PDGs": (+2,),
    "momenta": ("p1"),
    "indices": (-1+4+1,),
    "vertex_id": -1
 },
 -3+4+1:{
    "PDGs": (-2,),
    "momenta": ("p2"),
    "indices": (-3+4+1,),
    "vertex_id": -1
 },
 -2+4+1:{
    "PDGs": (22,),
    "momenta": ("p3"),
    "indices": (-2+4+1,),
    "vertex_id": -2
 },
 -4+4+1:{
    "PDGs": (22,),
    "momenta": ("p4"),
    "indices": (-4+4+1,),
    "vertex_id": -2
 },
 1+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("k1-p1","-k1","p1"),
    "indices": (4+4+1,1+4+1,-1+4+1),
    "vertex_id": 0
    },
 2+4:{
    "PDGs": (-2,21,+2),
    "momenta": ("p2","k1","-k1-p2"),
    "indices": (-3+4+1,2+4+1,5+4+1),
    "vertex_id": 0
    },
 3+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1-p1+p4","-p4","-k1+p1"),
    "indices": (8+4+1,-4+4+1,3+4+1),
    "vertex_id": 0
    },
 4+4:{
    "PDGs": (-2,22,+2),
    "momenta": ("k1+p2","-p3","-k1+p1-p4"),
    "indices": (6+4+1,-2+4+1,7+4+1),
    "vertex_id": 0
    }
},
"overall_factor": "+1"
}
)


