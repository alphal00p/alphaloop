import yaml
import os
import sys
import math
import yaml
import itertools

pjoin = os.path.join

file_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(file_path,os.pardir))
sys.path.insert(0, pjoin(file_path,os.pardir,os.pardir))

import vectors
import ltd_commons
from ltd_utils import TopologyGenerator, TopologyCollection
#from generate_1L_mom import *
import numpy as np
import generate_1L_mom

def topology_pentabox(qs,ms,topology_name):
    points = len(qs)
    props = len(ms)
    
    factory = TopologyGenerator([
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])
    
    print(qs[0]+qs[1])
    print([np.sqrt(generate_1L_mom.sp2(q)) for q in qs] )
    print(sum(q for q in qs))
     
    return factory.create_loop_topology(
            topology_name,
            ext_mom={'q%d' % (i+1): qs[i] for i in range(points)},
            mass_map={'p%d' % (i+1): ms[i] for i in range(props)},
            loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing
            analytic_result=None
         )

def Pentabox_physical_massive(topolgy_name):
    topology_name = "T1_Pentabox_physical_bis_massive"
    rescaling = 1./3.
    # This topology correspond to a two-loop pentabox with physical kinematics
    q1 = vectors.LorentzVector(
        [0.149500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,    0.149165176901313E+01]
    )*rescaling
    q2 = vectors.LorentzVector(
        [0.150500000000000E+01,    0.000000000000000E+00,    0.000000000000000E+00,   -0.149165176901313E+01]
    )*rescaling
    q3 = vectors.LorentzVector(
        [-0.126041949101381e+01,    -0.452362952912639e+00,    -0.101350243653045e+01,   0.516563513332600e+00]
    )*rescaling
    q4 = vectors.LorentzVector(
        [-0.105098730574850e+01,   0.489324061520790e-01,   0.928212188578101e+00,    -0.283905035967510e+00]
    )*rescaling
    q5 = -q4-q3-q2-q1
    qs = [q1,q2,q3,q4,q5]
    points = len(qs)
    print(qs)
    print(np.sqrt((q1+q2).square()))
    print(np.sqrt((q4+q5).square()))
    print(np.sqrt(q1.square()))
    print(np.sqrt(q2.square()))
    print(np.sqrt(q3.square()))
    print(np.sqrt(q4.square()))
    print(np.sqrt(q5.square()))
    print(q1+q2+q3+q4+q5)
    factory = TopologyGenerator([
        ('p1', 1, 6), ('p2', 6, 7), ('p3', 7, 2), ('p4', 2, 1),
        ('p5', 7, 3), ('p6', 3, 4), ('p7', 4, 5), ('p8', 5, 6),
        ('q1', 101,1), ('q2', 102,2), ('q3', 103,3), ('q4', 104,4), ('q5', 105,5),
    ])
    return factory.create_loop_topology(
            topology_name,
            ext_mom={'q1':q1, 'q2': q2 , 'q3': q3, 'q4': q4, 'q5': q5},
            mass_map={'p1': 0.35, 'p2': 0.35, 'p3': 0.35, 'p4': 0.35, 'p5': 0.35, 'p6': 0.35, 'p7': 0.35, 'p8': 0.35},
            loop_momenta_names=('p4','p5'), # If not specified an arbitrary spanning tree will be used for momentum routing
            analytic_result=None
         )



def ddAA_amp(topolgy_name, amplitude_name):
    tree_factor = params['alpha_ew'] * params['q_d']**2 * (4.0 * np.pi)
    # Initialize
    amp = Amplitude(topolgy_name,                # name of topology
                    'qqbar_photons',             # process type
                    zip(["u", "vbar", "a", "a"], # polarizations
                        ["+", "-", "+", "+"]),
                    1,                           # uv_pos
                    91.188)                      # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / amp.sij['s23']
    amp.add_born([6, 5, 7], born_factor)
    # Diagrams and vectors
    factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
    amp.create_amplitude(
        amplitude_name,
        1,
        [qqbar_diagram("D1", [0, 3, 4], [1, 1, 1], [6, 5, -1, 3, 7, 4, -1], factor/amp.sij["s23"], False, True,),
         qqbar_diagram("D2", [0, 2, 3], [1, 1, 1], [-1, 2, 6, 3, -1, 5, 7], factor/amp.sij["s23"], False, True),
         qqbar_diagram("D3", [0, 2, 3, 4], [1, 1, 1, 1], [-1, 2, 6, 3, 7, 4, -1], factor, False, False,),
         qqbar_diagram("D4", [0, 3], [1, 1], [6, 5, -1, 3, -1, 5, 7], factor/amp.sij["s23"]/amp.sij["s23"], False, True),
         qqbar_diagram("IR", [0, 2, 4], [1, 1, 1], [-1, 2, 6, 5, 7, 4, -1], factor / amp.sij["s23"], True, True),
         ],
        # Vectors [loopmomenta, externals]
        [
            [[-1], -amp.ps[0]],
            [[-1], -amp.ps[0]-amp.ps[1]],
            [[-1], amp.ps[3]],
            [[-1], vectors.LorentzVector([0, 0, 0, 0])],
            [[], -amp.ps[1]-amp.ps[2]],
        ],
    )
    return amp

def ddAAA_amp(topology_name, amplitude_name):
    tree_factor = params['alpha_ew']**1.5*params['q_d']**3 * (4.0 * np.pi)**1.5
    # Initialize
    amp = Amplitude(topology_name,                    # name of topology
                    'qqbar_photons',                  # process type
                    zip(["u", "vbar", "a", "a", "a"],  # polarizations
                        ["+", "-", "+", "+", "+"]),
                    1,                                # uv_pos
                    91.188)                           # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / amp.sij['s23'] / amp.sij['s15']
    amp.add_born([8, 6, 9, 7, 10], born_factor)
    # Diagrams and vectors
    factor = -1j * params['C_F'] * tree_factor\
        * params['alpha_s'] * (4.0 * np.pi)
    amp.create_amplitude(
        amplitude_name,
        1,
        [qqbar_diagram("D1", [0, 3, 4, 5], [1, 1, 1, 1], [8, 6, -1, 3, 9, 4, 10, 5, -1], factor/amp.sij["s23"], False, False),
         qqbar_diagram("D2", [0, 4, 5], [1, 1, 1], [8, 6, 9, 7, -1, 4, 10, 5, -1], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D3", [0, 4], [1, 1], [8, 6, 9, 7, -1, 4, -1, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"]**2, False, True,),
         qqbar_diagram("D4", [0, 3, 4], [1, 1, 1], [8, 6, -1, 3, 9, 4, -1, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D5", [0, 2, 3, 4], [1, 1, 1, 1], [-1, 2, 8, 3, 9, 4, -1, 7, 10], factor/amp.sij["s15"], False, False,),
         qqbar_diagram("D6", [0, 2, 3], [1, 1, 1], [-1, 2, 8, 3, -1, 6, 9, 7, 10], factor/amp.sij["s23"]/amp.sij["s15"], False, True,),
         qqbar_diagram("D7", [0, 3], [1, 1], [8, 6, -1, 3, -1, 6, 9, 7, 10], factor/amp.sij["s23"]**2/amp.sij["s15"], False, True,),
         qqbar_diagram("D8", [0, 2, 3, 4, 5], [1, 1, 1, 1, 1], [-1, 2, 8, 3, 9, 4, 10, 5, -1], factor, False, False,),
         qqbar_diagram("IR", [0, 2, 5], [1, 1, 1], [-1, 2, 8, 6, 9, 7, 10, 5, -1], factor/amp.sij["s23"]/amp.sij["s15"], True, True,),
         ],
        # Vectors [loopmomenta, externals]
        [
            [[-1], -amp.ps[0]],
            [[-1], -amp.ps[0]-amp.ps[1]],
            [[-1], amp.ps[3]+amp.ps[4]],
            [[-1], amp.ps[4]],
            [[-1], vectors.LorentzVector([0, 0, 0, 0])],
            [[], -amp.ps[1]-amp.ps[2]],
            [[], amp.ps[0]+amp.ps[4]],
        ],
    )
    return amp


def topology_uv(qs, ms, topology_name):
    points = len(qs)
    # Create the Graph for the topology
    # pi: propagators, qi: externals
    graph = [('p%d' % (i+1), i+1, ((i+1) % points)+1) for i in range(points)]
    graph.extend([('q%d' % (i+1), i+101, i+1) for i in range(points)])

    top_uv = TopologyGenerator(graph)

    return top_uv.create_loop_topology(
        topology_name,
        ext_mom={'q%d' % (i+1): qs[i] for i in range(points)},
        mass_map={'p%d' % (i+1): ms[i] for i in range(points)},
        # If not specified an arbitrary spanning tree will be used for momentum routing
        loop_momenta_names=('p%d' % points,),
        # For triangle and box one-loop topology, the analytic result is automatically computed
        analytic_result=None,
        fixed_deformation=None,
    )

# with open("../topologies.yaml", 'r') as stream:
#    try:
#        print("Importing topologies from yaml")
#        topologies = {top['name']: top for top in yaml.safe_load(stream)}
#    except yaml.YAMLError as exc:
#        print(exc)
#


#topo_base_name = "ddAAA_amplitude"
#amp_base_name = "ddAAA"
#generate_case = "ddAAA"
#ms = [0., 1e2, .0, .0, .0, .0]
topo_base_name = "ddAA"
amp_base_name = "ddAA"
generate_scan = "ddAA"
ms = [0., 1e2, .0, .0, .0]

topo_base_name = "PBPM"
generate_scan = "pentabox_physical_massive"
ms=[.35]*8


if len(sys.argv) == 2:
    root_dir = sys.argv[1]
else:
    root_dir = ''

generate_amplitudes = False
#Pentabox_physical_massive("T1")
topology_collection = TopologyCollection()
topo_names = []
madgraph_input = ""
pysecdec_input = ""
file_mg5 = open("mg5_inputs.txt", "w")
file_psd = open("pySecDec_inputs.txt", "w")

#for seed in range(10):
#    topo_name = topo_base_name + "_S%d" % seed
#    topo_names += [topo_name]
#    qs = generate_moms_random([0.]*5, 3, seed)
#    qs[2:] = list(itertools.permutations(qs[2:]))[np.random.randint(5)]
#    # Store input to evaluate with MG5
#    madgraph_input += "#ID{}\n".format(seed)
#    madgraph_input += output_moms_string(qs, "fortran")
#    qs.insert(1, np.array([0.]*4))
#    print(qs)
#    topology_collection.add_topology(topology_uv(qs, ms, topo_name), entry_name=topo_name)
if generate_scan == "pentabox_physical_massive":
    top_id = -1
    for s45 in np.linspace(0,1,27)[1:-1]:
        for theta in np.linspace(0,np.pi,27)[1:-1]:
            top_id+=1
            topo_name = topo_base_name + "_S%d"%top_id
            topo_names += [topo_name]
            qs = generate_1L_mom.generate_moms_2to3_scan([(i/30.0)**2 for i in range(1,6)],s45, theta, 0, np.pi/6.0,3)
            #Store input to evaluate with MG5
            if qs is None:
                print(topo_name)
            else:
                qs = 3.0 * np.array(qs)
                pysecdec_input += "#ID{}\n".format(top_id)
                pysecdec_input += generate_1L_mom.output_moms_string(qs,"python")
                #print(qs)
                topology_collection.add_topology(topology_pentabox(qs,ms,topo_name), entry_name = topo_name)
elif generate_scan == "ddAAA":
    top_id = 0
    for s45 in np.linspace(0,1,27)[1:-1]:
        for theta in np.linspace(0,np.pi,27)[1:-1]:
            topo_name = topo_base_name + "_S%d"%top_id
            topo_names += [topo_name]
            qs = generate_1L_mom.generate_moms_2to3_scan([0.]*5,s45, theta, 0 *np.pi/3.0, np.pi/6.0,3)
            #Store input to evaluate with MG5
            if qs is None:
                print(topo_name)
                #continue        
            madgraph_input += "#ID{}\n".format(top_id)
            madgraph_input += generate_1L_mom.output_moms_string(qs,"fortran")
            #Add UV prop to the topology
            qs.insert(1,np.array([0.]*4))
            print(qs)
            topology_collection.add_topology(topology_uv(qs,ms,topo_name), entry_name = topo_name)
            top_id+=1
elif generate_scan == "ddAA":
    top_id = 0
    for theta in np.linspace(0,np.pi,102)[1:-1]:
        topo_name = topo_base_name + "_S%d"%top_id
        topo_names += [topo_name]
        qs = [
            [0.5,0.,0.0,0.5],
            [0.5,0.,0.0,-0.5],
            [-0.5,0.,0.5*np.sin(theta),0.5*np.cos(theta)],
            [-0.5,0.,-0.5*np.sin(theta),-0.5*np.cos(theta)],
        ] 
        madgraph_input += "#ID{}\n".format(top_id)
        madgraph_input += generate_1L_mom.output_moms_string(qs,"fortran")
        #Add UV prop to the topology
        qs.insert(1,np.array([0.]*4))
        topology_collection.add_topology(topology_uv(qs,ms,topo_name), entry_name = topo_name)
        top_id+=1
else:
    print("No such option for generating the scan! (%s)"%generate_scan)
    sys.exit(1)

file_mg5.write(madgraph_input)
file_mg5.close()
file_psd.write(pysecdec_input)
file_psd.close()
topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))

if generate_amplitudes:
    from amplitudes import params, AmplitudesCollection,Amplitude, qqbar_diagram
    amplitudes_collection = AmplitudesCollection()
    for (i,topo_name) in enumerate(topo_names):
        amp_name = amp_base_name + "_S%d"%i
        if generate_scan == "2to3":
            amplitudes_collection.add(ddAAA_amp(topo_name, amp_name))
        elif generate_scan == "2to2":
            amplitudes_collection.add(ddAA_amp(topo_name, amp_name))
    amplitudes_collection.export_to(os.path.join(root_dir, 'amplitudes.yaml'))


'''
if output_path is not None:
            open(output_path, 'w').write(
                yaml.dump(flat_record, Dumper=Dumper))
        else:
            return yaml.dump(flat_record, Dumper=Dumper)

for top in topologies:
    if top['name'] == topology:
        mytop = top
        break
else:
    raise AssertionError("Could not find topology %s" % topology)
self.topology_name = topology
self.type = amp_type
self.ps = [vectors.LorentzVector(p)
           for p in mytop['external_kinematics']]
self.uv_pos = uv_pos
   
'''
