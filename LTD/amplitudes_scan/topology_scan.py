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
from amplitudes import params, AmplitudesCollection,Amplitude, qqbar_diagram
from generate_1L_mom import *



def ddAAA_amp(topology_name, amplitude_name):
    tree_factor = params['alpha_ew']**1.5*params['q_d']**3 * (4.0 * np.pi)**1.5
    # Initialize
    amp = Amplitude(topology_name,                    # name of topology
                    'qqbar_photons',                  # process type
                    zip(["u", "vbar", "a", "a", "a"], # polarizations
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
 
def topology_uv(qs,ms,topology_name):
    points = len(qs)
    #Create the Graph for the topology
    # pi: propagators, qi: externals
    graph = [('p%d'%(i+1),i+1,((i+1) % points)+1) for i in range(points)]
    graph.extend([('q%d'%(i+1),i+101,i+1) for i in range(points)])

    top_uv = TopologyGenerator(graph)

    return top_uv.create_loop_topology(
            topology_name, 
            ext_mom={'q%d'%(i+1): qs[i] for i in range(points)}, 
            mass_map={'p%d'%(i+1): ms[i] for i in range(points)}, 
            loop_momenta_names= ('p%d'%points,), # If not specified an arbitrary spanning tree will be used for momentum routing 
            analytic_result=None, # For triangle and box one-loop topology, the analytic result is automatically computed
            fixed_deformation = None,
         )

#with open("../topologies.yaml", 'r') as stream:
#    try:
#        print("Importing topologies from yaml")
#        topologies = {top['name']: top for top in yaml.safe_load(stream)}
#    except yaml.YAMLError as exc:
#        print(exc)
#

topo_base_name = "ddAAA_amplitude"
amp_base_name = "ddAAA"
if len(sys.argv) == 2:
    root_dir = sys.argv[1]
else:
    root_dir = ''


topology_collection = TopologyCollection()
amplitudes_collection = AmplitudesCollection()
topo_names = []
ms = [0.,1e2,.0,.0,.0,.0]

for seed in range(10):
    topo_name = topo_base_name + "_S%d"%seed
    topo_names += [topo_name]
    qs = generate_moms([0.]*5,3,seed)
    qs[2:] = list(itertools.permutations(qs[2:]))[np.random.randint(5)]
    qs.insert(1,np.array([0.]*4))
    print(qs)
    topology_collection.add_topology(topology_uv(qs,ms,topo_name), entry_name = topo_name)
topology_collection.export_to(os.path.join(root_dir,'topologies.yaml'))

for topo_name in topo_names:
    amp_name = amp_base_name + topo_name[-3:]
    amplitudes_collection.add(ddAAA_amp(topo_name,amp_name))
amplitudes_collection.export_to(os.path.join(root_dir,'amplitudes.yaml'))


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