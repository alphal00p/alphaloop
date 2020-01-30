#!/usr/bin/env python3

# This file contains the description of some aplitudes
# import topologies
import vectors
import os
import sys
import numpy as np
import copy
import yaml
from ltd_utils import TopologyCollection

pjoin = os.path.join
file_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, pjoin(file_path, "amplitudes/"))

from amplitudes_topologies import amplitude_topologies, plot_edge_info

with open("topologies.yaml", 'r') as stream:
    try:
        print("Importing topologies from yaml for the creation of the amplitudes...")
        topologies = {top['name']: top for top in yaml.safe_load(stream)}
        print("done")
    except yaml.YAMLError as exc:
        print(exc)


params = {
    'g_f':   1.166390e-05,
    'alpha_s':   1.180000e-01,
    'alpha_ew':   1./1.325070e+02,
    'C_F':   4./3.,
    'q_u': 2. / 3.,
    'q_d': -1. / 3.,
}


class dihiggs_diagram(object):
    def __init__(self, name, propagators, pows, factor, ct, uv, coeffs=[]):
        """
            name:     given name of the diagram, usefull for debugging
            dens:     list with all the denominator for this diagram
            pows:     power of each propagator
            chain:    gamma chain contraction
            coeffs:   coefficient list for the tensor structure of the numerator
            position: position of the vectors that depend on the loop momentum
            sing:     sign rising from taking the derivative wrt k0
            uv:       True/False if needs uv counterterms
            ct:       True/False if it is a counterterm and therefore comes with a minus sing
        """
        self.name = name
        self.dens = propagators
        self.pows = pows
        self.chain = []
        self.coeffs = coeffs
        self.positions = []
        self.signs = 0
        print(factor)
        self.factor = complex(factor)
        print(self.factor)
        self.uv = uv
        self.ct = ct

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        res = {}

        res["name"] = str(self.name)
        res["denominators"] = [list(ll) for ll in self.dens]
        #res["denominators"] = self.dens
        res["pows"] = self.pows
        res["chain"] = self.chain
        res["tensor_coefficients_split"] = [[float(c.real), float(c.imag)] for c in self.coeffs]
        res["positions"] = self.positions
        res["loop_signature"] = self.signs
        res["factor"] = [float(self.factor.real), float(self.factor.imag)]
        res["ct"] = self.ct
        return res

class qqbar_diagram(object):
    def __init__(self, name, propagators, pows, chain, factor, ct, uv, positions=[], coeffs=[]):
        """
            name:     given name of the diagram, usefull for debugging
            dens:     list with all the denominator for this diagram
            pows:     power of each propagator
            chain:    gamma chain contraction
            coeffs:   coefficient list for the tensor structure of the numerator
            position: position of the vectors that depend on the loop momentum
            sing:     sign rising from taking the derivative wrt k0
            uv:       True/False if needs uv counterterms
            ct:       True/False if it is a counterterm and therefore comes with a minus sing
        """
        self.name = name
        self.dens = propagators
        self.pows = pows
        self.chain = chain
        self.coeffs = coeffs
        self.positions = positions
        self.signs = 0
        self.factor = complex(factor)
        self.uv = uv
        self.ct = ct

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        res = {}

        res["name"] = str(self.name)
        res["denominators"] = [list(ll) for ll in self.dens]
        #res["denominators"] = self.dens
        res["pows"] = self.pows
        res["chain"] = self.chain
        res["tensor_coefficients_split"] = self.coeffs
        res["positions"] = self.positions
        res["loop_signature"] = self.signs
        res["factor"] = [float(self.factor.real), float(self.factor.imag)]
        res["ct"] = self.ct
        return res



class AmplitudesCollection(dict):
    def add(self, amplitude, entry_name=None):
        assert(isinstance(amplitude, Amplitude))
        if entry_name == None:
            if amplitude.name == None:
                raise BaseException(
                    "Specify the name of the hard-coded amplitude to add to the collection.")
            self[amplitude.name] = amplitude
        else:
            self[entry_name] = amplitude

    def export_to(self, output_path, format='yaml'):
        """ Exports this amplitude to a given format."""

        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Amplitude can only be exported in the following formats: %s" % (
                ', '.join(allowed_export_format)))

        if export_format == 'yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException(
                    "Install yaml python module in order to export amplitudes to yaml.")

        flat_record = [self[k].to_flat_format() for k in sorted(self.keys())]
        if output_path is not None:
            open(output_path, 'w').write(
                yaml.dump(flat_record, Dumper=Dumper))
        else:
            return yaml.dump(flat_record, Dumper=Dumper)


def hard_coded_ddAAA(amp_name, topo_name):
    tree_factor = params['alpha_ew']**1.5*params['q_d']**3 * (4.0 * np.pi)**1.5
    # Initialize
    amp = Amplitude(topo_name,                        # name of topology
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
        amp_name,
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
    return(amp)

def hard_coded_ddAA(amp_name, topo_name):
    tree_factor = params['alpha_ew'] * params['q_d']**2 * (4.0 * np.pi)
    # Initialize
    amp = Amplitude(topo_name,                   # name of topology
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
        amp_name,
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
    return(amp)

def hard_coded_dd4A(amp_name, topo_name):
    tree_factor = params['alpha_ew']**2 * params['q_d']**4 * (4.0 * np.pi)**2
    # Initialize
    amp = Amplitude(topo_name,                   # name of topology
                    'qqbar_photons',             # process type
                    zip(["u", "vbar", "a", "a", "a", "a"],  # polarizations
                        ["+", "-", "+", "+", "+", "+"]),
                    1,                           # uv_pos
                    91.188)                      # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / (amp.qSQ(2)*amp.qSQ(3)*amp.qSQ(4))
    amp.add_born([10, 7, 11, 8, 12, 9, 13], born_factor)
    # Diagrams and vectors
    factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
    
    vector_list=[ [[-1], sum(-amp.ps[i] for i in range(i+1))] for i in range(5) ] # add k+p slashed from prop
    vector_list.extend([[[-1], vectors.LorentzVector([0, 0, 0, 0])]])
    vector_list.extend([ [[], sum(-amp.ps[i] for i in range(1,i+1))] for i in range(2,5) ]) # add k+p when k is soft
    
    amp.create_amplitude(
       amp_name,
       1,
       generate_qqbar_diags(amp,4,factor),
       # Vectors [loopmomenta, externals]
       vector_list)
    return(amp)
   
def generate_qqbar_photons(external_momenta, ms, amp_name, topo_name, amp_type, amp_top=None):
    # Initialize
    n = len(external_momenta)-3
    pols_name = ["u","vbar"]
    pols_name.extend(["a"]*n)
    pols_type = ["+","-"]
    pols_type.extend(["+"]*n)
    
    if amp_type == 'Nf_qqbarphotonsNNLO':
        amp_top = amplitude_topologies.create(amp_type, external_momenta, 1e2, topo_name, fixed_deformation=False)
        amp = generate_qqbar_photons(amp_top.qs[:-2], amp_top.ms[:-2], amp_name, topo_name, 'qqbarphotonsNLO', amp_top=amp_top)
        
        
        # Overwrite topology infos
        amp_init = Amplitude(topo_name, amp_type, external_momenta, 1e2)
        amp.amp_top = amp_init.amp_top
        amp.amp_top.get_edge_info()
        
        #Add Bubble ant its UV counterterms
        diags = copy.deepcopy(amp.diags)
        diagsCT = copy.deepcopy(amp.diags)
        
        new_den_n = len(external_momenta) - 1
        
        for diag in diags:
            if diag.name == 'born':
                continue
            factor = 1
            diag.dens +=[amp.amp_top.edge_info['p%d' % (new_den_n+1)]['loopline'],
                         amp.amp_top.edge_info['p%d' % (new_den_n+3)]['loopline']]
            diag.pows += [1, 1]
            diag.factor *= factor
        for diag in diagsCT:
            if diag.name == 'born':
                continue
            factor = 1
            diag.dens +=[amp.amp_top.edge_info['p%d' % new_den_n]['loopline']]
            diag.pows += [2]
            diag.ct = not diag.ct
            diag.name = "UV2" + diag.name
            diag.factor *= factor
        
        for subset in amp.sets:
            subset += ['UV2'+s for s in subset]

        amp.diags = diags + diagsCT
#        for diag in amp.diags:
#            print(diag.to_flat_format())

        return(amp)

    elif amp_type == 'qqbarphotonsNLO':
        amp = Amplitude(topo_name,                 # name of topology
                    'qqbarphotonsNLO',           # process type
                    external_momenta,
                    ms,
                    zip(pols_name,pols_type),  # polarizations
                    0,                         # uv_pos
                    91.188,
                    amp_top=amp_top)                    # mu_r_sq

        print(external_momenta)
        print(ms)
        # Vectors
        vector_list_loop_dep = []  # momenta flowing in each propagator 
        vector_list_loop_ind = []  # same but at Tree-Level
        # With this shift it's possible to recover the momentum flowing in the Tree-Level diagram
        mom_shift = [prop.q for prop in amp.mytop.loop_lines[0].propagators if prop.name == 'p1'][0]
        for ll in amp.mytop.loop_lines:
            for i, p in enumerate(ll.propagators):
                if i == 0 or p.name not in ['p%d' %i for i in range(len(external_momenta)+1)]: # UV prop
                    continue
                vector_list_loop_dep += [[[-1 * s for s in ll.signature], -p.q]]
                if n+2 > i > 2 :
                    vector_list_loop_ind += [[[], -p.q + mom_shift]]
        
        vector_list = vector_list_loop_dep + vector_list_loop_ind
    
        #Factors
        tree_factor = params['alpha_ew']**(n/2.) * params['q_d']**n * (4.0 * np.pi)**(n/2.)
        NLO_factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
        
        # Born Level Diagram
        born_factor = tree_factor
        for _, v in vector_list_loop_ind:
            born_factor /= v.square()
    
        amp.add_born(generate_qqbar_born_chain(n), born_factor)
        
        # NLO Level Diagrams
        amp.create_amplitude(
           amp_name,
           0,
           generate_qqbar_diags(amp,n,NLO_factor),
           vector_list)
        
        return(amp)
    else:
        raise Exception("NO SUCH AMPLITUDE TYPE! %s" % amp_type)
    

class Amplitude(object):
    def __init__(self, topology, amp_type, moms, ms, polarizations=None, uv_pos=-1, mu_r_sq=1e2, amp_top=None):
        """
            uv_pos:       position of the uv propagator
            mu_r_sq:  mu renormalization for the CT
        """

        #with open("topologies.yaml", 'r') as stream:
        #    try:
        #        topologies = yaml.safe_load(stream)
        #    except yaml.YAMLError as exc:
        #        print(exc)
        #for top in topologies:
        #    if top['name'] == topology:
        #        mytop = top
        #        break
        #else:
        #    raise AssertionError("Could not find topology %s" % topology)
        if amp_top is None:
            try:
                #amplitude_topologies.build_topology
                self.amp_top = amplitude_topologies.create(amp_type, moms, ms, topology, fixed_deformation=False)
                self.amp_top.get_edge_info()
                #plot_edge_info(self.amp_top.edge_info)
                self.mytop = self.amp_top.topology
            except:
                raise AssertionError("Could not create AmplitudeTopology for %s" %amp_type)
        else:
            self.amp_top = amp_top
            self.amp_top.get_edge_info()
            self.mytop = self.amp_top.topology

        self.topology_name = topology
        self.n_loops = self.mytop.n_loops
        self.type = amp_type
        self.ps = moms.copy()
        self.uv_pos = uv_pos
        if uv_pos >= 0:
            del self.ps[uv_pos+1]
        self.compute_invariants(polarizations)
        self.mu_r_sq = mu_r_sq
        self.diags = []

    def compute_invariants(self, polarizations):
        self.sij = {}
        for i, pi in enumerate(self.ps):
            for j, pj in enumerate(self.ps):
                self.sij['s%d%d' % (i+1, j+1)] = (pi+pj).square()

        # Polaizatoin vector contains a vector with all the polarization kinds
        # None, [vbar,+/-], [v,+/-], [ubar,+/-], [u,+/-], [a,+/-],
        if polarizations == None:
            self.polarizations = []
        else:
            self.polarizations = []
            for pol in polarizations:
                if any(pol[0] == p for p in ["u", "ubar", "v", "vbar", "a"]):
                    if any(pol[1] == s for s in ["+", "-"]):
                        self.polarizations += [pol[0]+pol[1]]
                    else:
                        raise BaseException(
                            "No polarization known with this sign %s" % str(pol[1]))
                else:
                    raise BaseException(
                        "No polarization known with this name %s" % pol[0])

    def qSQ(self, n):
        if self.type == 'qqbarphotonsNLO':
            return sum(p for p in self.ps[1:n+1]).square()
        else:
            raise BaseException(
                "Not knonw diagram structure for %s" % self.type)

    def add_born(self, chain, factor):
        if self.type == 'qqbarphotonsNLO':
            diag = qqbar_diagram("born", [], [], chain, factor, False, False)
            diag.positions = []
            diag.signs = 0
            self.diags += [diag]
        else:
            raise BaseException(
                "Not knonw diagram structure for %s" % self.type)

    def create_amplitude(self, name, loops, diags, vectors):
        self.name = name
        self.vectors = vectors
        # Find position in chain of loop dependent pieces
        loop_dependent = [
            i+1 for (i, v) in enumerate(vectors) if len(v[0]) > 0]
        
        for diag in diags:
            (positions, v_ids) = np.array(
                [(pos,v_id) for (pos, v_id) in enumerate(diag.chain) if any(v_id == i for i in loop_dependent)]
                ).transpose()
            diag.positions = positions.tolist()
            diag.signs = vectors[v_ids[0]-1][0][0]
        
        # Add new diagram to the amplitude
        self.diags.extend(diags)

        if self.type == 'DiHiggsTopologyLO':
            # Define diagrams
            dens = [ [0, 1, 4, 5], [0, 1, 2, 5], [0, 2, 3, 5]]
            coeffs = np.load('amplitudes/dihiggs_coeffs.npy')
            diags = []
            for i in range(3):
                diags += [dihiggs_diagram("D%d" % i ,dens[i],[1,1,1,1], 1, False,False, coeffs=coeffs[i])]
            self.diags = diags
            # Map denominator to LoopLines tuples
            for diag in self.diags:
                diag.dens = [self.amp_top.edge_info['p%d' % (d+1)]['loopline'] for d in diag.dens]

            # Global to the amplitude
            self.sets = [['D0', 'D1', 'D2'], ['D0', 'D1', 'D2']]

            vector_default = [[], [0,0,0,0]]
            self.vectors = [vector_default]



        elif self.type == 'qqbarphotonsNLO':
            # Create UV CT
            for diag in diags:
                if diag.uv and self.uv_pos >= 0:
                    # Remove Leading UV pole
                    uv_diag = copy.deepcopy(diag)
                    uv_diag.name = "UV" + uv_diag.name
                    uv_diag.dens = [self.uv_pos]
                    uv_diag.ct = not uv_diag.ct
                    # Rise power of the unique propagator
                    uv_diag.pows = [len(diag.pows)]
                    for i in uv_diag.positions:
                        uv_diag.chain[i] = self.uv_pos+1
                    self.diags += [uv_diag]

                    # Check if a bubble
                    if len(diag.dens) == 2:
                        pair = [uv_diag.chain[uv_diag.positions[0]],
                                uv_diag.chain[uv_diag.positions[0]+2]]
                        rpair = [i for i in reversed(pair)]
                        for (pp, s) in zip([pair, rpair], ["a", "b"]):
                            uv_log_diag = copy.deepcopy(uv_diag)
                            uv_log_diag.name += "_LOG"+s
                            uv_log_diag.ct = not uv_log_diag.ct
                            # Insert momenta in the gamma chain
                            for v in pp:
                                uv_log_diag.chain.insert(
                                    uv_log_diag.positions[0], v)
                            # Update positions
                            if s == "b":
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0], uv_log_diag.positions[0]+2]
                            else:
                                # uv_log_diag.positions = [
                                #    uv_log_diag.positions[0]+1, uv_log_diag.positions[0]+2]
                                continue
                            # Rise power of the unique propagator
                            uv_log_diag.pows = [len(diag.pows)+1]

                            self.diags += [uv_log_diag]

            # Create UV Approximations

            # Number of diagrams for the approximation
            n_UV_LO_diags = 0
            for diag in diags:
                if self.uv_pos >= 0:
                    # Approximate Leading UV pole
                    uv_diag = copy.deepcopy(diag)
                    uv_diag.name += "UV_LO"
                    uv_diag.dens = [self.uv_pos+1]
                    uv_diag.pows = [len(diag.pows)]
                    self.diags += [uv_diag]
                    n_UV_LO_diags += 1

                    # Check if a bubble
                    if len(diag.dens) == 2:
                        pair = [self.uv_pos+1,
                                uv_diag.chain[uv_diag.positions[0]+2]]
                        rpair = [i for i in reversed(pair)]
                        for (pp, s) in zip([pair, rpair], ["a", "b"]):
                            uv_log_diag = copy.deepcopy(uv_diag)
                            uv_log_diag.name += "_LOG"+s
                            uv_log_diag.ct = not uv_log_diag.ct
                            # Insert momenta in the gamma chain
                            for v in pp:
                                uv_log_diag.chain.insert(
                                    uv_log_diag.positions[0], v)
                            # Update positions
                            if s == "b":
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0], uv_log_diag.positions[0]+2]
                            else:
                                uv_log_diag.positions = [
                                    uv_log_diag.positions[0]+1, uv_log_diag.positions[0]+2]
                            # Rise power of the unique propagator
                            uv_log_diag.pows = [len(diag.pows)+1]

                            self.diags += [uv_log_diag]
                            n_UV_LO_diags += 1

            # Set of diagrams to compute the amplitude and to do the same using the UV approximation
            set_amp = [
                d.name for d in self.diags if not "UV_LO" in d.name and not "born" in d.name]
            set_uv = [d.name for d in self.diags if "UV" in d.name]
            self.sets = [set_amp, set_uv]
            
            # Map denominator to LoopLines tuples
            for diag in self.diags:
                diag.dens = [self.amp_top.edge_info['p%d' % (d+1)]['loopline'] for d in diag.dens]
            
            # Print all diagrams
            #for diag in self.diags:
            #    print(diag.to_flat_format())
        else:
            print("NO SUCH AMPLITUDE! %s" % name)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""

        res = {}

        res['name'] = self.name
        res['amp_type'] = self.type
        res['topology'] = self.topology_name
        res['n_loops'] = self.n_loops
        res['diagrams'] = [diag.to_flat_format() for diag in self.diags]
        res['vectors'] = [
            [l, [float(v) for v in vec]] for (l, vec) in self.vectors]
        res['ps'] = [[float(v) for v in vec] for vec in self.ps]
        res['pols_type'] = self.polarizations
        res['sets'] = self.sets
        res['mu_r_sq'] = self.mu_r_sq

        return res



def hard_coded_dd6A(amp_name, topo_name):
    tree_factor = params['alpha_ew']**(6./2.) * params['q_d']**4 * (4.0 * np.pi)**(6./2.)
    # Initialize
    amp = Amplitude(topo_name,                   # name of topology
                    'qqbar_photons',             # process type
                    zip(["u", "vbar", "a", "a", "a", "a", "a", "a"],  # polarizations
                        ["+", "-", "+", "+", "+", "+", "+", "+"]),
                    1,                           # uv_pos
                    91.188)                      # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / (amp.qSQ(2)*amp.qSQ(3)*amp.qSQ(4)*amp.qSQ(5)*amp.qSQ(6))
    amp.add_born(generate_qqbar_born_chain(6), born_factor)
    # Diagrams and vectors
    factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
    
    vector_list=[ [[-1], sum(-amp.ps[i] for i in range(i+1))] for i in range(7) ] # add k+p slashed from prop
    vector_list.extend([[[-1], vectors.LorentzVector([0, 0, 0, 0])]])
    vector_list.extend([ [[], sum(-amp.ps[i] for i in range(1,i+1))] for i in range(2,7) ]) # add k+p when k is soft
    
    amp.create_amplitude(
       amp_name,
       1,
       generate_qqbar_diags(amp,6,factor),
       # Vectors [loopmomenta, externals]
       vector_list)
    return(amp)
    
 
def hard_coded_dd10A(amp_name, topo_name):
    tree_factor = params['alpha_ew']**(10./2.) * params['q_d']**10 * (4.0 * np.pi)**(10./2.)
    # Initialize
    amp = Amplitude(topo_name,  # name of topology
                    'qqbar_photons',             # process type
                    zip(["u", "vbar", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a"],  # polarizations
                        ["+", "-", "+", "+", "+", "+", "+", "+", "+", "+", "+", "+"]),
                    1,                           # uv_pos
                    91.188)                      # mu_r_sq
    # Born Level Diagram
    born_factor = tree_factor / (amp.qSQ(2)*amp.qSQ(3)*amp.qSQ(4)*amp.qSQ(5)*amp.qSQ(6)*amp.qSQ(7)*amp.qSQ(8)*amp.qSQ(9)*amp.qSQ(10))
    amp.add_born(generate_qqbar_born_chain(10), born_factor)
    # Diagrams and vectors
    factor = -1j*params['C_F']*tree_factor * params['alpha_s'] * (4.0 * np.pi)
    # Vectors [loopmomenta, externals]
    vector_list=[ [[-1], sum(-amp.ps[i] for i in range(i+1))] for i in range(11) ] # add k+p slashed from prop
    vector_list.extend([[[-1], vectors.LorentzVector([0, 0, 0, 0])]])
    vector_list.extend([ [[], sum(-amp.ps[i] for i in range(1,i+1))] for i in range(2,11) ]) # add k+p when k is soft
    amp.create_amplitude(
       amp_name,
       1,
       generate_qqbar_diags(amp,10,factor),
       vector_list)
    return(amp)

def generate_qqbar_born_chain(n):
    #Polarization position
    pol=(n+2)+n
    
    chain=[]
    for i in range(n):
        chain+=[pol+i]
        chain+=[i+3+n]
    del chain[-1]
    return(chain)

def generate_qqbar_diags(amp,n,factor):
    #Polarization position
    pol=(n+2)+n
    
    #Define the propagators
    max_props=[1]
    max_props.extend(range(2,n+3))
    max_den_pows=[1]*(n+2)
    max_factor = factor
    
    #Create the numerator for the largest loop
    max_chain=[-1]
    for i in range(n+1):
        max_chain+=[i+2]
        max_chain+=[pol+i]
    max_chain[-1]=-1
    
    #Define max diagram
    max_diag=(max_props,max_den_pows,max_chain,max_factor)
    diags=[max_diag]

    diags=right_left_reduce(amp,diags[-1],n)
        
    complete_diags=[]
    diag_n = 1
    for diag in diags:
        #print("D%i"%diag_n,diag[0],diag[1],diag[2],diag[3],False,len(diag[0])<4)
        complete_diags+=[qqbar_diagram("D%i"%diag_n,diag[0],diag[1],diag[2],diag[3],False,len(diag[0])<4)]
        diag_n += 1
    
    #ADD IR
    ir_chain=[-1,2]
    ir_chain.extend(generate_qqbar_born_chain(n))
    ir_chain.extend([n+2,-1])

    ir_factor=factor
    for i in range(2,n+1):
        ir_factor/=amp.qSQ(i)
    complete_diags+=[qqbar_diagram("IR",[1,2,n+2],[1,1,1],ir_chain, ir_factor,True,True)]

    
    return(complete_diags)

def right_reduce(amp,diag,n):
    chain = diag[2];  dens=diag[0]; factor=diag[3]

    #position of right -1
    pos = chain.index(-1)
    pos = chain.index(-1,pos+1)
    #Move -1 across this pair
    pair = chain[pos-2:pos]
    
    #Define new diagram
    #Change from loop propagator to externals propagator
    #Swap position of the polarization vector pair[0]
    new_chain = chain.copy()
    new_chain[pos-2]=-1
    new_chain[pos-1], new_chain[pos] = pair[1]+(n-1), pair[0]
    new_dens = dens[:-1]
    new_pows = [1]*(len(new_dens))
    new_factor = factor/amp.qSQ(pair[1]-2)
    return (new_dens,new_pows,new_chain,new_factor)

def left_reduce(amp,diag,n):
    chain = diag[2];  dens=diag[0]; factor=diag[3]

    #position of right -1
    pos = chain.index(-1)
    #Move -1 across this pair
    pair = chain[pos+1:pos+3]
    #Define new diagram
    #Change from loop propagator to externals propagator
    #Swap position of the polarization vector pair[1]
    new_chain = chain.copy()
    new_chain[pos+2]=-1
    new_chain[pos+1], new_chain[pos] = pair[0]+(n+1), pair[1]
    new_dens = dens.copy(); del new_dens[1]
    new_pows = [1]*(len(new_dens))
    new_factor = factor/amp.qSQ(pair[0])
    
    return (new_dens,new_pows,new_chain,new_factor)

def right_left_reduce(amp,diag,n):
    diags=[]
    #Do all the reductions from the left
    l_diags=[diag]
    for i in range(n-1):
        l_diags += [left_reduce(amp,l_diags[-1],n)]
    #Do all the reductions from the right
    for diag in l_diags:
        diags +=[diag]
        #Positions of the two -1
        pos1=diag[2].index(-1)
        pos2=diag[2].index(-1,pos1+1)
        while len(diags[-1][0])!=2 and (pos1+4 != 2*n-1) and (pos2-4 != 0):
            diags += [right_reduce(amp,diags[-1],n)]
            pos1=diags[-1][2].index(-1)
            pos2=diags[-1][2].index(-1,pos1+1)
    return(diags)


if __name__ == "__main__":
    hard_coded_topology_collection = TopologyCollection()
    amplitudes_collection = AmplitudesCollection()

    amp_type = 'DiHiggsTopologyLO'
    amp_name = 'gghh'
    topo_name = 'gghh_topo'
    qs =[
        vectors.LorentzVector([5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00,  5.0000000000000000e+02]),
        vectors.LorentzVector([5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.0000000000000000e+02]),
        vectors.LorentzVector([4.9999999999999977e+02,  1.0740197658523471e+02,  4.3070555989194781e+02, -1.9321629357729731e+02]),
        vectors.LorentzVector([5.0000000000000011e+02, -1.0740197658523471e+02, -4.3070555989194793e+02,  1.9321629357729719e+02]),
    ]
    ms = []

    amp_top = amplitude_topologies.create(amp_type, qs, ms, topo_name, fixed_deformation=True)
    amp_top.get_edge_info()
    #plot_edge_info(amp_top.edge_info)
    hard_coded_topology_collection.add_topology(amp_top.topology)
    # Create Amplitude ggHH
    amp = Amplitude(topo_name, amp_type, qs, ms)
    amp.create_amplitude(amp_name, 0, [], [])
    amplitudes_collection.add(amp)

   
    # add amplitude ddNA
    #amplitudes_collection.add(generate_qqbar_photons("dd4A","dd4A_topo",4))
    
    amp_type = 'qqbarphotonsNLO'
    amp_name = 'dd4A'
    topo_name = 'dd4A_topo'

    qs=[
        vectors.LorentzVector([0.5, -0.21677420292274124, 0.41425105700535525, -0.17721457817334357]),
        vectors.LorentzVector([0.0, 0.0, 0.0, 0.0]),
        vectors.LorentzVector([0.5, 0.21677420292274124, -0.41425105700535525, 0.17721457817334357]),
        vectors.LorentzVector([-0.14240531681379026, -0.11789358532625369, 0.029636836115461997, -0.07417570182224474]),
        vectors.LorentzVector([-0.07211206369149568, 0.00215409767279241, -0.048933159215160185, 0.05292499902962857]),
        vectors.LorentzVector([-0.42614282294148587, 0.34834743945157864, -0.18960367542549586, -0.15589167148122599]),
        vectors.LorentzVector([-0.3593397965532281, -0.23260795179811736, 0.20889999852519406, 0.17714237427384216]),
    ]
    ms = [1e2,0.0,0.0,0.0,0.0,0.0,0.0]
    
    amp_top = amplitude_topologies.create(amp_type, qs, ms, topo_name, fixed_deformation=True)
    amp_top.get_edge_info()
#    plot_edge_info(amp_top.edge_info)
    hard_coded_topology_collection.add_topology(amp_top.topology)
    amplitudes_collection.add(generate_qqbar_photons(qs, ms, amp_name=amp_name,topo_name=topo_name, amp_type=amp_type))
    
    # Export amplitudes
    root_dir = ''
    if "full" in sys.argv[1]:
        hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
        print("Synchronised topologies.yaml")
    amplitudes_collection.export_to(os.path.join(root_dir, 'amplitudes.yaml'))
    print("Synchronised amplitudes.yaml")
    sys.exit()
     
    amp_type = 'Nf_qqbarphotonsNNLO'
    amp_name = 'dd4A_NF'
    topo_name = 'dd4A_NF_topo'
    qs=[
        vectors.LorentzVector([0.5, -0.21677420292274124, 0.41425105700535525, -0.17721457817334357]),
        vectors.LorentzVector([0.5, 0.21677420292274124, -0.41425105700535525, 0.17721457817334357]),
        vectors.LorentzVector([-0.14240531681379026, -0.11789358532625369, 0.029636836115461997, -0.07417570182224474]),
        vectors.LorentzVector([-0.07211206369149568, 0.00215409767279241, -0.048933159215160185, 0.05292499902962857]),
        vectors.LorentzVector([-0.42614282294148587, 0.34834743945157864, -0.18960367542549586, -0.15589167148122599]),
        vectors.LorentzVector([-0.3593397965532281, -0.23260795179811736, 0.20889999852519406, 0.17714237427384216]),
    ]
    m_uv = 1e2
   
    amp_top = amplitude_topologies.create(amp_type, qs, m_uv, topo_name, fixed_deformation=False)
    amp_top.get_edge_info()
    #plot_edge_info(amp_top.edge_info, file_name='graph.pdf')
    hard_coded_topology_collection.add_topology(amp_top.topology)
    #amplitudes_collection.add(generate_qqbar_photons(qs, ms, amp_name=amp_name,topo_name=topo_name, amp_type=amp_type))

   # qs=[
   # vectors.LorentzVector([0.5, 0.1725782925959797, -0.07382815036215355, 0.46342867535227]),
   # vectors.LorentzVector([0.0, 0.0, 0.0, 0.0]),
   # vectors.LorentzVector([0.5, -0.1725782925959797, 0.07382815036215355, -0.46342867535227]),
   # vectors.LorentzVector([-0.14240531681379026, 0.05261451019605886, -0.1316847117088254, -0.013043169700204691]),
   # vectors.LorentzVector([-0.07645578602199649, -0.043069675054158645, 0.050010445882643195, 0.0385933363365972]),
   # vectors.LorentzVector([-0.034165480044262904, -0.020696767511008327, -0.013341175959563878, -0.023684105751878234]),
   # vectors.LorentzVector([-0.3745152436018937, -0.2634293072621643, 0.08221240043720716, -0.25319515986668134]),
   # vectors.LorentzVector([-0.3724581735180563, 0.2745812396312724, 0.012803041348538933, 0.25132909898216704]),
   # ]
   # ms = [1e2,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

   # amp_type = 'qqbarphotonsNLO'
   # amp_name = 'dd5A'
   # topo_name = 'dd5A_topo'

   # amp_top = amplitude_topologies.create(amp_type, qs, ms, topo_name, fixed_deformation=True)
   # amp_top.get_edge_info()
#  #  plot_edge_info(amp_top.edge_info)
   # hard_coded_topology_collection.add_topology(amp_top.topology)
   # amplitudes_collection.add(generate_qqbar_photons(qs, ms, amp_name=amp_name,topo_name=topo_name))

   # qs=[
   #     vectors.LorentzVector([0.5, -0.05113984432357776, 0.32101129713180626, 0.37991639006012284]),
   #     vectors.LorentzVector([0.0, 0.0, 0.0, 0.0]),
   #     vectors.LorentzVector([0.5, 0.05113984432357776, -0.32101129713180626, -0.37991639006012284]),
   #     vectors.LorentzVector([-0.14240531681379026, -0.13804204305188122, -0.013672853666337794, -0.03219816267300758]),
   #     vectors.LorentzVector([-0.07566099323396622, 0.04880913154311299, 0.03703210345748521, 0.044394570486394884]),
   #     vectors.LorentzVector([-0.03429822242922356, -0.010560879730620403, -0.02034249356928476, -0.025515070772064767]),
   #     vectors.LorentzVector([-0.0027507309565222744, 0.0005772669796523299, -0.0015766584134506466, 0.002178860224215476]),
   #     vectors.LorentzVector([-0.4038059382134179, 0.28245986275916485, -0.28845495918070846, -0.008330557637290945]),
   #     vectors.LorentzVector([-0.3410787983530794, -0.1832433384994286, 0.2870148613722964, 0.019470360371752933]),
   # ]
   # ms = [1e2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

   # amp_type = 'qqbarphotonsNLO'
   # amp_name = 'dd6A'
   # topo_name = 'dd6A_topo'

   # amp_top = amplitude_topologies.create(amp_type, qs, ms, topo_name, fixed_deformation=True)
   # amp_top.get_edge_info()
   # plot_edge_info(amp_top.edge_info)
   # hard_coded_topology_collection.add_topology(amp_top.topology)
   # amplitudes_collection.add(generate_qqbar_photons(qs, ms, amp_name=amp_name,topo_name=topo_name))

 
    #amplitudes_collection.add(generate_qqbar_photons("dd5A","dd5A_topo",5))
    #amplitudes_collection.add(generate_qqbar_photons("dd6A","dd6A_topo",6))
    #amplitudes_collection.add(generate_qqbar_photons("dd7A","dd7A_topo",7))
    #amplitudes_collection.add(generate_qqbar_photons("dd8A","dd8A_topo",8))
    #amplitudes_collection.add(generate_qqbar_photons("dd9A","dd9A_topo",9))
    #amplitudes_collection.add(generate_qqbar_photons("dd10A","dd10A_topo",10))
    #amplitudes_collection.add(generate_qqbar_photons("dd6A_topo","dd6A",6))
    #amplitudes_collection.add(generate_qqbar_photons("dd10A_topo","dd10A",10))

    # Export amplitudes
    root_dir = ''
    if "full" in sys.argv[1]:
        hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
        print("Synchronised topologies.yaml")
    amplitudes_collection.export_to(os.path.join(root_dir, 'amplitudes.yaml'))
    print("Synchronised amplitudes.yaml")
