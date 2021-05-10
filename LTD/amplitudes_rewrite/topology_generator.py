import numpy
import os
import sys

if True:
    src_path = os.path.dirname(os.path.realpath(__file__))    
    sys.path.append(src_path)
    sys.path.append(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..'))
    from ltd_utils import *
    from amplitudes_rewrite.amplitudes import Amplitude


class SquaredTopologyGeneratorForAmplitudes(TopologyGeneratorFromPropagators):
    def __init__(self, incoming_momenta, outgoing_momenta, masses, propagators, name, spinors={}, polarsations={}, color_struc='color(1)',extra_calls=[]):
        self.name = name
        self.incoming_momenta = incoming_momenta
        self.outgoing_momenta = outgoing_momenta
        self.masses = masses
        self.propagators = propagators
        self.n_loops_subgraph = len(self.propagators[0]["loop_signature"])
        assert(all([len(propagator["loop_signature"]) ==
                    self.n_loops_subgraph for propagator in self.propagators]))
        self.n_outgoing = len(self.outgoing_momenta)
        self.n_incoming = len(self.incoming_momenta)
        # don't change, this is an member of TopologyGeneratorFromPropagators that will be used later
        self.n_loops = None
        self.spinor_v = spinors.get("v", [])
        self.spinor_u = spinors.get("u", [])
        self.spinor_ubar = spinors.get("ubar", [])
        self.spinor_vbar = spinors.get("vbar", [])
        self.cpol = polarsations.get("cpol", [])
        self.pol = polarsations.get("pol", [])
        self.color_struc = color_struc
        self.extra_calls = extra_calls

    @classmethod
    def from_amplitude(cls, amplitude):
        assert(isinstance(amplitude, Amplitude))
        name = amplitude.name
        # amplitude.get('name')
        external_data = amplitude.external_data
        incoming_momenta = [['p%i' % (i+1), numpy.array(p)]
                            for i, p in enumerate(external_data["in_momenta"])]
        outgoing_momenta = [['q%i' % (i+1), numpy.array(q)]
                            for i, q in enumerate(external_data["out_momenta"])]
        spinors = {}
        spinors["v"] = external_data.get("spinor_v", [])
        spinors["u"] = external_data.get("spinor_u", [])
        spinors["vbar"] = external_data.get("spinor_vbar", [])
        spinors["ubar"] = external_data.get("spinor_ubar", [])
        polarsations = {}
        polarsations["cpol"] = external_data.get("cpol", [])
        polarsations["pol"] = external_data.get("pol", [])
        color_struc = amplitude.color_strucs[0]
        if external_data["n_out"] == len(external_data["out_momenta"])+1:
            q_sum = numpy.sum([p[1] for p in incoming_momenta]+[-q[1]
                                                                for q in outgoing_momenta], axis=0)
            outgoing_momenta += [['q%i' % external_data["n_out"], q_sum]]
        diagram_list = amplitude.diagram_list

        all_propagators_dict = {}
        for diagram in diagram_list:
            diag_props = diagram["propagators"]
            prop_id = 1
            for prop in diag_props:
                # assign a unique id to the same props
                combined_sig = prop["loop_signature"]+prop["outgoing_signature"] + prop["incoming_signature"]
                unique_id = {tuple(combined_sig),tuple(-s for s in combined_sig),prop["mass"]}

                # do the merging
                all_propagators_dict.update({tuple(unique_id): (
                    # enforce positive signatures
                    [abs(s) if sum(abs(sig) for sig in combined_sig) == 1 else s for s in prop["loop_signature"]],
                    prop["outgoing_signature"],
                    prop["incoming_signature"],
                    prop["mass"]
                )})

        propagators = [
            {
                "loop_signature": prop[0],
                "outgoing_signature": prop[1],
                "incoming_signature": prop[2],
                "mass": prop[3],
                "power": 1,  # irrelevant for detection of E-surfaces
                "name": 'l%i' % (i+1)
            }
        for i, prop in enumerate(all_propagators_dict.values()) ]
            
        masses = amplitude.masses
        
        additional_options = amplitude.additional_options
        if additional_options["integrand_per_diag"]:
            extra_calls = [[i+1,0] for i in range(len(diagram_list)-1)]
        else:
            extra_calls = []
        return SquaredTopologyGeneratorForAmplitudes(incoming_momenta, outgoing_momenta, masses, propagators, name, spinors=spinors, polarsations=polarsations, color_struc=color_struc, extra_calls=extra_calls)

    def get_supergraph_topology(self):
        supergraph_loop_lines = self.get_supergraph_loop_lines()
        # for ll in self.supergraph_loop_lines:
        #	print(ll.to_flat_format())
        assert(self.n_loops == None)
        # required to call get_cut_structure
        self.n_loops = len(supergraph_loop_lines[0].signature)
        #supergraph_cut_structure = self.get_cut_structures(
        #    supergraph_loop_lines)
        pure_signatures = [[int(value) for value in line] for line in numpy.identity(self.n_loops)]
        dummy_supergraph_cut_structure = [[1 if ll.signature in pure_signatures else 0 for ll in supergraph_loop_lines]]
        assert(sum(dummy_supergraph_cut_structure[0])==self.n_loops)
        #print(dummy_supergraph_cut_structure)
        supergraph_topology = LoopTopology(
            dummy_supergraph_cut_structure,
            supergraph_loop_lines,
            2*[in_mom[1] for in_mom in self.incoming_momenta],
            n_loops=self.n_loops,
            name=self.name,
            constant_deformation={"alpha": [[0.0, 0.0, 0.0, 1.0]], "beta": [
                [0.0, 0.0, 0.0, 0.0]]},  # is default
        )
        self.n_loops = None
        return supergraph_topology

    def get_supergraph_loop_lines(self):
        supergraph_loop_lines = []
        all_signatures = [tuple(propagator["loop_signature"]+propagator["outgoing_signature"])
                          for propagator in self.propagators]
        # maybe have to remove the ones where -1 * signature1 = signature2
        unique_signatures = set(all_signatures)
        for signature in unique_signatures:
            ll_propagators = []
            for propagator in self.propagators:
                if tuple(propagator["loop_signature"]+propagator["outgoing_signature"]) == signature:
                    ll_propagators += [Propagator(
                        numpy.sum([sign*in_mom[1] for sign, in_mom in zip(
                            propagator["incoming_signature"], self.incoming_momenta)], axis=0),
                        self.masses[propagator["mass"]]**2 if isinstance(
                            propagator["mass"], str) else propagator["mass"]**2,  # needed
                        power=propagator["power"],  # needed
                        name=propagator["name"],  # needed
                    )]
            supergraph_loop_lines += [LoopLine(list(signature),
                                               ll_propagators, 0, 0)]
        outgoing_loop_lines = []
        for i in range(self.n_outgoing):
            if i != self.n_outgoing-1:
                dummy_signature = [0]*self.n_loops_subgraph + \
                    [0 if i != j else 1 for j in range(self.n_outgoing-1)]
                dummy_q = numpy.array([0, 0, 0, 0])
            else:
                dummy_signature = [0]*self.n_loops_subgraph + \
                    [-1]*(self.n_outgoing-1)
                dummy_q = sum(in_mom[1] for in_mom in  self.incoming_momenta)
            dummy_propagators = [Propagator(
                dummy_q,
                # mass of cut momenta not needed as it is automatically
                # determined from fixed values #float(self.outgoing_momenta[i][1].square()),
                0.,
                name=self.outgoing_momenta[i][0])]
            outgoing_loop_lines += [LoopLine(dummy_signature,
                                             dummy_propagators, 0, 0)]
        supergraph_loop_lines += outgoing_loop_lines
        return supergraph_loop_lines

    def get_left_graph_loop_lines(self):
        left_graph_loop_lines = []
        all_signatures = [tuple(propagator["loop_signature"])
                          for propagator in self.propagators]
        unique_signatures = set(all_signatures)
        for signature in unique_signatures:
            ll_propagators = []
            for propagator in self.propagators:
                if tuple(propagator["loop_signature"]) == signature:
                    ll_propagators += [Propagator(
                        # q is always zero because parametric_shift is supplied
                        numpy.array([0., 0., 0., 0.]),
                        self.masses[propagator["mass"]]**2 if isinstance(
                            propagator["mass"], str) else propagator["mass"]**2,  # needed
                        # parametric_shift signature is in cut basis, i.e. [out1,..,outm-1,0] (last entry is the dependent cut momentum)
                        # parametric_shift incomming signature in terms of "in1,..inl,0,..,0" (last entries are l-times 0)
                        parametric_shift=[propagator["outgoing_signature"]+[0],
                                          propagator["incoming_signature"]+[0]*self.n_incoming],
                        power=propagator["power"],  # needed
                        name=propagator["name"],  # needed
                    )]
            left_graph_loop_lines += [LoopLine(list(signature),
                                               ll_propagators, 0, 0)]
        return left_graph_loop_lines

    def get_left_graph_topology(self):
        left_graph_loop_lines = self.get_left_graph_loop_lines()
        assert(self.n_loops == None)
        self.n_loops = self.n_loops_subgraph
        left_graph_cut_structure = self.get_cut_structures(
            left_graph_loop_lines)
        loop_momentum_map_signatures = [[1 if i == j else 0 for i in range(
            self.n_loops_subgraph)]+[0]*(self.n_outgoing-1) for j in range(self.n_loops_subgraph)]
        loop_momentum_map = [[signature, [
            0]*(2*self.n_incoming)] for signature in loop_momentum_map_signatures]
        left_graph_topology = LoopTopology(
            left_graph_cut_structure,
            left_graph_loop_lines,
            [],  # do not provide external kinematics
            n_loops=self.n_loops,
            name=self.name+"_0",
            constant_deformation={"alpha": [[0.0, 0.0, 0.0, 1.0]], "beta": [
                [0.0, 0.0, 0.0, 0.0]]},  # is default
            loop_momentum_map=loop_momentum_map
        )
        self.n_loops = None
        return left_graph_topology

    def get_cuts(self):
        outgoing_momenta_signatures = []
        for i in range(len(self.outgoing_momenta)):
            # loop_signature in terms of "k1,..kn, out1, .. outm-1"
            # shift signature in terms of "in1,..inl,0,..,0" (l-times 0)
            if i != len(self.outgoing_momenta)-1:
                loop_signature = [0]*self.n_loops_subgraph + \
                    [0 if i != j else 1 for j in range(self.n_outgoing-1)]
                shift_signature = [0]*(2*self.n_incoming)
            else:
                # defined by momentum conservation
                loop_signature = [0]*self.n_loops_subgraph + \
                    [-1]*(self.n_outgoing-1)
                shift_signature = [1]*self.n_incoming+[0]*self.n_incoming
            outgoing_momenta_signatures += [[loop_signature, shift_signature]]
        cuts = [{
                "m_squared": 0.,  # not used
                "name": outgoing_momentum[0],
                "particle_id": 666,  # not used
                "power": 1,  # not used
                "sign": 1,  # not used
                "signature": [signature[0], signature[1]]
                } for i, (outgoing_momentum, signature) in enumerate(zip(self.outgoing_momenta, outgoing_momenta_signatures))]
        return cuts

    def get_diagram_sets(self, cuts):
        left_graph_topology = self.get_left_graph_topology()
        right_graph_topology = LoopTopology([[]], [], [],
                                            n_loops=0,
                                            name=self.name+"_1",
                                            constant_deformation={"alpha": [[0.0, 0.0, 0.0, 1.0]], "beta": [
                                                [0.0, 0.0, 0.0, 0.0]]},  # is default
                                            )
        left_diag = {"conjugate_deformation": False,
                     "graph": left_graph_topology.to_flat_format()}
        right_diag = {"conjugate_deformation": True,
                      "graph": right_graph_topology.to_flat_format()}
        # cm_to_lmb is the inverse of the matrix M of the trsf
        # [out1,.. out(m-1), k1, .. kn] = M [k1, .. , kn,out1,..out(m-1)]
        # M = [loop_signature_cut1, .., loop_signature_cutm-1,[1,...,0,0...,0], ... [0,...,1,0,...,0] #n unit matrix and m-1 zeros at the end #]
        matrix = numpy.array([cut["signature"][0] for cut in cuts[:-1]] +
                             [[1 if i == j else 0 for i in range(self.n_loops_subgraph)]+[0]*(self.n_outgoing-1) for j in range(self.n_loops_subgraph)])
        cb_to_lmb = numpy.linalg.inv(matrix)
        diagram_sets = [{"cb_to_lmb": [int(i) for i in cb_to_lmb.flatten()],
                         "diagram_info": [left_diag, right_diag],
                         "id":0,
                         "numerator_tensor_coefficients_sparse":[],
                         "uv_spinney": []
                         }]
        return diagram_sets

    def get_cross_section_set(self):
        cross_section_set = {}
        cross_section_set["color_struc"] = self.color_struc
        cross_section_set["external_data"] = {"cpol": self.cpol,
                                              "in_momenta": [[float(v) for v in p[1]] for p in self.incoming_momenta],
                                              "out_momenta": [[float(v) for v in p[1]] for p in self.outgoing_momenta[:-1]],
                                              "n_in": self.n_incoming,
                                              "n_out": self.n_outgoing,
                                              "pol": self.pol,
                                              "spinor_u": self.spinor_u,
                                              "spinor_ubar": self.spinor_ubar,
                                              "spinor_v": self.spinor_v,
                                              "spinor_vbar": self.spinor_vbar,
                                              }
        cross_section_set["name"] = self.name
        cross_section_set["topologies"] = [
            {"additional_LMBs": [], "multiplicity": 1, "name": self.name+'_0'}]
        return cross_section_set

    def export_yaml(self, output_path):
        out = {}
        out["topo"] = self.get_supergraph_topology().to_flat_format()
        cuts = self.get_cuts()
        diagram_sets = self.get_diagram_sets(cuts)
        out["cutkosky_cuts"] = [{"cuts": cuts, "diagram_sets": diagram_sets}]
        out["FORM_integrand"] = {
            "call_signature": {"extra_calls": self.extra_calls, "id": 0}}
        out["FORM_numerator"] = {"call_signature": {"id": 0}}
        out["MG_numerator"] = {}
        out["color_struc"] = self.color_struc
        out["default_fixed_cut_momenta"] = [[], []]
        # [[[float(v) for v in q[1]] for q in self.incoming_momenta], #incoming
        # [[float(v) for v in q[1]] for q in self.outgoing_momenta]] #outgoing
        tot_in = numpy.sum([p[1] for p in self.incoming_momenta], axis=0)
        out["e_cm_squared"] = float(tot_in[0]**2-tot_in[1:].dot(tot_in[1:]))
        out["edge_PDGs"] = []  # not necessary
        out["edge_signatures"] = {}  # not necessary
        out["external_momenta"] = [
            [float(v) for v in p[1]] for p in 2*self.incoming_momenta]
        # basis_signatures are unit matrix
        basis_signatures = [[1 if i == j else 0 for i in range(
            self.n_loops_subgraph)] for j in range(self.n_loops_subgraph)]
        # find the propagators that have the basis siganature and no shift
        basis_loop_momentum_name = [[[prop["name"]
                                      for prop in ll["propagators"]
                                      #with [0] select first element (in case two props of the same loop line have no shifts)
                                      if all([value == 0 for value in prop["parametric_shift"][0]+prop["parametric_shift"][1]])][0]
                                     for ll in diagram_sets[0]["diagram_info"][0]["graph"]["loop_lines"]
                                     if ll["signature"] == basis_signature]
                                    for basis_signature in basis_signatures]                
        basis_loop_momentum_name = [str(name) for name in numpy.array(
            basis_loop_momentum_name).flatten()]
        assert(len(basis_loop_momentum_name) == self.n_loops_subgraph)
        out["loop_momentum_basis"] = basis_loop_momentum_name+[cut["name"]
                                                               for cut in cuts[:-1]]  # has to be of form k1...kn,out1,..outm-1
        out["n_incoming_momenta"] = self.n_incoming
        out["n_loops"] = self.n_loops_subgraph + self.n_outgoing - 1
        out["name"] = self.name
        out["overall_numerator"] = 1.0
        out["topo_edges"] = []  # not necessary

        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException(
                "Install yaml python module in order to import topologies from yaml.")
        open(output_path+'/'+self.name+'.yaml', 'w').write(yaml.dump(
            self.get_cross_section_set(), Dumper=Dumper, default_flow_style=None))
        open(output_path+'/'+self.name+'_0.yaml',
             'w').write(yaml.dump(out, Dumper=Dumper, default_flow_style=None))


if __name__ == "__main__":
    # example explicit
    p1 = numpy.array([1, 0, 0, 1])
    p2 = numpy.array([-1, 0, 0, 1])
    q1 = numpy.array([0.0, 2.3, 0.2, 1.1])
    q2 = numpy.array([0.0, 0.2, 1.3, 3.123])
    q3 = p1+p2-q1-q2
    incoming_momenta = [["p1", p1],
                        ["p2", p2]]
    outgoing_momenta = [["q1", q1],
                        ["q2", q2],
                        ["q3", q3]]
    masses = {"m1": 1, "m2": 1, "m3": 1, "m4": 1, "m5": 1}
    propagators = []
    propagators += [{"loop_signature": [1],
                     "outgoing_signature": [0, 0],
                     "incoming_signature": [1, 0],
                     "mass": "m1",
                     "power": 1,
                     "name": "l1"}]
    propagators += [{"loop_signature": [1],
                     "outgoing_signature": [0, 0],
                     "incoming_signature": [1, 1],
                     "mass": "m2",
                     "power": 1,
                     "name": "l2"}]
    propagators += [{"loop_signature": [1],
                     "outgoing_signature": [-1, 0],
                     "incoming_signature": [1, 1],
                     "mass": "m3",
                     "power": 1,
                     "name": "l3"}]
    propagators += [{"loop_signature": [1],
                     "outgoing_signature": [-1, -1],
                     "incoming_signature": [1, 1],
                     "mass": "m4",
                     "power": 1,
                     "name": "l4"}]
    propagators += [{"loop_signature": [1],
                     "outgoing_signature": [0, 0],
                     "incoming_signature": [0, 0],
                     "mass": "m5",
                     "power": 1,
                     "name": "l5"}]
    name = "test"
    topo = SquaredTopologyGeneratorForAmplitudes(
        incoming_momenta, outgoing_momenta, masses, propagators, name)
    topo.export("./"+topo.name+".yaml")

    # example from runcard
    runcard_path = "/Users/Dario/Desktop/PhD/Code/alphaloop/LTD/amplitudes_rewrite/runcard_template_scalar.yaml"
    topo = SquaredTopologyGeneratorForAmplitudes.from_runcard(runcard_path)
    topo.export("./"+topo.name+".yaml")
