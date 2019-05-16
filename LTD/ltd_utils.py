import os
import vectors
import math
import itertools
from collections import defaultdict
from pprint import pformat
import numpy as numpy
import numpy.linalg

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

#############################################################################################################
# HyperParameters
#############################################################################################################
class HyperParameters(dict):

    def __init__(self, *args, **opts):
        super(HyperParameters, self).__init__(*args, **opts)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        # For now the hyperparameters dict is already supposed to be directly exportable to yaml
        return dict(self)

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        # Directly get HyperParameters from the flat_dict
        return HyperParameters(flat_dict)

    def export_to(self, output_path, format='yaml'):
        """ Exports these hyperparameters to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export hyperparameters to yaml.")

        if output_path is not None:
            open(output_path,'w').write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
        else:
            return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)


    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import hyperparameters from yaml.")
 
        if '\n' in input_path:
            return HyperParameters.from_flat_format(yaml.load(input_path, Loader=Loader))
        else:
            return HyperParameters.from_flat_format(yaml.load(open(input_path,'r'), Loader=Loader))


#############################################################################################################
# Topology Generator
#############################################################################################################

class TopologyGenerator(object):
    def __init__(self, edge_map_lin):
        self.edge_map_lin = edge_map_lin
        self.edge_name_map = {name: i for (
            i, (name, _, _)) in enumerate(edge_map_lin)}
        self.edges = [(v1, v2) for (name, v1, v2) in edge_map_lin]
        vertices = [y for x in self.edges for y in x]

        self.num_vertices = len(set(vertices))

        self.ext = [i for i, x in enumerate(self.edges) if vertices.count(
            x[0]) == 1 or vertices.count(x[1]) == 1]
        self.loop_momenta = None
        self.propagators = None
        self.n_loops = None

    def loop_momentum_bases(self):
        trees = []
        self.spanning_trees(trees, tree={self.edges[0][0]})
        self.n_loops = len(self.edge_map_lin) - len(trees[0])
        return [[i for i in range(len(self.edge_map_lin)) if i not in tree] for tree in trees]

    def spanning_trees(self, result, tree={1}, accum=[]):
        # find all edges that connect the tree to a new node
        edges = [(i, e) for i, e in enumerate(
            self.edges) if len(tree & set(e)) == 1]

        if len(edges) == 0:
            # no more new edges, so we are done
            s = list(sorted(accum))
            if s not in result:
                result.append(s)
        else:
            for i, e in edges:
                self.spanning_trees(result, tree | set(e), accum + [i])

    def find_path(self, start, dest):
        # find all paths from source to dest
        loop = start == dest
        start_check = 1 if loop else 0
        paths = [[(start, True)]]  # store direction
        if not loop:
            paths.append([(start, False)])
        res = []
        while True:
            newpaths = []
            for p in paths:
                if len(p) > 1 and p[-1][0] == dest:
                    res.append(p)
                    continue
                last_vertex = self.edges[p[-1][0]
                                         ][1] if p[-1][1] else self.edges[p[-1][0]][0]
                for i, x in enumerate(self.edges):
                    if all(i != pp[0] for pp in p[start_check:]):
                        if loop and i == start:
                            # if we need a loop, we need to enter from the right direction
                            if x[0] == last_vertex:
                                newpaths.append(p + [(i, True)])
                            continue

                        if x[0] == last_vertex and all(x[1] not in self.edges[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, True)])
                        if x[1] == last_vertex and all(x[0] not in self.edges[pp[0]] for pp in p[start_check:-1]):
                            newpaths.append(p + [(i, False)])
            paths = newpaths
            if len(paths) == 0:
                break
        return res

    def generate_momentum_flow(self, loop_momenta):
        """Generate momentum flow where `loop_momenta` are a
        list of all edge indices that go like 1/k^2
        """
        self.loop_momenta = loop_momenta

        # first we fix the loop momentum flows
        flows = []
        for l in loop_momenta:
            paths = self.find_path(l, l)
            # make sure other loop momenta propagators are not in the path
            paths = [x for x in paths if all(
                y[0] not in self.loop_momenta for y in x[1:-1])]
            # TODO: take shortest?
            flows.append(paths[0][:-1])

        # now route the external loop_momenta to the sink
        ext_flows = []
        for e in self.ext[:-1]:
            paths = self.find_path(e, self.ext[-1])
            paths = [x for x in paths if all(
                y[0] not in self.loop_momenta for y in x[1:-1])]
            ext_flows.append(paths[0])

        # propagator momenta
        mom = []
        s = []
        for i, x in enumerate(self.edges):
            if i in self.ext:
                mom.append(())
                continue
            if i in loop_momenta:
                mom.append([(i, True)])
                s.append("1/k{}^2".format(i))
            else:
                newmom = []
                s1 = "1/("
                for j, y in enumerate(flows):
                    for yy in y:
                        if yy[0] == i:
                            s1 += ("+" if yy[1] else "-") + \
                                "k{}".format(loop_momenta[j])
                            newmom.append((loop_momenta[j], yy[1]))
                            break
                for j, y in enumerate(ext_flows):
                    for yy in y:
                        if yy[0] == i:
                            s1 += ("+" if yy[1] else "-") + \
                                "p{}".format(self.ext[j])
                            newmom.append((self.ext[j], yy[1]))
                            break
                mom.append(tuple(newmom))
                s.append(s1 + ")^2")

        self.propagators = mom
        print("Constructed topology: {}".format('*'.join(s)))

    def evaluate_residues(self, allowed_systems, close_contour):
        residues = []
        for sigmas in itertools.product(*([[1, -1]]*self.n_loops)):
            for allowed_system in allowed_systems:
                residue = self.get_thetas_for_residue(
                    sigmas, allowed_system, close_contour)
                residues += [residue]
        residues.sort()
        return residues

    def get_thetas_for_residue(self, sigmas, allowed_system, close_contour):
        # close_contour is a list, 0 means contour closed on lower half plane, 1 upper half plane for every int variable
        n_down = sum([-1 for i in close_contour if i == 0])
        # normalize sign with n_loops
        contour_sign = (-1)**n_down*(-1)**self.n_loops
        sign = numpy.prod(sigmas)*numpy.linalg.det(allowed_system[2])
        residue = [[allowed_system[0], sigmas, sign*contour_sign]]
        thetas = []
        for n_iter in range(self.n_loops):
            theta = [0 for i in range(self.n_loops)]
            sub_matrix = numpy.array(
                [[allowed_system[2][i][j] for j in range(n_iter+1)] for i in range(n_iter+1)])
            if n_iter == 0:
                theta[allowed_system[1][n_iter]] = 1. / \
                    numpy.linalg.det(sub_matrix) * \
                    sigmas[allowed_system[1][n_iter]]
                #theta[allowed_system[1][n_iter]] *= (-1.)**close_contour[n_iter]
            else:
                for r in range(n_iter+1):
                    subsub_matrix = numpy.array([[allowed_system[2][i][j] for j in range(
                        n_iter)] for i in range(n_iter+1) if i != r])
                    theta[allowed_system[1][r]] = (-1.)**(n_iter+r+1)*numpy.linalg.det(
                        subsub_matrix)/numpy.linalg.det(sub_matrix)
                    theta[allowed_system[1][r]] *= sigmas[allowed_system[1][r]]
                    #theta[allowed_system[1][r]] *= (-1.)**close_contour[n_iter]
            theta = [th*(-1.)**close_contour[n_iter] for th in theta]
            thetas += [theta]
        residue += [thetas]
        return residue

    def evaluate_thetas(self, residues):
        evaluated_theta_resides = []
        for residue in residues:
            for i, theta in enumerate(residue[1]):
                if all(x <= 0. for x in theta):
                    add = False
                    break
                else:
                    add = True
                    residue[1] = [theta for theta in residue[1]
                                  if not all(x >= 0. for x in theta)]
            if add:
                evaluated_theta_resides += [residue]
        return evaluated_theta_resides

    def remove_cancelling_residues(self, residues):
        for index_i, residue_i in enumerate(list(residues)):
            cancelling_residue = [
                [residue_i[0][0], residue_i[0][1], -residue_i[0][2]], residue_i[1]]
            if cancelling_residue in residues:
                residues.remove(cancelling_residue)
                residues.remove(residue_i)
        if len(residues) != residues[-1][0][0] + 1:
            print('Error: thetas left in residue: %s' % residues)
        return residues

    def find_allowed_systems(self, momentum_bases, signatures):
        # get transformation matrices: cut_mometum basis <=> loop momentum basis
        s_matrices = []
        for i, momentum_basis in enumerate(momentum_bases):
            s_matrix = numpy.array([signatures[i] for i in momentum_basis])
            s_matrices += [s_matrix]

        allowed_systems = []
        n_iter = 1
        for res_index, s_matrix in enumerate(s_matrices):
            permuted_s_matrices = numpy.array(
                list(itertools.permutations(s_matrix)))
            permuted_energies = numpy.array(
                list(itertools.permutations(range(self.n_loops))))
            for perm_index, permuted_s_matrix in enumerate(permuted_s_matrices):
                for n_iter in range(self.n_loops):
                    sub_matrix = numpy.array(
                        [[permuted_s_matrix[i][j] for j in range(n_iter+1)] for i in range(n_iter+1)])
                    if numpy.linalg.det(sub_matrix) == 0:
                        allowed = False
                        break
                    else:
                        allowed = True
                if allowed:
                    allowed_systems += [[res_index,
                                         permuted_energies[perm_index], permuted_s_matrix]]
                else:
                    continue
        return allowed_systems

    def get_cut_structures(self, loop_lines, contour_closure=None):
        # create the LTD graph with the external momenta removed
        graph = TopologyGenerator(
            [(None, loop_line.start_node, loop_line.end_node) for loop_line in loop_lines])
        momentum_bases = graph.loop_momentum_bases()

        signature_matrix = [loop_line.signature for loop_line in loop_lines]

        if contour_closure is None:
            # close from below by default
            contour_closure = [0] * self.n_loops

        allowed_systems = self.find_allowed_systems(
            momentum_bases, signature_matrix)
        residues = self.evaluate_residues(allowed_systems, contour_closure)
        residues = self.evaluate_thetas(residues)
        residues = self.remove_cancelling_residues(residues)

        # generate the cut structure
        cut_stucture = []
        for residue, momentum_basis in zip(residues, momentum_bases):
            cut_struct_iter = iter(residue[0][1])
            cut_stucture += [[next(cut_struct_iter)
                              if i in momentum_basis else 0 for i in range(len(loop_lines))]]
        return cut_stucture

    def create_loop_topology(self, name, ext_mom, mass_map={}, loop_momenta_names=None, contour_closure=None):
        if loop_momenta_names is None:
            loop_momenta = self.loop_momentum_bases()[0]
        else:
            self.n_loops = len(loop_momenta_names)
            loop_momenta = [self.edge_name_map[edge_name]
                            for edge_name in loop_momenta_names]

        print("Creating topology with momentum basis %s" % ', '.join([self.edge_map_lin[i][0] for i in loop_momenta]))

        self.generate_momentum_flow(loop_momenta)

        # collect all loop lines and construct the signature
        loop_line_map = defaultdict(list)
        loop_line_vertex_map = defaultdict(list)

        for prop, (edge_name, v1, v2) in zip(self.propagators, self.edge_map_lin):
            if prop == ():
                # external momentum
                continue

            mass = 0. if edge_name not in mass_map else mass_map[edge_name]

            # construct the signature
            signature = [0]*len(loop_momenta)
            q = vectors.LorentzVector([0., 0., 0., 0.])
            for (mom, sign) in prop:
                s = 1 if sign else -1
                if mom not in self.ext:
                    signature[loop_momenta.index(mom)] = s
                else:
                    q += ext_mom[self.edge_map_lin[mom][0]] * s

            # flip the sign if the inverse exists
            alt_sig = tuple(s * -1 for s in signature)
            if alt_sig in loop_line_map:
                print('warning: changing sign of propagator %s: %s -> %s' % (edge_name, tuple(signature), alt_sig) )
                loop_line_map[alt_sig].append((-q, mass))
                loop_line_vertex_map[alt_sig] += [(v2, v1)]
            else:
                loop_line_map[tuple(signature)].append((q, mass))
                loop_line_vertex_map[tuple(signature)] += [(v1, v2)]

        # fuse vertices
        for sig, vertices in loop_line_vertex_map.items():
            # find the extermal vertices
            start = next((v[0] for v in vertices if not any(v[0] == vv[1] for vv in vertices)), None)
            end =  next((v[1] for v in vertices if not any(v[1] == vv[0] for vv in vertices)), None)
            loop_line_vertex_map[sig] = (start, end)

        ll = [LoopLine(
            start_node=loop_line_vertex_map[signature][0],
            end_node=loop_line_vertex_map[signature][1],
            signature=signature,
            propagators=tuple(
                Propagator(q=q, m_squared=mass)
                for (q, mass) in propagators)) for signature, propagators in loop_line_map.items()]

        cs = self.get_cut_structures(ll, contour_closure)

        # TODO: the external kinematics are given in a random order!
        external_kinematics = list(ext_mom.values())
        external_kinematics.append(-sum(external_kinematics))
        return LoopTopology(name=name, n_loops=len(loop_momenta), external_kinematics=external_kinematics,
                            ltd_cut_structure=cs, loop_lines=ll)


#############################################################################################################
# Define topology structures
#############################################################################################################

class LoopTopology(object):
    """ A simple container for describing a loop topology."""

    def __init__(self, ltd_cut_structure, loop_lines, external_kinematics, n_loops=1, name=None, analytical_result=None, **opts):
        """
            loop_lines          : A tuple of loop lines instances corresponding to each edge of the directed
                                  graph of this topology.
            ltd_cut_structure   : A tuple of tuples instructing what cuts result from the iterative
                                  application of LTD. Example: ((0,1,-1),(1,-1,0),...)
                                  Each entry of the tuple is either 0 (NO_CUT), -1 (NEGATIVE_CUT sol.)
                                  or +1 (POSTIVIVE_CUT sol.).
            n_loops             : Number of loops in this topology
            name                : Name of this topology
        """
        self.external_kinematics = external_kinematics
        self.loop_lines          = loop_lines
        self.ltd_cut_structure   = ltd_cut_structure
        self.n_loops             = n_loops 
        self.name                = name
        if callable(analytical_result):
            self.analytical_result   = analytical_result(self.external_kinematics)
        else:
            self.analytical_result   = analytical_result

    def __str__(self):
        return pformat(self.to_flat_format())
       
    def get_com_energy(self):
        """ Returns c.o.m energy of the external kinematics."""
        incoming_momenta_sum = vectors.LorentzVector()
        for v in self.external_kinematics:
            if v[0] > 0.:
                incoming_momenta_sum += v
        return math.sqrt(abs(incoming_momenta_sum.square()))

    def print_topology(self):
        """ Print the topology using graphviz."""

        try:
            from graphviz import Digraph
        except ImportError:
            print "The print function of the LoopTopology requires the package graphviz."

        dot = Digraph(comment=self.name if self.name else 'Topology',
                        graph_attr={'fontsize':'8','sep':'10', 'splines':'true'},
                        edge_attr={'fontsize':'8'},
                        node_attr={'fontsize':'8'})

        # Define node and 
        node_properties = {'width':'0.1','height':'0.1'}
        edge_properties = {'constraint':'true', 'arrowsize':'.5'}
        
        node_added = []
        for i_ll, loop_line in enumerate(self.loop_lines):
            if loop_line.start_node not in node_added:
                node_added.append('%d'%loop_line.start_node)
                dot.node('%d'%loop_line.start_node,'%d'%loop_line.start_node, style='filled', **node_properties)
            last_node = '%d'%loop_line.start_node
            for i,propagator in enumerate(loop_line.propagators[:-1]):
                new_node = 'LL%dP%d'%(i_ll,i+1)
                node_added.append(new_node)
                dot.node(new_node,'', **node_properties)
                dot.edge(last_node, new_node, label=str(propagator.signature)+('+pi' if propagator.q!=zero_lv else ''),**edge_properties)
                last_node = new_node
            if loop_line.end_node not in node_added:
                node_added.append('%d'%loop_line.end_node)
                dot.node('%d'%loop_line.end_node,'%d'%loop_line.end_node, style='filled', **node_properties)
            dot.edge(last_node, '%d'%loop_line.end_node, label= str(loop_line.propagators[-1].signature) +
                     ('+pi' if loop_line.propagators[-1].q!=zero_lv else ''),**edge_properties)

        filename = 'topology%s'%('_%s'%self.name if self.name else '')
        dot.render(filename, view=True)

    def export_to(self, output_path, format='yaml'):
        """ Exports this topology to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Topology can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export topologies to yaml.")

        if output_path is not None:
            open(output_path,'w').write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
        else:
            return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['name'] = self.name
        res['ltd_cut_structure'] = [list(cs) for cs in self.ltd_cut_structure]
        res['n_loops'] = self.n_loops
        res['loop_lines'] = [ll.to_flat_format() for ll in self.loop_lines]
        res['external_kinematics'] = [ [float(v) for v in vec] for vec in self.external_kinematics]
        res['analytical_result_real'] = float(self.analytical_result.real) if self.analytical_result else 0.
        res['analytical_result_imag'] = float(self.analytical_result.imag) if self.analytical_result else 0. 

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return LoopTopology(
            name                =   flat_dict['name'],
            ltd_cut_structure   =   tuple([tuple(cs) for cs in flat_dict['ltd_cut_structure']]),
            n_loops             =   flat_dict['n_loops'],
            loop_lines          =   tuple([LoopLine.from_flat_format(ll) for ll in flat_dict['loop_lines']]),
            external_kinematics =   vectors.LorentzVectorList([vectors.LorentzVector(v) for v in flat_dict['external_kinematics']]),
            analytical_result   =   (None if (flat_dict['analytical_result_real']==0. and flat_dict['analytical_result_imag']==0.) 
                                     else complex(flat_dict['analytical_result_real'],flat_dict['analytical_result_imag']))
        ) 

    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Topology can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import topologies from yaml.")
 
        if '\n' in input_path:
            return LoopTopology.from_flat_format(yaml.load(input_path, Loader=Loader))
        else:
            return LoopTopology.from_flat_format(yaml.load(open(input_path,'r'), Loader=Loader))

class LoopLine(object):
    """ A simple container for describing a loop line."""

    NO_CUT          = 0
    POSITIVE_CUT    = 1
    NEGATIVE_CUT    = -1

    def __init__(self, signature, propagators, start_node, end_node, **opts):
        """
            signature           : The signature of the loop line specifies its dependence on the loop
                                  momenta. It is a tuple of length n_loops and with entries which are
                                  either 0 (no dependence), +1 (positive dependence) or -1 (negative dependence).
                                  For instance (1,0,-1) means a loop momenta dependency reading 1*k + 0*l -1*m.
            propagators         : List of Propagator instances specifying all propagators making up this loop line.
            start_node          : Integer labeling the starting node of this directed edge.
            end_node            : Integer labeling the end node of this directed edge (where the arrow points).
        """
        self.signature      = signature
        self.propagators    = propagators
        # Automatically forward the signature attribute to the loop propagators if not specified by the user.
        for propagator in self.propagators:
            propagator.set_signature(self.signature, force=False)
        self.start_node     = start_node
        self.end_node       = end_node

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return LoopLine(
            signature   =   tuple(flat_dict['signature']),
            propagators =   tuple([Propagator.from_flat_format(p) for p in flat_dict['propagators']]),
            start_node  =   flat_dict['start_node'],
            end_node    =   flat_dict['end_node']             
        ) 


    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['signature'] = list(self.signature)
        res['propagators'] = [p.to_flat_format() for p in self.propagators]
        res['start_node'] = self.start_node
        res['end_node'] = self.end_node

        return res

class Propagator(object):
    """ A simple container for describing a loop propagator."""

    def __init__(self, q, m_squared, signature = None, **opts):
        self.q          = q
        self.m_squared  = m_squared
        # Note that this signature member is not strictly speaking necessary as it should always be the
        # same as the signature attribute of the LoopLine containing it. We forward it here for convenience however.
        self.signature  = signature

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['q'] = [float(v) for v in self.q]
        res['m_squared'] = self.m_squared

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return Propagator(
            q           =   vectors.LorentzVector(flat_dict['q']),
            m_squared   =   flat_dict['m_squared'] 
        ) 

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse propagator with the provided list loop momenta, given as a list of LorentzVector.""" 
        return (sum(wgt*loop_momentum for wgt, loop_momentum in zip(self.signature, loop_momenta)) + self.q).square() - self.m_squared

    def set_signature(self, signature, force=True):
        """ Overwrite (if force=True) the signature attribute of this propagator."""
        if (self.signature is None) or force:
            self.signature = signature

class TopologyCollection(dict):
    """ Collection of hard-coded topologies."""

    def add_topology(self, topology_to_add, entry_name=None):

        assert(isinstance(topology_to_add,LoopTopology))
        if entry_name is None:
            if topology_to_add.name is None:
                raise BaseException("Specify the name of the hard-coded topology to add to the collection.")
            self[topology_to_add.name] = topology_to_add
        else:
            self[entry_name] = topology_to_add

    def export_to(self, output_path, format='yaml'):
        """ Exports this topology to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Topology can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export topologies to yaml.")

        flat_record = [self[k].to_flat_format() for k in sorted(self.keys())]
        if output_path is not None:
            open(output_path,'w').write(yaml.dump(flat_record, Dumper=Dumper))
        else:
            return yaml.dump(flat_record, Dumper=Dumper)

    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Topology can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import topologies from yaml.")
        
        if '\n' in input_path:
            flat_record = yaml.load(input_path, Loader=Loader)
        else:
            flat_record = yaml.load(open(input_path,'r'), Loader=Loader)

        result = TopologyCollection()
        for topology in flat_record:
            result.add_topology(LoopTopology.from_flat_format(topology))
        
        return result

if __name__ == "__main__":
    triangle = TopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3),
                                ('p3', 2, 3), ('q2', 3, 4), ('q3', 2, 5)])  # triangle
    doubletriangle = TopologyGenerator([('q', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 3),
                                        ('p4', 3, 4), ('p5', 2, 4), ('-q', 4, 5)])  # double-triangle

    loop_topology = doubletriangle.create_loop_topology("DoubleTriangle", ext_mom={'q': vectors.LorentzVector(
        [0.1, 0.2, 0.3, 0.4])}, mass_map={'p1': 1.0, 'p2': 2.0, 'p3': 3.0}, loop_momenta_names=('p1', 'p5'))

    print(loop_topology)