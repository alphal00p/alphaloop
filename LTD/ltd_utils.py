import os
import sys

import os
import sys
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)

import vectors
import pprint
import copy
import math
import itertools
from collections import defaultdict
from pprint import pformat
import numpy as numpy
import numpy.linalg
from itertools import product, combinations, permutations
import multiprocessing
import signal
import time
from sympy import Matrix, S, nsimplify, eye

try:
    import cvxpy
except:
    print("Error: could not import package cvxpy necessary for building the fixed deformation. Make sure it is installed.")

zero_lv = vectors.LorentzVector([0.,0.,0.,0.])

class Colour:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

#===============================================================================
# Low-level muter that works with f2py as well
# From: http://code.activestate.com/recipes/577564/
#===============================================================================
class Silence:
    """Context manager which uses low-level file descriptors to suppress
    output to stdout/stderr, optionally redirecting to the named file(s)."""

    def __init__(self, stdout=os.devnull, stderr=os.devnull, mode='w', active=True):
        self.active = active
        if not self.active: return
        self.outfiles = stdout, stderr
        self.combine = (stdout == stderr)
        self.mode = mode
        
    def __enter__(self):
        if not self.active: return
        self.sys = sys
        # save previous stdout/stderr
        self.saved_streams = saved_streams = sys.__stdout__, sys.__stderr__
        self.fds = fds = [s.fileno() for s in saved_streams]
        self.saved_fds = map(os.dup, fds)
        # flush any pending output
        for s in saved_streams: s.flush()

        # open surrogate files
        if self.combine: 
            null_streams = [open(self.outfiles[0], self.mode)] * 2
            if self.outfiles[0] != os.devnull:
                # disable buffering so output is merged immediately
                sys.stdout, sys.stderr = map(os.fdopen, fds, ['w']*2, [0]*2)
        else:
            null_streams = [open(f, self.mode) for f in self.outfiles]
        self.null_fds = null_fds = [s.fileno() for s in null_streams]
        self.null_streams = null_streams
        
        # overwrite file objects and low-level file descriptors
        map(os.dup2, null_fds, fds)

    def __exit__(self, *args):
        if not self.active: return
        sys = self.sys
        # flush any pending output
        for s in self.saved_streams: s.flush()
        # restore original streams and file descriptors
        map(os.dup2, self.saved_fds, self.fds)
        sys.stdout, sys.stderr = self.saved_streams
        # clean up
        for s in self.null_streams: s.close()
        for fd in self.saved_fds: os.close(fd)
        return False

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

    def __init__(self, edge_map_lin, powers=None):
        #if len(set(e[0] for e  in edge_map_lin)) != len(edge_map_lin):
        #    raise AssertionError("Every edge must have a unique name. Input: ", edge_map_lin)
        self.edge_map_lin = edge_map_lin
        self.edge_name_map = {name: i for (
            i, (name, _, _)) in enumerate(edge_map_lin)}
        self.edges = [(v1, v2) for (name, v1, v2) in edge_map_lin]
        vertices = [y for x in self.edges for y in x]

        self.num_vertices = len(set(vertices))

        self.ext = [i for i, x in enumerate(self.edges) if vertices.count(
            x[0]) == 1 or vertices.count(x[1]) == 1]
        self.vertices = vertices # with repetition
        # set the power for every propagator
        self.powers = {e[0]: powers[e[0]] if powers is not None and e[0] in powers
                else 1 for i, e in enumerate(edge_map_lin) if i not in self.ext}
        self.loop_momenta = None
        self.propagators = None
        self.n_loops = None
    
    def loop_momentum_bases(self):
        trees = []
        self.spanning_trees(trees, tree={self.edges[0][0]})
        self.n_loops = len(self.edge_map_lin) - len(trees[0])
        return [[i for i in range(len(self.edge_map_lin)) if i not in tree] for tree in trees]

    def merge_duplicate_edges(self, loop_momenta_names=None):
        """Create a new topology, merging edges with the same or opposite signature.
        If it is opposite, one of the lines will be flipped and the power will be increased.
        """
        if self.propagators is None:
            self.generate_momentum_flow(loop_momenta_names)

        duplicates = {}
        new_powers = {}
        if loop_momenta_names:
            duplicates = {tuple(self.propagators[self.edge_name_map[l]]): l for l in loop_momenta_names}
            new_powers = {e: p for e, p in self.powers.items() if e in loop_momenta_names}

        fuse_verts = set()
        new_edge_map_lin = []
        for p, e in zip(self.propagators, self.edge_map_lin):
            min_p = tuple([(x[0], not x[1]) for x in p])

            if p != () and min_p in duplicates:
                print('Flipping edge {} in deduplication process {} -> {}'.format(e[0], e[0], duplicates[min_p]))
                new_powers[duplicates[min_p]] += 1
                fuse_verts.add((e[1], e[2]))
            elif p != () and tuple(p) in duplicates and duplicates[tuple(p)] != e[0]:
                # an edge is a duplicate if it is internal with repeated signatures and it is not one of the loop momenta
                new_powers[duplicates[tuple(p)]] += 1
                fuse_verts.add((e[1], e[2]))
            else:
                new_edge_map_lin.append(e)
                if p != () and tuple(p) not in duplicates:
                    duplicates[tuple(p)] = e[0]
                    new_powers[e[0]] = self.powers[e[0]]

        for v1, v2 in fuse_verts:
            new_edge_map_lin = [(name, o1 if o1 != v2 else v1, o2 if o2 != v2 else v1) for name, o1, o2 in new_edge_map_lin]

        return TopologyGenerator(new_edge_map_lin, powers=new_powers)

    def spanning_trees(self, result, tree={1}, accum=[]):
        # find all edges that connect the tree to a new node
        edges = [(i, e) for i, e in enumerate(
            self.edges) if e[0] != e[1] and len(tree & set(e)) == 1]

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

    def find_cutkosky_cuts(self, n_jets, incoming_particles, final_state_particle_ids, particle_ids):
        spanning_trees = []
        self.spanning_trees(spanning_trees)
        cutkosky_cuts = set()
        for spanning_tree in spanning_trees:
            # now select the extra cut
            for edge_index in spanning_tree:
                if edge_index in self.ext:
                    continue

                # verify that the graph is correctly split in two, with one of the subgraphs having all incoming particles
                # and the other having all outgoing particles
                cut_tree = TopologyGenerator([e for i, e in enumerate(self.edge_map_lin) if i in set(spanning_tree) - {edge_index,}])
                sub_tree_indices = []
                cut_tree.spanning_trees(sub_tree_indices)
                sub_tree = TopologyGenerator([cut_tree.edge_map_lin[i] for i in sub_tree_indices[0]])
                is_incoming_sub_tree = any(e[0] in incoming_particles for e in sub_tree.edge_map_lin)

                if set(set([sub_tree.edge_map_lin[i][0] for i in sub_tree.ext]) & set([self.edge_map_lin[i][0] for i in self.ext])) == set(incoming_particles):
                    # identify which of the n+1 cuts are the cuts that separate the original diagram in two
                    cutkosky_cut = tuple(sorted(cutkosky_edge for cutkosky_edge in set(self.edge_map_lin).difference(set(cut_tree.edge_map_lin))
                            if len(set(sub_tree.vertices) & set([cutkosky_edge[1], cutkosky_edge[2]]))==1))

                    # filter cuts with not enough jets and cuts that do not contain all desired final state particles
                    cutkosky_particles = tuple(sorted(particle_ids[e] for (e, _, _) in cutkosky_cut if e in particle_ids))
                    if len(cutkosky_cut) - len(final_state_particle_ids) < n_jets or \
                            not all(final_state_particle_ids.count(p) <= cutkosky_particles.count(p) for p in final_state_particle_ids):
                        continue

                    # check if the cut is incoming our outgoing
                    cutkosky_cut = tuple((name, 1 if is_incoming_sub_tree and sub_tree.vertices.count(v2) == 0 or
                                        not is_incoming_sub_tree and sub_tree.vertices.count(v1) == 0 else -1)
                                    for (name, v1, v2) in cutkosky_cut)
                    cutkosky_cuts.add(cutkosky_cut)
        return list(sorted(cutkosky_cuts))

    def split_graph(self, cutkosky_cut, incoming_momenta):
        """Split the graph in two using a cutkosky cut. The first has the incoming momenta, and the second the outgoing."""
        cut_tree = TopologyGenerator([e for e in self.edge_map_lin if e[0] not in cutkosky_cut])
        sub_tree_indices = []
        cut_tree.spanning_trees(sub_tree_indices)

        # get the vertices from the spanning tree and use it to create the two subgraphs
        sub_tree_vertices = set(v for e in sub_tree_indices[0] for v in cut_tree.edge_map_lin[e][1:])
        left_tree = [e for e in cut_tree.edge_map_lin if e[1] in sub_tree_vertices or e[2] in sub_tree_vertices]
        right_tree = [e for e in cut_tree.edge_map_lin if e not in left_tree]

        left_propagators = set(e[0] for e in left_tree)
        right_propagators = set(e[0] for e in right_tree)

        # add the cutkosky cut edges with renamed external vertices to prevent accidental loops
        highest_vertex = max(v for e in self.edge_map_lin for v in e[1:]) + 1
        for t in (left_tree, right_tree):
            for c in cutkosky_cut:
                (_, cut_v1, cut_v2) = next(e for e in self.edge_map_lin if e[0] == c)
                if any(e[1] == cut_v1 or e[2] == cut_v1 for e in t):
                    assert(not any (e[1] == cut_v2 or e[2] == cut_v2 for e in t))
                    t.append((c, cut_v1, highest_vertex))
                    highest_vertex += 1
                else:
                    assert(not any (e[1] == cut_v1 or e[2] == cut_v1 for e in t))
                    t.append((c, highest_vertex, cut_v2))
                    highest_vertex += 1

        l = TopologyGenerator(left_tree, powers={e: p for e, p in self.powers.items() if e in left_propagators})
        r = TopologyGenerator(right_tree, powers={e: p for e, p in self.powers.items() if e in right_propagators})

        if any(e[0] in incoming_momenta for e in left_tree):
            return (l, r)
        else:
            return (r, l)

    def get_signature_map(self):
        # TODO: the order of self.ext should be consistent with Rust
        edge_map = {}
        for i, (prop, (edge_name, _, _)) in enumerate(zip(self.propagators, self.edge_map_lin)):
            signature = [[0]*len(self.loop_momenta), [0]*len(self.ext)]

            if prop == ():
                signature[1][self.ext.index(i)] = 1 # FIXME: is it always +?

            for (mom, sign) in prop:
                s = 1 if sign else -1
                if mom not in self.ext:
                    signature[0][self.loop_momenta.index(mom)] = s
                else:
                    signature[1][self.ext.index(mom)] = s
            edge_map[edge_name] = signature
        return edge_map

    def build_proto_topology(self, sub_graph, cuts, loop_momenta_names=None):
        sub_graph.generate_momentum_flow()

        # for each edge, determine the momentum map from full graph to subgraph and the shift
        # the shift will in general also have loop momentum dependence
        # TODO: the order of self.ext should be consistent with Rust
        edge_map = self.get_signature_map()
        sub_graph_edge_map = sub_graph.get_signature_map()
        loop_momentum_map = [[]]*sub_graph.n_loops
        param_shift = {}
        for (prop, (edge_name, _, _)) in zip(sub_graph.propagators, sub_graph.edge_map_lin):
            if prop == ():
                continue

            # determine the loop momentum map from the full loop momentum basis to the one of the subgraph
            if len(prop) == 1 and prop[0][0] not in sub_graph.ext:
                (mom, sign) = prop[0]
                s = 1 if sign else -1
                loop_momentum_map[sub_graph.loop_momenta.index(mom)] = [[s * a for a in i] for i in edge_map[edge_name]]

            # determine the shift in the subgraph in the cut basis
            signature = [[0]*len(cuts), [0]*len(self.ext)]

            # map the external momenta back to cuts and the external momenta of the full graph
            for ext_index, s in enumerate(sub_graph_edge_map[edge_name][1]):
                mom = sub_graph.edge_map_lin[sub_graph.ext[ext_index]][0]
                if s != 0:
                    if mom in cuts:
                        signature[0][cuts.index(mom)] = s
                    else:
                        edge_index = next(i for i, e in enumerate(self.edge_map_lin) if e[0] == mom)
                        signature[1][self.ext.index(edge_index)] = s
            param_shift[edge_name] = signature

        return loop_momentum_map, param_shift

    def generate_momentum_flow(self, loop_momenta=None):
        if loop_momenta is None:
            loop_momenta = self.loop_momentum_bases()[0]
        else:
            self.n_loops = len(loop_momenta)
            if isinstance(loop_momenta[0], str):
                loop_momenta = [self.edge_name_map[edge_name]
                                for edge_name in loop_momenta]
            else:
                loop_momenta = loop_momenta

        self.loop_momenta = loop_momenta

        flows = []
        for l in loop_momenta:
            paths = self.find_path(l, l)
            # make sure other loop momenta propagators are not in the path
            paths = [x for x in paths if all(
                y[0] not in loop_momenta for y in x[1:-1])]
            # TODO: take shortest?
            flows.append(paths[0][:-1])

        # now route the external loop_momenta to the sink
        ext_flows = []
        for e in self.ext[:-1]:
            paths = self.find_path(e, self.ext[-1])
            paths = [x for x in paths if all(
                y[0] not in loop_momenta for y in x[1:-1])]
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
                    overall_sign = 1 if y[0][1] else -1
                    for yy in y:
                        if yy[0] == i:
                            prop_sign = 1 if yy[1] else -1
                            s1 += ("+" if overall_sign * prop_sign == 1 else "-") + \
                                "p{}".format(self.ext[j])
                            newmom.append((self.ext[j], True if prop_sign * overall_sign == 1 else False))
                            break
                mom.append(tuple(newmom))
                s.append(s1 + ")^2")

        self.propagators = mom
        #print("Constructed topology: {}".format('*'.join(s)))

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
        contour_sign = numpy.prod([-1 for i in close_contour if i == 0])
        cut_struct_sign = numpy.prod(sigmas)
        det_sign = numpy.linalg.det(allowed_system[2])
        residue = [[allowed_system[0], sigmas, det_sign*cut_struct_sign*contour_sign]]
        thetas = []
        for n_iter in range(self.n_loops):
            theta = [0 for i in range(self.n_loops)]
            sub_matrix = numpy.array(
                [[allowed_system[2][i][j] for j in range(n_iter+1)] for i in range(n_iter+1)])
            if n_iter == 0:
                theta[allowed_system[1][n_iter]] = 1. / \
                    numpy.linalg.det(sub_matrix) * \
                    sigmas[allowed_system[1][n_iter]]
            else:
                for r in range(n_iter+1):
                    subsub_matrix = numpy.array([[allowed_system[2][i][j] for j in range(
                        n_iter)] for i in range(n_iter+1) if i != r])
                    theta[allowed_system[1][r]] = (-1.)**(n_iter+r)*numpy.linalg.det(
                        subsub_matrix)/numpy.linalg.det(sub_matrix)
                    theta[allowed_system[1][r]] *= sigmas[allowed_system[1][r]]
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
        if not all(all(theta==[] for theta in residue[1]) for residue in residues):
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
        if self.n_loops == 0:
            return []

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

    def create_loop_topology(self, name, ext_mom, mass_map={}, loop_momenta_names=None,
            contour_closure=None, analytic_result=None, fixed_deformation=None, constant_deformation=None,
            loop_momentum_map=None, shift_map=None, numerator_tensor_coefficients=None):
        if loop_momentum_map is None:
            self.generate_momentum_flow(loop_momenta_names)


        # collect all loop lines and construct the signature
        loop_line_map = defaultdict(list)
        loop_line_vertex_map = defaultdict(list)

        # since edges could be flipped, we create an new shift map
        new_shift_map = copy.deepcopy(shift_map)

        for prop, (edge_name, v1, v2) in zip(self.propagators, self.edge_map_lin):
            if prop == ():
                # external momentum
                continue

            mass = 0. if edge_name not in mass_map else mass_map[edge_name]

            # construct the signature
            signature = [0]*len(self.loop_momenta)
            q = vectors.LorentzVector([0., 0., 0., 0.])
            for (mom, sign) in prop:
                s = 1 if sign else -1
                if mom not in self.ext:
                    signature[self.loop_momenta.index(mom)] = s
                else:
                    q += ext_mom[self.edge_map_lin[mom][0]] * s

            # we keep the direction of the loop momenta of the graph in the loop graph, so we need to flip
            # all edges that depend on one loop momentum and have negative sign
            should_flip = len([e for e in signature if e != 0]) == 1 and next(e == -1 for e in signature if e != 0)

            # flip the sign if the inverse exists
            # TODO: should this generate a minus sign when there is an odd-power numerator?
            alt_sig = tuple(s * -1 for s in signature)
            if len(signature) > 0 and alt_sig in loop_line_map or should_flip:
                print('warning: changing sign of propagator %s: %s -> %s' % (edge_name, tuple(signature), alt_sig) )
                loop_line_map[alt_sig].append((edge_name, -q, mass))
                loop_line_vertex_map[alt_sig] += [(v2, v1)]

                if shift_map is not None:
                    new_shift_map[edge_name] = [[s * -1 for s in shift_map[edge_name][0]], [s * -1 for s in shift_map[edge_name][1]]]
            else:
                loop_line_map[tuple(signature)].append((edge_name, q, mass))

                loop_line_vertex_map[tuple(signature)] += [(v1, v2)]

        # fuse vertices
        fuse_map = {}
        for sig, edges in loop_line_vertex_map.items():
            # fuse vertices with repeating edges
            vertlist = [v for vs in edges for v in vs]
            assert(len([1 for v in vertlist if vertlist.count(v) > 2]) == 0)
            double_verts = set(v for v in vertlist if vertlist.count(v) == 2)

            for v in double_verts:
                # keep one edge in the one-loop case
                if len(edges) == 1:
                    break
                other = [(1, e[1]) if e[0] == v else (0, e[0]) for e in edges if v in e]

                if len(sig) == 0:
                    # if there is no loop momentum dependence, the edge may appear in both
                    # directions. Make sure we can fuse by adding the other option too
                    other += [(0, e[1]) if e[0] == v else (1, e[0]) for e in edges if v in e]

                # create a new edge with the correct orientation
                newedge = tuple( next(v[1] for v in other if v[0] == i) for i in range(0, 2))
                d = [e for e in edges if v in e]
                edges[edges.index(d[0])] = newedge
                del edges[edges.index(d[1])]

            # if we have multiple sets left, they are disconnected parts in the graph
            # we shrink the vertices of all but the first group, that will be the representative
            for edge in edges[1:]:
                assert(edge[0] not in fuse_map)
                fuse_map[edge[0]] = edge[1]

            loop_line_vertex_map[sig] = (edges[0][0], edges[0][1])

        # vertices that are fused may again be fused with another vertex
        def multifuse(v):
            if v in fuse_map:
                return multifuse(fuse_map[v])
            else:
                return v

        # now fuse the vertices in the map
        for sig, edges in loop_line_vertex_map.items():
            loop_line_vertex_map[sig] = tuple(multifuse(v) for v in edges)

        loop_line_list = list(loop_line_map.items())
        ll = [LoopLine(
            start_node=loop_line_vertex_map[signature][0],
            end_node=loop_line_vertex_map[signature][1],
            signature=signature,
            propagators=tuple(
                Propagator(q=q, m_squared=mass**2,
                            power=self.powers[edge_name],
                            parametric_shift=new_shift_map[edge_name] if shift_map is not None else None,
                            name=edge_name)
                for (edge_name, q, mass) in propagators)) for signature, propagators in loop_line_list]

        cs = self.get_cut_structures([l for l in ll if any(s != 0 for s in l.signature)], contour_closure)

        # pad the cut structure with the removed edges
        map_to_full_graph = [i for i,  (signature, _) in enumerate(loop_line_list) if not all(sig == 0 for sig in signature)]
        full_cs = []
        for c in cs:
            full_c = [0]*len(loop_line_list)
            for i, cc in enumerate(c):
                full_c[map_to_full_graph[i]] = cc
            full_cs.append(full_c)
        cs = full_cs

        if cs == []:
            # if we have a tree graph, still add a cut structure so that the loop lines are evaluated
            cs = [[0]*len(loop_line_list)]

        # the external kinematics are stored in the right order
        try:
            external_kinematics = [ext_mom["q%d"%n] for n in sorted([int(qi.replace("q","")) for qi in ext_mom.keys()])]
        except ValueError:
            print("External kinematics are not labels as q1,...,qn. The order may be random.")
            external_kinematics = list(ext_mom.values())
        
        loop_topology = LoopTopology(name=name, n_loops=len(self.loop_momenta), external_kinematics=external_kinematics,
            ltd_cut_structure=cs, loop_lines=ll, analytic_result = analytic_result, fixed_deformation = fixed_deformation,
            constant_deformation = constant_deformation, loop_momentum_map=loop_momentum_map)

        if analytic_result is None:
            loop_topology.analytic_result = self.guess_analytical_result(self.loop_momenta, ext_mom, mass_map)

        loop_topology.numerator_tensor_coefficients = numerator_tensor_coefficients
       
        if (fixed_deformation is None) or (fixed_deformation is True):
            # TODO generate the source coordinates automatically with cvxpy if 
            # not specified.
            loop_topology.build_fixed_deformation(force=(fixed_deformation is True))
        
        if constant_deformation is None:
            loop_topology.build_constant_deformation()

        return loop_topology

    def guess_analytical_result(self, loop_momenta, ext_mom, masses):
        """ Try and guess the analytic value for this particular loop topology."""

        # We do not have a guess yet beyond one loop.
        if self.n_loops!=1:
            return None

        loop_edge =self.edges[loop_momenta[0]]
        sorted_external_momenta = []
        start_vertex = loop_edge[0]
        curr_vertex = loop_edge[0]
        curr_edge_index = loop_momenta[0]
        next_vertex = loop_edge[1]
        found_error = False
        ordered_external_momenta = []
        ordered_masses = []
        while True:

            # External momenta, incoming
            next_ext_mom = vectors.LorentzVector()
            for i_ext in self.ext:
                if curr_vertex in self.edges[i_ext]:
                    if curr_vertex==self.edges[i_ext][1]:
                        next_ext_mom += ext_mom[self.edge_map_lin[i_ext][0]]
                    else:
                        next_ext_mom -= ext_mom[self.edge_map_lin[i_ext][0]]
            ordered_external_momenta.append(next_ext_mom)
             
            # Add mass
            ordered_masses.append(masses.get(self.edge_map_lin[curr_edge_index][0],0.))
            
            if next_vertex == start_vertex:
                break

            # Find the next edge
            found_it = False
            for i_edge, edge in enumerate(self.edges):
                if (i_edge not in self.ext) and \
                   (i_edge != curr_edge_index) and \
                   (next_vertex in edge):
                    found_it = True
                    curr_edge_index = i_edge
                    curr_vertex = next_vertex
                    if next_vertex==edge[0]:
                        next_vertex = edge[1]
                    else:
                        next_vertex = edge[0]
                    break

            if not found_it:
                found_error = True
                print("Error in guessing one-loop analytical formula.")
                break
        
        if found_error:
            return None

        # Now compute the corresponding
        try:
            import AVHOneLOopHook.pyAVH_OneLOop_hook as one_loop
        except ImportError:
            print("Run make in LTD/AVHOneLOopHook to generate the python bindings for OneLoop.")
            return None
        
        #print(ordered_external_momenta,ordered_masses)
        # Order masses so as to match OneLOop conventions
        ordered_masses = ordered_masses[1:]+[ordered_masses[0],]
        result = None

        if len(ordered_external_momenta)==3:
            with Silence(active=True):            
                result = one_loop.compute_one_loop_triangle(
                    ordered_external_momenta[0].square(),
                    ordered_external_momenta[1].square(),
                    ordered_external_momenta[2].square(),
                    ordered_masses[0]**2,
                    ordered_masses[1]**2,
                    ordered_masses[2]**2,
                )

        elif len(ordered_external_momenta)==4:
            with Silence(active=True):            
                result = one_loop.compute_one_loop_box(
                    ordered_external_momenta[0].square(),
                    ordered_external_momenta[1].square(),
                    ordered_external_momenta[2].square(),
                    ordered_external_momenta[3].square(),
                    (ordered_external_momenta[0]+ordered_external_momenta[1]).square(),
                    (ordered_external_momenta[1]+ordered_external_momenta[2]).square(),
                    ordered_masses[0]**2,
                    ordered_masses[1]**2,
                    ordered_masses[2]**2,
                    ordered_masses[3]**2,                
                ) 
        
        #print(complex(result[0])*(complex(0.,1.)/(16.0*math.pi**2)))
        if result is None:
           return complex(0.0)
        else:
           return complex(result[0])*(complex(0.,1.)/(16.0*math.pi**2))

#############################################################################################################
# Define topology structures
#############################################################################################################

def solve_constraint_problem(id_and_constraints):
    """ Solve a constraint problem using cvxpy. The input is a tuple (id, constraints). On success, this function
    returns the id. Otherwise, it returns `None`."""
    p = cvxpy.Problem(cvxpy.Minimize(1), id_and_constraints[1])
    p.solve()
    if p.status == cvxpy.OPTIMAL:
        return id_and_constraints[0]

def solve_center_problem(id_and_problem):
    """ Find the center in a region of overlap. Note that all source_coordinates must be used in the constraints."""
    ret_id, extra_constraints, source_coordinates, radius, overlap, ellipsoid_list, delta_param = id_and_problem

    constraints = []

    # do a basis change of this surface to the constrainted basis
    if len(extra_constraints) > 0:
        M = Matrix([cs[:-1] for cs in extra_constraints])
        # pad the constraints with the nullspace to form a basis
        nullspace = nsimplify(M, rational=True).nullspace()
        total_m = M.tolist()
        total_m.extend(r for x in nullspace for r in x.T.tolist())
        total_m = Matrix(total_m)
        # this is the matrix suitable for mapping signature vectors
        inv = total_m.T**-1

        # construct the rhs, note the minus sign
        rhs = [-numpy.copy(cs[-1]) for cs in extra_constraints]
        num_free_momenta = len(nullspace)
        rhs.extend(source_coordinates[:num_free_momenta])
    else:
        num_free_momenta = len(source_coordinates)
        rhs = [x for x in source_coordinates]
        inv = eye(num_free_momenta)

    for overlap_ellipse_ids in overlap:
        (overall_sign, foci) = ellipsoid_list[overlap_ellipse_ids][2]

        for direction in range(len(extra_constraints), len(source_coordinates)):
            # go through all cartesian directions in 3*L space
            for dir_sign in [-3,3,-2,2,-1,1]:
                new_rhs = [x for x in rhs]
                dir_vector = [0, 0, 0]
                dir_vector[abs(dir_sign) - 1] = -1 if dir_sign < 0 else 1
                new_rhs[direction] = new_rhs[direction] + cvxpy.Constant(dir_vector) * radius

                expr = cvxpy.Constant(0)
                is_direction_used = False
                for (sign, delta_index, shift) in foci:
                    (mom_dep, shift1, m1) = delta_param[delta_index]
                    mom = numpy.copy(shift1)
                    for (loop_index, mom_sign) in mom_dep:
                        vec = Matrix([0]*len(source_coordinates))
                        vec[loop_index] = mom_sign

                        for i, (a,b) in enumerate(zip(inv * vec, new_rhs)):
                            if a != 0.:
                                is_direction_used |= i == direction
                                if a == 1.:
                                    mom += b
                                elif a == -1:
                                    mom += (-b)
                                else:
                                    mom += float(a) * b

                    if int(sign) == 1:
                        expr += cvxpy.norm(cvxpy.hstack([m1, mom]), 2) + (-shift)
                    elif int(sign) == -1:
                        expr -= cvxpy.norm(cvxpy.hstack([m1, mom]), 2) + shift
                    else:
                        assert(False)

                if not is_direction_used:
                    # the radius does not appear in this constraints, so we can skip it
                    continue

                if overall_sign < 0:
                    constraints.append(expr >= 0)
                else:
                    constraints.append(expr <= 0)

    # Example of how to force the z-component of the sources to be zero
    #for c in source_coordinates:
    #    constraints.append(c[2]==0.)

    # The radius is not used in any of the constraints, so this problem has no solution
    if constraints == []:
        return (ret_id, overlap, None, None)

    objective = cvxpy.Maximize(radius)
    p = cvxpy.Problem(objective, constraints)
    try:
        p.solve()
        print('ECOS solved problem status: %s and radius %s' % (p.status, radius.value))
    except cvxpy.SolverError:
        print('Solving failed. Trying again with SCS solver')
        p.solve(solver=cvxpy.SCS, eps=1e-9)
        print('SCS solved problem status: %s and radius %s' % (p.status, radius.value))
    except Exception as e:
        print('Error: %s' % e)
        #print("Could not solve system, it should have a solution", p)
        raise

    # rotate sources back
    new_rhs = [x for x in rhs[:len(extra_constraints)]] + [numpy.array([float(c.value[0]), float(c.value[1]), float(c.value[2])]) for c in source_coordinates[:num_free_momenta]]

    res = [numpy.array([0., 0., 0.]) for _ in source_coordinates]
    for i in range(len(source_coordinates)):
        for j in range(len(source_coordinates)):
            # take the transpose of the inverse matrix for the signatures
            res[i] += new_rhs[j] * float(inv[(j,i)])

    return (ret_id, overlap, float(radius.value), [[0., float(x[0]), float(x[1]), float(x[2])] for x in res])

class LoopTopology(object):
    """ A simple container for describing a loop topology."""

    _cvxpy_threshold = 1.0e-12
    # The existence condition itself is inaccurate because of the kinematics being not exactly on-shell
    # and squaring + squareroots decreases the accuracy down to less than 16/2 digits.
    _existence_threshold = 1.0e-7
    def __init__(self, ltd_cut_structure, loop_lines, external_kinematics, n_loops=1, name=None, analytic_result=None,
        fixed_deformation=None, constant_deformation=None, maximum_ratio_expansion_threshold=None, loop_momentum_map=None,
        numerator_tensor_coefficients=None, **opts):
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
        if callable(analytic_result):
            self.analytic_result   = analytic_result(self.external_kinematics)
        else:
            self.analytic_result   = analytic_result

        self.fixed_deformation = fixed_deformation
        self.constant_deformation = constant_deformation
        self.maximum_ratio_expansion_threshold = maximum_ratio_expansion_threshold
        self.loop_momentum_map = loop_momentum_map
        self.numerator_tensor_coefficients = numerator_tensor_coefficients

    def evaluate(self, loop_momenta):
        """ Evaluates Loop topology with the provided list loop momenta, given as a list of LorentzVector."""

        result = 1.0
        for loop_line in self.loop_lines:
            result *= loop_line.evaluate_inverse(loop_momenta)
        return 1.0/result

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
            print("The print function of the LoopTopology requires the package graphviz.")

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
        allowed_export_format = ['yaml', 'mathematica']
        if export_format not in allowed_export_format:
            raise BaseException("Topology can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format == 'yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export topologies to yaml.")

            if output_path is not None:
                open(output_path,'w').write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
            else:
                return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)

        # convert the dual integrands to a mathematica expression
        if export_format == 'mathematica':
            names = ['k', 'l', 'm' ,'n']

            formatted_str = self.name.replace('_','') + '['
            formatted_str += ','.join('%sx_,%sy_,%sz_' % (name, name, name) for name in names[:self.n_loops])
            formatted_str += "]:= Module[{delta},\n\tdelta = Table[0, {i, %s}];\n" % sum(len(ll.propagators) for ll in self.loop_lines)

            # compute all deltas
            prop_count = 1
            for ll in self.loop_lines:
                for p in ll.propagators:
                    formatted_str += '\tdelta[[%s]] = Sqrt[Total[(' % prop_count
                    for sig, name in zip(p.signature, names):
                        if sig == 1:
                            formatted_str += '+{%sx,%sy,%sz}' % (name, name, name)
                        elif sig == -1:
                            formatted_str += '-{%sx,%sy,%sz}' % (name, name, name)
                    formatted_str += '+{%s,%s,%s}' % (p.q[1], p.q[2], p.q[3])
                    formatted_str += ')^2]+' + str(p.m_squared) + '];\n'
                    prop_count += 1

            surface_equations = []
            surface_ids = []

            prop_to_tuple = []
            for (ll_index, ll) in enumerate(self.loop_lines):
                for p_index, p in enumerate(ll.propagators):
                    prop_to_tuple.append((ll_index, p_index))

            # construct all dual integrands
            formatted_str += '\t(-2 Pi I)^%s (1/(2 Pi)^4)^%s(' % (self.n_loops, self.n_loops)
            for cs in self.ltd_cut_structure:
                for cut in product(*[[(c, i, p) for i, p in enumerate(ll.propagators)] if c != 0 else [(0,-1,None)] for c, ll in zip(cs, self.loop_lines)]):
                    # construct the cut basis to loop momentum basis mapping
                    mat = []
                    cut_info = []
                    cut_prop_count = 1
                    for ll, (cut_sign, cut_prop_index_in_ll, cut_prop) in zip(self.loop_lines, cut):
                        if cut_sign != 0:
                            mat.append(ll.signature)
                            cut_info.append((cut_sign, cut_prop_count + cut_prop_index_in_ll, cut_prop.q[0]))
                        cut_prop_count += len(ll.propagators)
                    nmi = numpy.linalg.inv(numpy.array(mat).transpose())

                    prop_count = 1
                    formatted_str += '+1/('
                    for ll_index, (ll, (_, cut_prop_index_in_ll, _)) in enumerate(zip(self.loop_lines, cut)):
                        cut_energy = ''
                        sig_map = nmi.dot(ll.signature)

                        surf_id = []
                        for (sig_sign, (cut_sign, index, shift)) in zip(sig_map, cut_info):
                            if sig_sign != 0:
                                if cut_sign * sig_sign == 1:
                                    cut_energy += '+delta[[%s]]' % index
                                else:
                                    cut_energy += '-delta[[%s]]' % index
                                if sig_sign == 1:
                                    cut_energy += '-(%s)' % shift # add parenthesis to prevent -- operator
                                else:
                                    cut_energy += '+%s' % shift
                                surf_id.append((prop_to_tuple[index - 1], int(cut_sign * sig_sign), -int(sig_sign)))

                        for p_index, p in enumerate(ll.propagators):
                            if cut_prop_index_in_ll == p_index:
                                formatted_str += '2*delta[[%s]]' % prop_count
                            else:
                                formatted_str += '((%s+%s)^2-delta[[%s]]^2)' % (cut_energy, p.q[0], prop_count)
                                # TODO: filter for existence?
                                surface_equations.append('S%s' % len(surface_equations)
                                    + '=%s+%s+delta[[%s]];' % (cut_energy, p.q[0], prop_count))
                                surface_equations.append('S%s' % len(surface_equations)
                                    + '=%s+%s-delta[[%s]];' % (cut_energy, p.q[0], prop_count))

                                surf_id_pos = surf_id + [((ll_index, p_index), +1, +1)]
                                surf_id_pos = tuple(sorted(surf_id_pos))
                                surface_ids.append(('SurfaceID[%s]=%s;' % (len(surface_ids), surf_id_pos)).replace('(','{').replace(')','}'))

                                surf_id_neg = surf_id + [((ll_index, p_index), -1, +1)]
                                surf_id_neg = tuple(sorted(surf_id_neg))
                                surface_ids.append(('SurfaceID[%s]=%s;' % (len(surface_ids), surf_id_neg)).replace('(','{').replace(')','}'))
                            prop_count += 1
                    formatted_str += ')'

                    surf_id = list(sorted(surf_id))
                    if surf_id[0][1] == -1:
                        surf_id = [(x, -a, -b) for x, a, b in surf_id]

            formatted_str += ')];'

            return formatted_str + '\n' + '\n'.join(surface_equations) + '\n'  + '\n'.join(surface_ids)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        res={}
        
        res['name'] = self.name
        res['ltd_cut_structure'] = [list(cs) for cs in self.ltd_cut_structure]
        res['n_loops'] = self.n_loops

        if self.fixed_deformation is not None and self.fixed_deformation:
            res['fixed_deformation'] = self.fixed_deformation
        res['constant_deformation'] = self.constant_deformation
        if self.loop_momentum_map is not None:
            res['loop_momentum_map'] = self.loop_momentum_map
        res['maximum_ratio_expansion_threshold'] = -1.0 if self.maximum_ratio_expansion_threshold is None else self.maximum_ratio_expansion_threshold 
        res['loop_lines'] = [ll.to_flat_format() for ll in self.loop_lines]
        res['external_kinematics'] = [ [float(v) for v in vec] for vec in self.external_kinematics]
        res['analytical_result_real'] = float(self.analytic_result.real) if self.analytic_result else 0.
        res['analytical_result_imag'] = float(self.analytic_result.imag) if self.analytic_result else 0.
        if self.numerator_tensor_coefficients:
            res['numerator_tensor_coefficients'] = [list(c) for c in self.numerator_tensor_coefficients]

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
            analytic_result   =   (None if (flat_dict['analytical_result_real']==0. and flat_dict['analytical_result_imag']==0.) 
                                     else complex(flat_dict['analytical_result_real'],flat_dict['analytical_result_imag'])),
            fixed_deformation =  flat_dict['fixed_deformation'] if 'fixed_deformation' in flat_dict else None,
            maximum_ratio_expansion_threshold = None if ('maximum_ratio_expansion_threshold' not in flat_dict or flat_dict['maximum_ratio_expansion_threshold']<0.) else flat_dict['maximum_ratio_expansion_threshold']
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

    def build_fixed_deformation(self, force=False, bottom_up=False):
        """ This function identifies the fixed deformation sources for the deformation field as well as a the list of
        surfaces ids to exclude for each."""

        # note that all variables need to be created before the pool starts to prevent duplicating variable ids
        source_coordinates = [cvxpy.Variable(3) for _ in range(self.n_loops)]
        radius = cvxpy.Variable(1, nonneg=True)

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = multiprocessing.Pool(1) # use all available cores
        signal.signal(signal.SIGINT, original_sigint_handler)

        ellipsoids, ellipsoid_param, delta_param, expansion_threshold = self.build_existing_ellipsoids(source_coordinates)

        print('Number of ellipsoids %s: %s' % (self.name, len(ellipsoids)))
        print('Maximum expansion threshold allowed: %.3g'%expansion_threshold)
        if self.maximum_ratio_expansion_threshold is None:
            self.maximum_ratio_expansion_threshold = expansion_threshold

        self.fixed_deformation = []

        if len(ellipsoids) == 0:
            return

        center_problems = []

        # We set all combinations of independent propagators on-shell from 0 to n-1 propagators
        for num_cuts in range(self.n_loops):
            for ll_combs in combinations(range(len(self.loop_lines)), num_cuts):
                for prop_combs in product(*[range(len(self.loop_lines[ll].propagators)) for ll in ll_combs]):
                    # construct the on-shell constraints
                    constraints = []
                    extra_constraints = []
                    for (ll_index, prop_index) in zip(ll_combs, prop_combs):
                        q = self.loop_lines[ll_index].propagators[prop_index].q
                        constraint = cvxpy.Constant(numpy.array([q[1], q[2], q[3]]))
                        extra_constraints.append(list(self.loop_lines[ll_index].signature) + [numpy.array([q[1], q[2], q[3]])])
                        for source, sign in zip(source_coordinates, self.loop_lines[ll_index].signature):
                            constraint += source * int(sign)
                        constraints.append(constraint == 0)

                    # check for the independence of the constraints
                    if len(extra_constraints) > 0:
                        _, s, _ = numpy.linalg.svd([cs[:-1] for cs in extra_constraints])
                        if numpy.sum(s > 1e-14) != len(extra_constraints):
                            continue

                    # Filter non-existing ellipsoids under the extra constraints
                    non_existing_ellipsoids = set()
                    for surf_id in list(ellipsoids):
                        e = ellipsoids[surf_id]
                        p = cvxpy.Problem(cvxpy.Minimize(e), [x for x in constraints])
                        result = p.solve()
                        if result > -self._cvxpy_threshold*self.get_com_energy():
                            non_existing_ellipsoids.add(surf_id)

                    if len(non_existing_ellipsoids) > 0:
                        print('Non-existing ellipsoids in limit:', non_existing_ellipsoids)

                    ellipsoid_list = [(surf_id, ellipsoids[surf_id], surf) for surf_id, surf in ellipsoid_param.items() if surf_id not in non_existing_ellipsoids]

                    all_overlaps = self.find_overlap_structure(ellipsoid_list, constraints, bottom_up, pool)
                    print('Overlap structure of %s with cuts %s: %s' % (self.name, list(zip(ll_combs, prop_combs)), all_overlaps))

                    for overlap in all_overlaps:
                        excluded_ellipsoids = list(non_existing_ellipsoids) + list(ellipsoid_list[i][0] for i in range(len(ellipsoid_list)) if i not in overlap)
                        center_problems.append(((tuple(zip(ll_combs, prop_combs)), excluded_ellipsoids), extra_constraints,
                            source_coordinates, radius, overlap, ellipsoid_list, delta_param))

        print('Determining centers of {} cases'.format(len(center_problems)))

        deformation_per_sub_source = defaultdict(list)
        for (prop_combs, excluded_ellipsoids), overlap, radius_value, sources in pool.imap_unordered(solve_center_problem, center_problems):
            if radius_value is None:
                print("No center for {}".format(tuple(prop_combs)))
                continue

            print("Found center for {} with radius {}".format(tuple(prop_combs), radius_value))

            # produce yaml-friendly deformation structure
            d = { 'deformation_sources': sources,
                  'excluded_surface_ids': [[[list(focus[0]), focus[1], focus[2]] for focus in surf_id] for surf_id in excluded_ellipsoids],
                  'radius': radius_value,
                  'overlap': list(overlap),
                }
            deformation_per_sub_source[tuple(prop_combs)].append(d)

        for ll_prop_combs, current_deformation in deformation_per_sub_source.items():
            self.fixed_deformation.append({'deformation_per_overlap': current_deformation,
                                        'excluded_propagators': [list(x) for x in ll_prop_combs]
                                        })

        pool.close()
        pool.join()

        print('Fixed deformation: %s' % pprint.pformat(self.fixed_deformation))

    def build_existing_ellipsoids(self, source_coordinates):

        # Store the shifts and mass from each propagator delta
        deltas = []
        delta_param = []
        for ll in self.loop_lines:
            for p in ll.propagators:
                mom = sum(source_coordinates[i_loop_momentum] * sig for i_loop_momentum, sig in enumerate(p.signature) if sig!=0)
                mom += p.q[1:]
                deltas.append(cvxpy.norm(cvxpy.hstack([math.sqrt(p.m_squared), mom]), 2))
                delta_param.append(([(i_loop_momentum, sig) for i_loop_momentum, sig in enumerate(p.signature) if sig!=0], p.q[1:], math.sqrt(p.m_squared)))

        prop_id_to_ll_map = {}
        counter = 0
        for ll_index, ll in enumerate(self.loop_lines):
            for p_index, p in enumerate(ll.propagators):
                prop_id_to_ll_map[counter] = (ll_index, p_index)
                counter +=1

        # construct all ellipsoid surfaces
        ellipsoids = {}
        ellipsoid_param = {}
        ellipsoid_shift_and_masses_in_cut_basis = {}
        for cs in self.ltd_cut_structure:
            for cut in product(*[[(c, i, p) for i, p in enumerate(ll.propagators)] if c != 0 else [(0,-1,None)] for c, ll in zip(cs, self.loop_lines)]):
                # construct the cut basis to loop momentum basis mapping
                mat = []
                cut_info = []
                cut_prop_count = 0
                for ll, (cut_sign, cut_prop_index_in_ll, cut_prop) in zip(self.loop_lines, cut):
                    if cut_sign != 0:
                        mat.append(ll.signature)
                        cut_info.append((cut_sign, cut_prop_count + cut_prop_index_in_ll, cut_prop.q[0], cut_prop.q.space(), cut_prop.m_squared))
                    cut_prop_count += len(ll.propagators)
                nmi = numpy.linalg.inv(numpy.array(mat).transpose())

                prop_count = 0
                for ll, (_, cut_prop_index_in_ll, _) in zip(self.loop_lines, cut):
                    sig_map = nmi.dot(ll.signature)
                    eq = 0.
                    shift_energy = 0.
                    shift_vector = vectors.Vector([0.,0.,0.])
                    ellipsoid_masses = []
                    eq_fn = []
                    is_ellipsoid = True
                    surf_sign = None

                    sig_sign_map = []
                    for (sig_sign, (cut_sign, index, shift, shift_v, cut_m_squared)) in zip(sig_map, cut_info):
                        if sig_sign != 0:
                            eq += cut_sign * sig_sign * deltas[index]
                            eq += (-sig_sign) * shift
                            shift_energy += (-sig_sign) * shift
                            shift_vector += (-sig_sign) * shift_v
                            ellipsoid_masses.append(cut_m_squared)
                            eq_fn.append([cut_sign * sig_sign, index, sig_sign * shift])
                            sig_sign_map.append((index, sig_sign))

                            if surf_sign is None:
                                surf_sign = cut_sign * sig_sign
                            if surf_sign != cut_sign * sig_sign:
                                is_ellipsoid = False

                    for p_index, p in enumerate(ll.propagators):
                        if cut_prop_index_in_ll != p_index and is_ellipsoid:
                            # create an tuple of tuples that identifies a surface:
                            # [(x_1, a_1, b_1)] where (x,a,b) means a*E_x+b*p_x^0
                            # the sorted list with a=1 is a unique ID
                            if surf_sign == -1.:
                                surf_id = tuple([(prop_id_to_ll_map[x[0]], -int(surf_sign), int(x[1])) for x in sig_sign_map] +
                                        [(prop_id_to_ll_map[prop_count], -int(surf_sign), -1)])
                            else:
                                surf_id = tuple([(prop_id_to_ll_map[x[0]], int(surf_sign), -int(x[1])) for x in sig_sign_map] +
                                        [(prop_id_to_ll_map[prop_count], int(surf_sign), 1)])

                            # make the sign it unique
                            surf_id = tuple(sorted(surf_id))

                            ellipsoid_param[surf_id] = [surf_sign, eq_fn + [(surf_sign, prop_count, -p.q[0])]]
                            if surf_sign == -1.:
                                ellipsoids[surf_id] = eq * -1.0 - p.q[0] + deltas[prop_count]
                                this_ellipsoid_shift_energy = shift_energy*-1.0 - p.q[0]
                                this_ellipsoid_shift_vector = shift_vector*-1.0 - p.q.space()
                            else:
                                ellipsoids[surf_id] = eq + p.q[0] + deltas[prop_count]
                                this_ellipsoid_shift_energy = shift_energy + p.q[0]
                                this_ellipsoid_shift_vector = shift_vector + p.q.space()
                            this_ellipsoid_masses = ellipsoid_masses + [p.m_squared,]

                            ellipsoid_shift_and_masses_in_cut_basis[surf_id] = (this_ellipsoid_shift_energy, this_ellipsoid_shift_vector, this_ellipsoid_masses)

                        prop_count += 1

        # We will compute an upper bound on the expansion threshold here
        expansion_threshold = 1.0
        # Filter non-existing ellipsoids
        for surf_id in list(ellipsoids):
            e = ellipsoids[surf_id]
            p = cvxpy.Problem(cvxpy.Minimize(e), [])
            result = p.solve()
            shift_E, shift_v, masses_squared = ellipsoid_shift_and_masses_in_cut_basis[surf_id]
            mass_term = sum(math.sqrt(m_sq) for m_sq in masses_squared)**2
            existence_equation = shift_E**2 - shift_v.square() - mass_term
            # Be lose in the condition below so as to really only keep true E-surface with non-zero volume
            if existence_equation <= (self._existence_threshold*self.get_com_energy())**2 or shift_E >= -self._existence_threshold*self.get_com_energy():
                if existence_equation < -(self._existence_threshold*self.get_com_energy())**2 and result < -self._cvxpy_threshold*self.get_com_energy():
                    print("WARNING: cvxpy detects the ellipsoid for the following E-surface to be existing even though it does not!")
                    print('> E-surface ID         = %s'%str(surf_id))
                    print('> E-surface params     = %s'%str(ellipsoid_param[surf_id]))
                    print("> cvxpy result         = %.6e < %.6e"%(result, -self._cvxpy_threshold*self.get_com_energy()))
                    print("> shift_E              = %.6e"%shift_E)
                    print("> shift_v              = (%s)"%(', '.join('%.6e'%v_i for v_i in shift_v)))
                    print("> masses_squared       = %s"%(', '.join('%.6e'%m_sq_i for m_sq_i in masses_squared)))
                    print("> existence condition: = %.6e < %.6e or shift_E > %.6e"%
                          (existence_equation, -(self._existence_threshold*self.get_com_energy())**2, -self._existence_threshold*self.get_com_energy()))

                # This is a non-existing E-surface, update the expansion threshold only if it is non-existing *despite*
                # the energy shift having the right sign
                # Be strict in the condition below so as to limit the expansion threshold only for non-existing E-surface which are physical.
                if shift_E < -self._existence_threshold*self.get_com_energy() and existence_equation<-(self._existence_threshold*self.get_com_energy())**2:
                    max_threshold_for_this_non_existing_ellipsoid = math.sqrt(1.0 - math.sqrt(
                            1.0 + existence_equation/(shift_v.square() + mass_term)
                        )
                    )
                    expansion_threshold = min(expansion_threshold, max_threshold_for_this_non_existing_ellipsoid)
                del ellipsoids[surf_id]
                del ellipsoid_param[surf_id]
                continue
            if (existence_equation > (self._existence_threshold*self.get_com_energy())**2) and \
               shift_E < -self._existence_threshold*self.get_com_energy() and \
               result > self._cvxpy_threshold*self.get_com_energy():
                print("WARNING: cvxpy detects the ellipsoid following E-surface to be non-existent even though it does exist!")
                print('E-surface ID         = %s'%str(surf_id))
                print('E-surface params     = %s'%str(ellipsoid_param[surf_id]))
                print("cvxpy result         = %.6e > %.6e "%(result, self._cvxpy_threshold*self.get_com_energy()))
                print("shift_E              = %.6e"%shift_E)
                print("shift_v              = (%s)"%(', '.join('%.6e'%v_i for v_i in shift_v)))
                print("masses_squared       = %s"%(', '.join('%.6e'%m_sq_i for m_sq_i in masses_squared)))
                print("existence condition: = %.6e > %.6e and shift_E < %.6e"%
                      (existence_equation, (self._existence_threshold*self.get_com_energy())**2, -self._existence_threshold*self.get_com_energy() ))

        return ellipsoids, ellipsoid_param, delta_param, expansion_threshold

    def find_overlap_structure(self, ellipsoid_list, extra_constraints, bottom_up, pool):
        if bottom_up:
            return self.find_overlap_structure_bottom_up(ellipsoid_list, extra_constraints)
        else:
            return self.find_overlap_structure_top_down(ellipsoid_list, extra_constraints, pool)

    def find_overlap_structure_top_down(self, ellipsoid_list, extra_constraints, pool):
        overlap_structure = []
        indices = list(range(len(ellipsoid_list)))

        # first collect basic overlap info for all pairs
        pair_non_overlap = [set() for _ in indices]
        pair_overlap = [{i,} for i in indices]
        for es in pool.imap_unordered(solve_constraint_problem, ((r, [ellipsoid_list[e][1] <= 0 for e in r] + extra_constraints) for r in combinations(indices, 2))):
            if es:
                pair_overlap[es[0]].add(es[1])
                pair_overlap[es[1]].add(es[0])

        for i in indices:
            pair_non_overlap[i] = set(indices) - pair_overlap[i]

        for n in range(len(ellipsoid_list), 0, -1):
            # Construct all plausible new subsets of length n
            ellipsoids_not_in_overlap = [e for e in indices if all(e not in o for o in overlap_structure)]

            # take one element from the ellipsoids that are not in an overlap yet and pad it with all possible
            # other ellipsoids that overlap with it
            suggested_options = set()
            for e in ellipsoids_not_in_overlap:
                compatible_completion = [i for i in indices if i != e and i not in pair_non_overlap[e]]
                suggested_options.update(tuple(sorted(set({e,}) | set(a))) for a in combinations(compatible_completion, n - 1))

            # construct all possible 2-tuples from elements that are in two overlaps but are not fully in any overlap structure
            if n > 1:
                disjoint_set_ellipsoids_options = set()
                for o1, o2 in combinations(overlap_structure, 2):
                    left_set = set(o1).difference(set(o2))
                    right_set = set(o2).difference(set(o1))

                    for a in set(product(*[left_set, right_set])):
                        # filter for compatiblity
                        if all(not set(a).issubset(o) for o in overlap_structure) and a[0] not in pair_non_overlap[a[1]]:
                            disjoint_set_ellipsoids_options.add(tuple(sorted(a)))

                for a in disjoint_set_ellipsoids_options:
                    compatible_completion = [i for i in indices if i != a[0] and i != a[1] and i not in pair_non_overlap[a[0]] and i not in  pair_non_overlap[a[1]]]
                    suggested_options.update(tuple(sorted(set(a) | set(b))) for b in combinations(compatible_completion, n - 2))

            # filter any subsets of the overlap structure that are still there
            options = set()
            for r in suggested_options:
                seen = False
                for y in overlap_structure:
                    if set(r).issubset(y):
                        seen = True
                        break
                if not seen:
                    # an option is only possible if all pairs of ellipsoids overlap
                    if all(len(pair_non_overlap[e] & set(r)) == 0 for e in r):
                        options.add(tuple(r))

            print('Progress: n={}, options to consider={} current structure={}'.format(n, len(options), overlap_structure))

            if len(options) == 0:
                continue

            id_and_constraints = [(x, [ellipsoid_list[e][1] <= 0 for e in x] + extra_constraints) for x in options]

            for r in pool.imap_unordered(solve_constraint_problem, id_and_constraints):
                if r and r not in overlap_structure:
                    overlap_structure.append(r)

        return list(sorted(tuple(i) for i in overlap_structure))

    def find_overlap_structure_bottom_up(self, ellipsoid_list, extra_constraints):
        # Find the maximal overlap between all ellipsoids
        # This algorithm is only fast enough for simple cases
        max_overlap = [(i,) for i in range(len(ellipsoid_list))]

        # Already group together ellipsoids that share a focus
        # For these matches, we do not need to test
        grouped = {}
        for i, (surf_id, _, _) in enumerate(ellipsoid_list):
            for focus, _, _ in surf_id:
                if focus not in grouped:
                    grouped[focus] = {i}
                else:
                    grouped[focus].add(i)

        done_overlaps = []
        for n in range(1, len(ellipsoid_list)):
            unused = set(max_overlap)
            new_max_overlap = []
            print('Processing overlap level %s with %s options' % (n, len(max_overlap)))

            current_indices = [i for o in max_overlap for i in o]
            # do a quick filter of the indices
            possible_indices = [i for i in set(current_indices) if current_indices.count(i) >= n]

            for x in combinations(possible_indices, n + 1):
                # if we know they are connected, simply add them
                connected = False
                for groups in grouped.values():
                    if set(x).issubset(groups):
                        connected = True
                        new_max_overlap.append(x)
                        for y in list(unused):
                            if set(y).issubset(x):
                                unused.remove(y)
                        break
                if connected:
                    continue

                # now check if it's possible by checking if all combinations of the
                # indices appear
                possible = True
                for y in combinations(x, n):
                    if y not in max_overlap:
                        possible = False
                        break

                if not possible:
                    continue

                p = cvxpy.Problem(cvxpy.Minimize(1), [ellipsoid_list[e][1] <= 0 for e in x] + extra_constraints)
                result = p.solve()
                if p.status == cvxpy.OPTIMAL:
                    new_max_overlap.append(x)
                    for y in list(unused):
                        if set(y).issubset(x):
                            unused.remove(y)

            for x in unused:
                done_overlaps.append(x)

            max_overlap = new_max_overlap

        return list(sorted(max_overlap + done_overlaps))

    def build_constant_deformation(self):
        # deform with constant norm(a * k + b)

        # TODO compute the right deformation vectors for n_legs <= 4
        # TODO this should come together with an automatic boosting in a proper frame for such topologies
        self.constant_deformation = {
            'alpha' : [[0., 0., 0., 1.],]*self.n_loops,
            'beta'  : [[0., 0., 0., 0.],]*self.n_loops
        }

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

    def evaluate_inverse(self, loop_momenta):
        """ Evaluates the inverse of this loop line with the provided list loop momenta, given as a list of LorentzVector."""

        result = 1.0
        for prop in self.propagators:
            result *= prop.evaluate_inverse(loop_momenta)
        return result

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

    def __init__(self, q, m_squared, power=1, signature = None, name = None, parametric_shift = None, **opts):
        self.name       = name
        self.q          = q
        self.m_squared  = m_squared
        self.power      = power
        # Note that this signature member is not strictly speaking necessary as it should always be the
        # same as the signature attribute of the LoopLine containing it. We forward it here for convenience however.
        self.signature  = signature
        self.parametric_shift = parametric_shift
        self.name = name

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        res={}

        res['q'] = [float(v) for v in self.q]
        res['m_squared'] = self.m_squared
        res['power'] = self.power

        if self.parametric_shift is not None:
            res['parametric_shift'] = self.parametric_shift

        if self.name is not None:
            res['name'] = self.name

        return res

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        return Propagator(
            q           =   vectors.LorentzVector(flat_dict['q']),
            m_squared   =   flat_dict['m_squared'],
            power       =   flat_dict['power'],
            name        =   flat_dict['name'],
            parametric_shift = flat_dict['parametric_shift']
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

class SquaredTopologyGenerator:
    def __init__(self, edges, name, incoming_momentum_names, n_cuts, external_momenta, final_state_particle_ids=(), loop_momenta_names=None, masses={}, powers=None, particle_ids={}):
        self.name = name
        self.topo = TopologyGenerator(edges, powers)
        self.topo.generate_momentum_flow(loop_momenta_names)
        self.external_momenta = external_momenta
        self.cuts = self.topo.find_cutkosky_cuts(n_cuts, incoming_momentum_names, final_state_particle_ids, particle_ids)
        self.masses = masses
        self.incoming_momenta = incoming_momentum_names
        self.topologies = [self.topo.split_graph([a[0] for a in c], incoming_momentum_names) for c in self.cuts]

        self.topo.generate_momentum_flow(loop_momenta_names)
        edge_map = self.topo.get_signature_map()

        self.loop_topologies = []
        self.cut_signatures = []
        self.cut_symmetry_factors = []
        for sub_graphs, c in zip(self.topologies, self.cuts):
            # determine the signature of the cuts
            self.cut_signatures.append(
                [copy.deepcopy(edge_map[cut_edge]) for cut_edge, _ in c]
            )

            loop_topos = []
            for i, s in enumerate(sub_graphs):
                (loop_mom_map, shift_map) = self.topo.build_proto_topology(s, [a[0] for a in c])
                loop_topos.append(
                    s.create_loop_topology(name + '_' + ''.join(a[0] for a in c) + '_' + ['l', 'r'][i],
                    # provide dummy external momenta
                    ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                    fixed_deformation=False,
                    mass_map=masses,
                    loop_momentum_map=loop_mom_map,
                    shift_map=shift_map))

            # compute external state symmetry factor
            cutkosky_particles = tuple(sorted(particle_ids[e] if e in particle_ids else 'NONE' for (e, _) in c))
            sym_factor = 1.
            for x in set(cutkosky_particles):
                sym_factor *= math.factorial(cutkosky_particles.count(x))

            self.cut_symmetry_factors.append(sym_factor)

            self.loop_topologies.append(loop_topos)

    def export(self, output_path):
        out = {
            'name': self.name,
            'n_loops': self.topo.n_loops,
            'n_incoming_momenta': len(self.incoming_momenta),
            'external_momenta': [self.external_momenta["q%d"%n] for n in sorted([int(qi.replace("q","")) for qi in self.external_momenta.keys()])],
            'topo': [list(x) for x in self.topo.edge_map_lin],
            'loop_momentum_basis': [self.topo.edge_map_lin[e][0] for e in self.topo.loop_momenta],
            'e_cm_squared': sum(self.external_momenta[e][0] for e in self.incoming_momenta)**2 - sum(x*x for x in (sum(self.external_momenta[e][i] for e in self.incoming_momenta) for i in range(1, 4))),
            'cutkosky_cuts': [
                {'cuts':
                    [{
                       'name': a[0],
                       'sign': a[1],
                       'signature': sig,
                       'power': self.topo.powers[a[0]],
                       'm_squared': self.masses[a[0]]**2 if a[0] in self.masses else 0.
                    }
                    for a, sig in zip(c, cut_sig)],
                'subgraph_left': ts[0].to_flat_format(),
                'subgraph_right': ts[1].to_flat_format(),
                'symmetry_factor': sym
                }

                for c, cut_sig, ts, sym in zip(self.cuts, self.cut_signatures, self.loop_topologies, self.cut_symmetry_factors)
            ]
        }

        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        open(output_path,'w').write(yaml.dump(out, Dumper=Dumper))

if __name__ == "__main__":
    #triangle = TopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3),
    #                            ('p3', 2, 3), ('q2', 3, 4), ('q3', 2, 5)])  # triangle
    #doubletriangle = TopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 3), ('p3', 2, 3),
    #                                    ('p4', 3, 4), ('p5', 2, 4), ('q2', 4, 5)])  # double-triangle
    #loop_topology = doubletriangle.create_loop_topology("DoubleTriangle", ext_mom={'q1': vectors.LorentzVector(
    #    [0.1, 0.2, 0.3, 0.4])}, mass_map={'p1': 1.0, 'p2': 2.0, 'p3': 3.0}, loop_momenta_names=('p1', 'p5'), analytic_result=None)

    # Construct a cross section
    # result is -2 Zeta[3] 3 Pi/(16 Pi^2)^3 = -5.75396*10^-6
    mercedes = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 6),
                                        ('p4', 6, 5), ('p5', 5, 1), ('p6', 2, 4), ('p7', 3, 4), ('p8', 4, 5), ('q2', 6, 7)], "M", ['q1'], 2,
                                        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
                                        loop_momenta_names=('p1', 'p2', 'p3'),
                                        particle_ids={'p%s' % i: i for i in range(9)})
    mercedes.export('mercedes_squared.yaml')

    # result is -5 Zeta[5] 4 Pi/(16 Pi^2)^4 = -1.04773*10^-7
    doublemercedes = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 7), ('p3', 7, 3), ('p4', 3, 6),
                                        ('p5', 6, 5), ('p6', 5, 1), ('p7', 2, 4), ('p8', 3, 4), ('p9', 4, 5), ('p10', 7, 4), ('q2', 6, 8)], "DM", ['q1'], 2,
                                        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
                                        loop_momenta_names=('p1', 'p2', 'p3', 'p4'),
                                        particle_ids={'p%s' % i: i for i in range(11)})
    doublemercedes.export('doublemercedes_squared.yaml')

    bubble = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 2), ('q2', 2, 3)], "B", ['q1'], 0,
    {'q1': [2., 0., 0., 0.], 'q2': [2., 0., 0., 0.]},
    masses={'p1': 0.24, 'p2': 0.24})
    bubble.export('bubble_squared.yaml')

    t1 = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 4, 3), ('p4', 4, 1), ('p5', 2, 4), ('q2', 3, 5)], "T", ['q1'], 2,
        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
        particle_ids={'p%s' % i: i for i in range(9)})
            #masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100})
    t1.export('t1_squared.yaml')

    bu = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 3, 2), ('p3', 4, 3),
                                        ('p4', 4, 1), ('p5', 2, 5), ('p6', 5, 4), ('p7', 3, 5), ('q2', 3, 6)], "BU", ['q1'], 2,
                                        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
                                        loop_momenta_names=('p2', 'p4', 'p7'),
                                        particle_ids={'p%s' % i: i for i in range(9)})
    bu.export('bu_squared.yaml')

    insertion = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 2, 3), ('p4', 3, 4), ('p5', 1, 4), ('q2', 4, 5)], "I", ['q1'], 3,
        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
        masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100})#, loop_momenta_names=('p4', 'p3'), powers={'p3': 2})
    insertion.export('insertion_squared.yaml')

    # TODO: whether it's t and tbar should be determined from the cutkosky cut direction and not hardcoded in the topology
    # NOTE: for 2 -> N, the first two entries need to be the two incoming momenta
    # the outgoing momenta will be set to the input momenta in the same order, i.e., q3=q1, q4=q2.
    tth = SquaredTopologyGenerator([('q1', 0, 1), ('q2', 6, 7), ('q3', 4, 5), ('q4', 10, 11), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 10),
        ('p5', 10, 9), ('p6', 9, 8), ('p7', 8, 7), ('p8', 1, 7), ('p9', 2, 8), ('p10', 3, 9), ], "TTH", ['q1', 'q2'], 0,
        {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        final_state_particle_ids=('t', 'tbar', 'H'), particle_ids={'p1': 't', 'p2': 't', 'p3': 't', 'p4': 'tbar', 'p5': 'tbar', 'p6': 'tbar', 'p7': 'tbar',
            'p8': 't', 'p9': 'g', 'p10': 'H'})
    tth.export('tth_squared.yaml')

    two_to_two = SquaredTopologyGenerator([('q1', 0, 1), ('q2', 7, 5), ('q3', 2, 8), ('q4', 4, 9), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5),
        ('p5', 5, 6), ('p6', 6, 1), ('p7', 6, 3), ], "two_to_two", ['q1', 'q2'], 3,
        {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        masses={'p1': 100, 'p2':100, 'p3': 100, 'p4': 100, 'p5': 100, 'p6': 100, 'p7': 100}, loop_momenta_names=('p1', 'p7'),)
    two_to_two.export('two_to_two_squared.yaml')

    twoI_twoF = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 2), ('q3', 4, 104), ('q4', 3, 103), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 1),], "twoI_twoF", ['q1', 'q2'], 2,
        {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        masses={'p1': 0., 'p2': 1., 'p3': 0., 'p4': 1.}, loop_momenta_names=('p1',))
    twoI_twoF.export('twoI_twoF_squared.yaml')

    two_to_three = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 2), ('q3', 6, 103), ('q4', 5, 104), ('p1', 2, 3), ('p2', 3, 4),
        ('p3', 3, 4), ('p4', 4, 5), ('p5', 5, 6), ('p6', 6, 1), ('p7', 1, 2)], "two_to_three", ['q1', 'q2'], 3,
        {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        masses={'p1': 1.0, 'p2': 1.0, 'p3': 1.0, 'p4': 1.0, 'p6': 1.0},
        particle_ids={'p2': 1, 'p3': 2, 'p6': 3})
    two_to_three.export('two_to_three_squared.yaml')
