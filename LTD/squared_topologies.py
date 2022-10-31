import os
import sys
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, os.path.abspath( os.path.join(root_path,os.path.pardir) ) )
sys.path.insert(0, root_path)

from ltd_utils import TopologyGenerator
import copy
import math
from itertools import combinations_with_replacement, product
import vectors
import numpy
from sympy import Matrix, diag
from pprint import pprint, pformat

class SquaredTopologyGenerator:
    def __init__(self, edges, name, incoming_momentum_names, n_jets, external_momenta, final_state_particle_ids=(),
        loop_momenta_names=None, loop_momenta_signs=None, masses={}, powers=None, particle_ids={}, jet_ids=None,
        subgraphs_info={},overall_numerator=1., numerator_structure={},
        cut_filter=set(), FORM_integrand={},
        vertex_weights={}, edge_weights={}, generation_options={},analytic_result=None,
        default_kinematics=None, uv_test=None):

        self.name = name
        self.topo = TopologyGenerator(edges, powers)
        self.topo.generate_momentum_flow(loop_momenta_names)
        self.external_momenta = external_momenta
        self.FORM_integrand = FORM_integrand
        self.subgraphs_info = subgraphs_info
        self.generation_options = generation_options

        self.default_kinematics = default_kinematics

        self.particle_ids = particle_ids
        # The edge #i of the LMB may not always carry k_i but sometimes -k_i.
        # This is supported by adjusting the cb to lmb rotation matrix to be applied
        # before calling the numerator.
        self.loop_momenta_signs = loop_momenta_signs
        # However, we no longer want to support this case and instead enforce momenta to always follow
        # the edge orientation. The numerator must therefore be modified upstream so as to satisfy this requirement.
        assert(loop_momenta_signs is None or all(lms==1 for lms in loop_momenta_signs))

        self.loop_topo = self.topo.create_loop_topology(name,
            external_momenta,
            loop_momenta_names=loop_momenta_names,
            fixed_deformation=False,
            mass_map=masses,
            shift_map=None,
            analytic_result=analytic_result
            )

        cutkosky_cuts = self.topo.find_cutkosky_cuts(n_jets, incoming_momentum_names, final_state_particle_ids, 
                                    particle_ids, PDGs_in_jet=jet_ids, fuse_repeated_edges=True)

        if len(cut_filter) > 0:
            new_cutkosky_cuts = [ c for c in cutkosky_cuts if tuple(sorted([e_name for e_name, _ in c])) in [ tuple(sorted(cf)) for cf in cut_filter ] ]
            if len(new_cutkosky_cuts)==0 and len(cutkosky_cuts)>0:
                print("WARNING: could not find selected cuts for %s.\nExisting cuts: %s\nSelected cuts: %s"%(
                    self.name, pformat( [ tuple(sorted([e_name for e_name, _ in c])) for c in cutkosky_cuts ] ),
                    pformat( [ tuple(sorted(cf)) for cf in cut_filter ] )
                ))
            cutkosky_cuts = new_cutkosky_cuts

        self.masses = copy.deepcopy(masses)
        self.overall_numerator = overall_numerator
        self.incoming_momenta = incoming_momentum_names

        edge_map = self.topo.get_signature_map()
        mu_uv = 2. * math.sqrt(sum(self.external_momenta[e][0] for e in self.incoming_momenta)**2 - 
                        sum(x*x for x in (sum(self.external_momenta[e][i] for e in self.incoming_momenta) for i in range(1, 4))))

        diagram_set_counter = 0
        graph_counter = 0
        self.cuts = []
        cut_infos = []
        for cutkosky_cut in cutkosky_cuts:
            # detect external bubbles
            duplicate_cut_edges_set = set()
            cut_powers = {}
            for c, _ in cutkosky_cut:
                sig = tuple(self.topo.propagators[self.topo.edge_name_map[c]])
                inv_sig = tuple([(x[0], not x[1]) for x in sig])
                bubble_props = [eml[0] for prop, eml in zip(self.topo.propagators, self.topo.edge_map_lin) if eml[0] != c and abs(particle_ids[eml[0]]) ==
                    particle_ids[c] and (tuple(prop) == sig or tuple(prop) == inv_sig)]
                cut_powers[c] = len(bubble_props) + 1
                duplicate_cut_edges_set |= set(bubble_props)

            # sort cutkosky cuts such that they align best with the lmb and such that cut edge with the highest power is last
            cutkosky_cut = tuple(sorted(cutkosky_cut, key=lambda x: cut_powers[x[0]] * 1000 + (self.topo.edge_name_map[x[0]] in self.topo.loop_momenta)*-10 + cutkosky_cut.index(x)))

            # split the graph along the Cutkosky cut
            graphs = self.topo.split_graph([a[0] for a in cutkosky_cut], incoming_momentum_names)
            graphs[0].conjugate = False
            graphs[1].conjugate = True

            # TODO: instruct the C code to drop the imaginary part of the bubble forest
            if len(duplicate_cut_edges_set) > 2:
                print("Iterated bubble: check multiplicity factor. Wrong result expected when one of the bubbles needs a deformation as im is not set to 0 yet!")

            # factorize external bubbles
            for d in duplicate_cut_edges_set:
                split_graphs = []
                for g in graphs:
                    if d in g.edge_name_map:
                        sg = g.split_graph([d], incoming_momentum_names)
                        sg[0].conjugate = g.conjugate
                        sg[1].conjugate = g.conjugate
                        split_graphs.extend(sg)
                    else:
                        split_graphs.append(g)

                graphs = split_graphs

            cut_info = {
                'cuts': [{
                        'edge': c[0],
                        'sign': c[1],
                        'power': cut_powers[c[0]],
                        'signature': copy.deepcopy(edge_map[c[0]]),
                        'particle_id': self.particle_ids[c[0]] if c[0] in self.particle_ids else 0,
                        }
                    for c in cutkosky_cut],
                'diagram_info': [{
                    'graph':  g,
                    'conjugate_deformation': g.conjugate,
                } for g in graphs]
            }

            cut_infos.append(cut_info)

        # determine supergraph collinear thresholds in the cmb
        self.sg_collinear = []
        for i1, c1 in enumerate(cut_infos):
            c1_edges = set(a['edge'] for a in c1['cuts'])
            for c2 in cut_infos[:i1]: # TDOD: changed order
                c2_edges = set(a['edge'] for a in c2['cuts'])
                # one overlap is required
                if len(c1_edges & c2_edges) == 0:
                    continue

                # check for crossing
                c1_left = set(self.topo.split_graph(c1_edges, incoming_momentum_names)[0].edge_name_map.keys())
                c2_left = set(self.topo.split_graph(c2_edges, incoming_momentum_names)[0].edge_name_map.keys())
                if not c1_left.issubset(c2_left) and not c2_left.issubset(c1_left):
                    continue


                prop_list = [name for i, (name, v1, v2) in enumerate(self.topo.edge_map_lin) if i not in self.topo.ext]

                col = []
                for col_edge in c1_edges.symmetric_difference(c2_edges):
                    edge_info = next(cc for cc in c1['cuts'] if cc['edge'] == col_edge) if col_edge in c1_edges else next(cc for cc in c2['cuts'] if cc['edge'] == col_edge)

                    col.append((prop_list.index(col_edge), 1 if col_edge in c1_edges else -1))
                col.sort()

                if col not in self.sg_collinear:
                    self.sg_collinear.append(col)

        pure_forest_counter = 0

        uv_representative_graphs = {}

        for cut_info in cut_infos:
            c = cut_info['cuts']
            cut_name = tuple(a['edge'] for a in c)

            # add the uv structure to the diagram
            for i, diag_info in enumerate(cut_info['diagram_info']):
                diag_info['graph'].inherit_loop_momentum_basis(self.topo)

                vw = {v: vertex_weights[v] if v in vertex_weights else 0 for e in diag_info['graph'].edge_map_lin for v in e[1:]}
                ew = {e: edge_weights[e] if e in edge_weights else -2 for e, _, _ in diag_info['graph'].edge_map_lin}
                
                uv_limits = diag_info['graph'].construct_uv_limits(vw, ew, particle_ids=particle_ids,
                            UV_min_dod_to_subtract=self.generation_options.get('UV_min_dod_to_subtract',0), sg_name=self.name)

                subgraph_counter = {}
                for uv_limit in uv_limits:
                    for uv_sg in uv_limit['uv_subgraphs']:
                        subgraph_internal_edges = [e[0] for i, e in enumerate(uv_sg['graph'].edge_map_lin) if i not in uv_sg['graph'].ext]

                        uv_sg['onshell'] = None
                        if len(uv_sg['graph'].edge_map_lin) - len(subgraph_internal_edges) == 2:
                            ext_edge = uv_sg['graph'].edge_map_lin[uv_sg['graph'].ext[0]][0]
                            if masses[ext_edge] > 0.:
                                uv_sg['onshell'] = [uv_sg['graph'].edge_map_lin[i][0] for i in uv_sg['graph'].ext]

                    # give every subdiagram a globally unique id
                    for uv_sg in uv_limit['uv_subgraphs']:
                        r = tuple(uv_sg['subgraph_momenta'])
                        if r in subgraph_counter:
                            uv_sg['first_occurrence_id'] = subgraph_counter[r]
                        else:
                            uv_sg['first_occurrence_id'] = graph_counter
                            subgraph_counter[r] = graph_counter

                        uv_sg['id'] = graph_counter
                        graph_counter += 1

                        uv_sg['integrated_ct_id'] = graph_counter
                        graph_counter += 1

                        for dgi, dg in enumerate(uv_sg['derived_graphs']):
                            dg['id'] = graph_counter
                            graph_counter += 1

                            if dg['derivatives'] < uv_sg['taylor_order']:
                                dg['soft_ct_id'] = graph_counter
                                graph_counter += 1

                                if uv_sg['onshell']:
                                    dg['onshell_ct_id'] = graph_counter
                                    graph_counter += 1
                    uv_limit['remaining_graph_id'] = graph_counter
                    graph_counter += 1

                diag_info.pop('graph')
                diag_info['uv'] = [{
                    'uv_subgraphs': uv_limit['uv_subgraphs'],
                    'uv_spinney': [[list(g), dod] for g, dod in uv_limit['spinney']],
                    'uv_vertices': [x for x in uv_limit['uv_vertices']],
                    'uv_propagators': [m for g, _ in uv_limit['spinney'] for m in g],
                    'remaining_graph': uv_limit['remaining_graph'],
                    'remaining_graph_id' : uv_limit['remaining_graph_id'],
                } for uv_limit in uv_limits]

            cut_info['id'] = diagram_set_counter
            diagram_set_counter += 1

            # construct a matrix from the cut basis to the loop momentum basis
            # this is useful if the numerator is specified in the loop momentum basis
            # the matrix will be padded with the loop momentum maps
            cut_to_lmb = [ cut_edge['signature'][0] for cut_edge in c[:-1]]

            for i, diag_info in enumerate(cut_info['diagram_info']):
                for uv_structure in diag_info['uv']:
                    # create the loop topo of the remaing graph
                    (loop_mom_map, shift_map) = self.topo.build_proto_topology(uv_structure['remaining_graph'], c, skip_shift=False)

                    if uv_structure['uv_spinney'] == []:
                        cut_to_lmb.extend([x[0] for x in loop_mom_map])

                    uv_structure['remaining_graph_loop_topo'] = uv_structure['remaining_graph'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                        # provide dummy external momenta
                        ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                        fixed_deformation=False,
                        mass_map=masses,
                        loop_momentum_map=loop_mom_map,
                        shift_map=shift_map,
                        check_external_momenta_names=False,
                        analytic_result=0)
                    uv_structure['remaining_graph_loop_topo'].external_kinematics = []

            lmb_to_cb_matrix = Matrix(cut_to_lmb)
            # The edge #i of the LMB may not always carry k_i but sometimes -k_i.
            # This is supported by adjusting the cb to lmb rotation matrix to be applied
            # before calling the numerator.
            if self.loop_momenta_signs is not None:
                assert(len(self.loop_momenta_signs)==len(cut_to_lmb[0]))
                assert(all(abs(s)==1 for s in self.loop_momenta_signs))
                lmb_to_cb_matrix = lmb_to_cb_matrix*diag(*[int(s) for s in self.loop_momenta_signs])
            cb_to_lmb_matrix = lmb_to_cb_matrix**-1

            cut_info['lmb_to_cb'] = [int(x) for x in lmb_to_cb_matrix]
            cut_info['cb_to_lmb'] = [int(x) for x in cb_to_lmb_matrix]

            # construct the forest matrix that maps the amplitude momenta in the cmb
            # to ones suitable for the spinney
            lmb_offset = len(c) - 1
            for i, diag_info in enumerate(cut_info['diagram_info']):
                amp_loops = diag_info['uv'][0]['remaining_graph'].n_loops

                diag_info['n_loops'] = amp_loops
                diag_info['propagators'] = []
                diag_info['thresholds'] = []
                for uv_structure in diag_info['uv']:
                    # create the LTD representation of the derived UV graph
                    forest_to_lmb = []
                    for uv_subgraph in uv_structure['uv_subgraphs']:
                        for di, d in enumerate(uv_subgraph['derived_graphs']):
                            (loop_mom_map, shift_map) = self.topo.build_proto_topology(d['graph'], c, skip_shift=True)

                            if len(uv_structure['uv_subgraphs']) > 1:
                                # get the representative graph and map the LTD basis momenta to this representation
                                rg, rg_lmmm = uv_representative_graphs[tuple(uv_subgraph['full_subgraph_momenta'])]
                                sm = rg.get_signature_map()
                                # get the signature of the LTD basis momenta in the basis of the representative
                                ltd_mom_sig = [sm[d['graph'].edge_map_lin[e][0]] for e in d['graph'].loop_momenta]

                                basis_shift_map = []
                                loop_mom_map = [] # overwrite the generated loop momentum map
                                for lmb_sig, shift_sig in ltd_mom_sig:
                                    # construct the loop momentum part that should not get expanded
                                    lmc = numpy.array([0]*self.topo.n_loops, dtype=int)
                                    for sign, (rg_lmb_sig, rg_shift_sig) in zip(lmb_sig, rg_lmmm):
                                        assert(all(s == 0 for s in rg_shift_sig))
                                        lmc += sign * numpy.array(rg_lmb_sig)
                                    loop_mom_map.append((lmc, [0]*len(self.topo.ext)))

                                    # construct the shift part that will get expanded
                                    # it is split between loop momenta external to the subgraph
                                    # and the external momenta of the supergraph
                                    lm_shift = numpy.array([0]*self.topo.n_loops, dtype=int)
                                    ext_shift = numpy.array([0]*len(self.topo.ext), dtype=int)

                                    assert(len(shift_sig) == len(rg.ext))
                                    for sign, ext in zip(shift_sig, rg.ext):
                                        # get the lmb dependence of the external edge of the rg graph from the supergraph
                                        sig = edge_map[rg.edge_map_lin[ext][0]]
                                        lm_shift += sign * numpy.array(sig[0], dtype=int)
                                        ext_shift += sign * numpy.array(sig[1], dtype=int)
                                    basis_shift_map.append((lm_shift, ext_shift))
                            else:
                                # the UV subgraph has no subgraphs of itself and therefore contains at least as many lmb
                                # edges as there are loops. We select the basis of this graph, which is a subset of lmb momenta
                                # as a representative for all UV subgraphs that involve the same momenta (including subgraph momenta)
                                basis_shift_map = [[[0]*len(lmp), shp] for lmp, shp in loop_mom_map] # should be all 0
                                uv_representative_graphs[tuple(uv_subgraph['full_subgraph_momenta'])] = (d['graph'], loop_mom_map)

                            if di == 0:
                                forest_to_lmb.extend([x[0] for x in loop_mom_map])

                                # strip external momenta for the graph and compute all propagators and thresholds
                                new_edges = [e for i, e in enumerate(d['graph'].edge_map_lin) if i not in d['graph'].ext]
                                stripped_topo = TopologyGenerator(new_edges)
                                stripped_topo.inherit_loop_momentum_basis(self.topo)
                                sigs = stripped_topo.get_signature_map()

                                # construct all propagators of the amplitude
                                new_props = []
                                for e, v1, v2 in new_edges:
                                    lmb_sig = numpy.array([0]*self.topo.n_loops, dtype=int)
                                    for lm, s in zip(loop_mom_map, sigs[e][0]):
                                        lmb_sig += s * numpy.array(lm[0], dtype=int)

                                    cmb_sig = [int(x) for x in cb_to_lmb_matrix.T @ lmb_sig]
                                    assert(all(s == 0 for s in cmb_sig[len(c) - 1:lmb_offset]) and all(s == 0 for s in cmb_sig[lmb_offset + amp_loops:]))

                                    amp_sig = numpy.array(cmb_sig[lmb_offset:lmb_offset + amp_loops], dtype=int)

                                    if all(s == 0 for s in amp_sig):
                                        continue
                                
                                    ext_shift_in_cmb = (numpy.array(list(cmb_sig[:lmb_offset]) + [0]*amp_loops + list(cmb_sig[lmb_offset + amp_loops:]), dtype=int),
                                            sum(-numpy.array(cc['signature'][1], dtype=int) * s for s, cc in zip(cmb_sig, c[:-1])))

                                    for uv in (False, True):
                                        new_prop = {
                                            'name': e,
                                            'id': -1,
                                            'lmb_sig': (lmb_sig.tolist(), [0]*len(self.topo.ext)),
                                            'amp_sig': amp_sig.tolist(),
                                            'amp_shift_in_cmb_sig': tuple(x.tolist() for x in ext_shift_in_cmb),
                                            'mass': 0. if uv else masses[e],
                                            'uv_mass': uv,
                                            'uv_only': True,
                                        }
                                        neg_sig = ([-x for x in lmb_sig], [0]*len(self.topo.ext))
                                        new_props.append(new_prop)
                                        for p in diag_info['propagators']:
                                            if (p['lmb_sig'] == new_prop['lmb_sig'] or p['lmb_sig'] == neg_sig) and p['mass'] == new_prop['mass'] and p['uv_mass'] == new_prop['uv_mass']:
                                                new_prop['id'] = p['id']
                                                break
                                        else:
                                            new_prop['id'] = len(diag_info['propagators'])
                                            diag_info['propagators'].append(new_prop)

                                # construct all thresholds
                                thresholds = d['graph'].find_thresholds(fuse_repeated_edges=False)
                                for t in thresholds:
                                    foci_no_uv = [next(p['id'] for p in new_props if p['name'] == f and not p['uv_mass']) for f in t[0]]
                                    foci_uv = [next(p['id'] for p in new_props if p['name'] == f and p['uv_mass']) for f in t[0]]
                                    shift = ([0]*self.topo.n_loops, [0]*len(self.topo.ext))

                                    # create a map into the focus basis
                                    m = Matrix([diag_info['propagators'][d]['amp_sig'] for d in foci[:len(foci) - 1]])
                                    dm = m.T * (m * m.T)**-1 # the right inverse should always work
                                    # dependent edge in focus basis, should be integer
                                    sig_in_fb = Matrix(diag_info['propagators'][foci[-1]]['amp_sig']).T * dm

                                    # add a +m and -m to the surfaces of every massive propagator that have a dependency on the external momentum
                                    if uv_subgraph['onshell'] and len(t[1]) == 1:
                                        for sign in (1.,-1.):
                                            threshold = {'foci': foci_no_uv, 'fb_to_cmb': [float(x) for x in dm], 'sig_in_fb': [int(x) for x in sig_in_fb], 'shift_in_lmb_sig': shift, 'mass_shift': sign * masses[uv_subgraph['onshell'][0]]}
                                            if threshold not in diag_info['thresholds']: diag_info['thresholds'].append(threshold)
                                    else:
                                        # not needed when dod = 0
                                        if len(uv_subgraph['derived_graphs']) > 1:
                                            threshold = {'foci': foci_no_uv, 'fb_to_cmb': [float(x) for x in dm], 'sig_in_fb': [int(x) for x in sig_in_fb], 'shift_in_lmb_sig': shift, 'mass_shift': 0.}
                                            if threshold not in diag_info['thresholds']: diag_info['thresholds'].append(threshold)

                                    threshold = {'foci': foci_uv, 'fb_to_cmb': [float(x) for x in dm], 'sig_in_fb': [int(x) for x in sig_in_fb], 'shift_in_lmb_sig': shift, 'mass_shift': 0.}
                                    if threshold not in diag_info['thresholds']: diag_info['thresholds'].append(threshold)

                            # note: the shift map signs may get swapped when edges switch orientation
                            # the basis_shift_map should be unaffected since defining edges are not swapped
                            loop_topo = d['graph'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                                # provide dummy external momenta
                                ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                                fixed_deformation=False,
                                mass_map=masses,
                                loop_momentum_map=loop_mom_map,
                                shift_map=shift_map,
                                check_external_momenta_names=False,
                                analytic_result=0)
                            loop_topo.external_kinematics = []

                            # when soft derivatives have to be taken we require that the propagators of every loop line have the
                            # same mass, such that raising a power in the loop line yields a unique topology
                            if d['derivatives'] > 0 and d['derivatives'] < uv_sg['taylor_order']:
                                if any(len(set(masses[p.name] for p in ll.propagators)) > 1 for ll in loop_topo.loop_lines):
                                    raise AssertionError('WARNING: In supergraph %s not all propagators have the same mass in a loop line for the derivative soft CT graph %s' % (self.name, str(uv_sg['graph'].edge_map_lin)))

                            d['loop_topo_orig_mass'] = copy.deepcopy(loop_topo)
                            # do not drop the shift for the 0th order of the on-shell expansion
                            if not uv_subgraph['onshell'] or d['derivatives'] > 0:
                                for ll in d['loop_topo_orig_mass'].loop_lines:
                                    for prop in ll.propagators:
                                        prop.parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]

                            # merge equal-mass and equal-shift propagators, needed for LTD derivative treatment
                            from collections import defaultdict
                            for ll in d['loop_topo_orig_mass'].loop_lines:
                                unique_props = defaultdict(list)
                                for p in ll.propagators:
                                    unique_props[(p.m_squared, (tuple(p.parametric_shift[0]), tuple(p.parametric_shift[1])))].append(p)

                                for props in unique_props.values():
                                    for p in props[1:]:
                                        d['loop_topo_orig_mass'].propagator_map[p.name] = props[0].name
                                        props[0].power += p.power
                                    ll.propagators = [p for p in ll.propagators if p not in props[1:]]

                            # take the UV limit of the diagram, add the mass and set the parametric shift to 0
                            # collect the parametric shifts of the loop lines such that it can be used to Taylor expand
                            # the UV subgraph
                            uv_loop_lines = []
                            for ll in loop_topo.loop_lines:
                                derived_ll_power = sum(pp.power for pp in ll.propagators)
                                # determine the power of the loop line of the non-derived graph
                                orig_ll_power = sum(uv_subgraph['graph'].powers[pp.name]
                                    + len([1 for _,mn in loop_topo.propagator_map.items() if mn == pp.name])
                                    for pp in ll.propagators)

                                if uv_structure['remaining_graph'].n_loops == 0 and len(uv_structure['uv_subgraphs']) > 0 and \
                                    uv_test == pure_forest_counter:
                                    # raise every loop line so that it has at least 2 propagators by raising the first prop
                                    if orig_ll_power == 1:
                                        prop_pow = 2
                                    else:
                                        prop_pow = 1
                                else:
                                    prop_pow = 1

                                # note: we assume that every propagator has power 1 and that only merging of identical edges raises it
                                uv_loop_lines.append((ll.signature, [(p.name, p.parametric_shift,
                                    (prop_pow if ii == 0 else 1) + len([1 for _,mn in loop_topo.propagator_map.items() if mn == p.name])) for ii, p in enumerate(ll.propagators)], derived_ll_power - orig_ll_power))

                                prop = ll.propagators[0]
                                prop.uv = True
                                prop.m_squared = mu_uv**2
                                prop.power = derived_ll_power + prop_pow - 1
                                prop.parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]
                                ll.propagators = [prop]

                            loop_topo.uv_loop_lines = (uv_loop_lines, basis_shift_map)

                            d['loop_topo_orig_mass'].uv_loop_lines = (uv_loop_lines, basis_shift_map)
                            d['loop_topo'] = loop_topo

                            # check if raised loop lines have a propagator with a shift, otherwise this configuration is impossible
                            d['skip_pf'] = False
                            for i, (ll_sig, propagators, raised_power) in enumerate(loop_topo.uv_loop_lines[0]):
                                if raised_power > 0 and all(all(s == 0 for s in param_shift[1]) for _, param_shift, _ in propagators):
                                    if all(all(x == 0 for y in lm_shift for x in y) for s, lm_shift in zip(ll_sig, loop_topo.uv_loop_lines[1]) if s != 0):
                                        d['skip_pf'] = True
                                        break

                        #  construct the normalizing tadpole for the integrated UV counterterm
                        s = uv_subgraph['graph']
                        lm = [s.edge_map_lin[i][0] for i in uv_subgraph['graph'].loop_momenta]
                        g_int = TopologyGenerator([(lm, i, i) for i, lm in enumerate(lm)],
                            powers={lm: 3 for lm in lm})
                        (loop_mom_map, shift_map) = self.topo.build_proto_topology(g_int, c, skip_shift=True)
                        loop_topo = g_int.create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i) + "_ict",
                            # provide dummy external momenta
                            ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                            fixed_deformation=False,
                            mass_map=masses,
                            loop_momentum_map=loop_mom_map,
                            shift_map=shift_map,
                            check_external_momenta_names=False,
                            analytic_result=0)
                        for ll in loop_topo.loop_lines:
                            ll.propagators[0].uv = True
                            ll.propagators[0].m_squared = mu_uv**2
                            ll.propagators[0].power = 3
                            ll.propagators[0].parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]

                        loop_topo.external_kinematics = []

                        uv_subgraph['integrated_ct_bubble_graph'] = loop_topo

                    new_props = []
                    for e, v1, v2 in uv_structure['remaining_graph'].edge_map_lin:
                        lmb_sig = edge_map[e]
                        cmb_sig = [int(x) for x in cb_to_lmb_matrix.T @ Matrix(lmb_sig[0])]
                        assert(all(s == 0 for s in cmb_sig[len(c) - 1:lmb_offset]) and all(s == 0 for s in cmb_sig[lmb_offset + amp_loops:]))

                        amp_sig = numpy.array(cmb_sig[lmb_offset:lmb_offset + amp_loops], dtype=int)

                        if all(s == 0 for s in amp_sig):
                            continue
                        
                        ext_shift_in_cmb = [numpy.array(list(cmb_sig[:lmb_offset]) + [0]*amp_loops + list(cmb_sig[lmb_offset + amp_loops:]), dtype=int),
                                sum(-numpy.array(cc['signature'][1], dtype=int) * s for s, cc in zip(cmb_sig, c[:-1])) + numpy.array(lmb_sig[1], dtype=int)]

                        # simplify ext_shift
                        ext_shift_in_cmb[1][:len(self.topo.ext) // 2] += ext_shift_in_cmb[1][len(self.topo.ext) // 2:]
                        ext_shift_in_cmb[1][len(self.topo.ext) // 2:] = [0]*(len(self.topo.ext) // 2)

                        new_prop = {
                            'name': e,
                            'id': -1,
                            'lmb_sig': tuple(edge_map[e]),
                            'amp_sig': amp_sig.tolist(),
                            'amp_shift_in_cmb_sig': tuple(x.tolist() for x in ext_shift_in_cmb),
                            'mass': masses[e],
                            'uv_mass': False,
                            'uv_only': False,
                        }

                        new_props.append(new_prop)
                        for p in diag_info['propagators']:
                            # TODO: check for an overall sign too
                            if (p['lmb_sig'], p['amp_sig'], p['amp_shift_in_cmb_sig'], p['mass'], p['uv_mass']) == \
                                (new_prop['lmb_sig'], new_prop['amp_sig'], new_prop['amp_shift_in_cmb_sig'], new_prop['mass'], new_prop['uv_mass']):
                                p['uv_only'] = False # also appears as normal propagator
                                new_prop['id'] = p['id']
                                break
                        else:
                            new_prop['id'] = len(diag_info['propagators'])
                            diag_info['propagators'].append(new_prop)

                    # construct all thresholds
                    thresholds = uv_structure['remaining_graph'].find_thresholds(fuse_repeated_edges=False)
                    for ii, t in enumerate(thresholds):
                        foci = [next(p['id'] for p in new_props if p['name'] == f) for f in t[0]]

                        lm_shift = numpy.array([0]*self.topo.n_loops, dtype=int)
                        ext_shift = numpy.array([0]*len(self.topo.ext), dtype=int)
                        for mom, sign in t[1]:
                            sig = edge_map[mom]
                            lm_shift += sign * numpy.array(sig[0], dtype=int)
                            ext_shift += sign * numpy.array(sig[1], dtype=int)

                        # simplify ext_shift
                        ext_shift[:len(self.topo.ext) // 2] += ext_shift[len(self.topo.ext) // 2:]
                        ext_shift[len(self.topo.ext) // 2:] = [0]*(len(self.topo.ext) // 2)

                        # create a map into the focus basis
                        m = Matrix([diag_info['propagators'][d]['amp_sig'] for d in foci[:len(foci) - 1]])
                        dm = m.T * (m * m.T)**-1 # the right inverse should always work
                        # dependent edge in focus basis, should be integer
                        sig_in_fb = Matrix(diag_info['propagators'][foci[-1]]['amp_sig']).T * dm

                        for sign in (1,-1):
                            threshold = {'foci': foci, 'fb_to_cmb': [float(x) for x in dm], 'sig_in_fb': [int(x) for x in sig_in_fb], 'shift_in_lmb_sig': ((sign * lm_shift).tolist(), (sign * ext_shift).tolist()), 'mass_shift': 0.}
                            if threshold not in diag_info['thresholds']:
                                diag_info['thresholds'].append(threshold)

                    forest_to_lmb.extend([x[0] for x in uv_structure['remaining_graph_loop_topo'].loop_momentum_map])

                    if forest_to_lmb != []:
                        fmb_in_cb = (Matrix(forest_to_lmb) * cb_to_lmb_matrix).tolist()
                        shift = Matrix([r[:len(c) - 1] for r in fmb_in_cb])
                        loops = [r[len(c) - 1:] for r in fmb_in_cb]
                        # filter out all columns with zeroes as these are cmb momenta belonging to another amplitude
                        loops = [[a for i, a in enumerate(r) if any(b[i] != 0 for b in loops) ] for r in loops]
                        loopsi = Matrix(loops)**-1
                        shifti = loopsi * shift
                        # we use that UV subgraph defining momenta do not have a shift wrt external momenta
                        remext  = [x[1] for x in uv_structure['remaining_graph_loop_topo'].loop_momentum_map]
                        extshift = Matrix([[0]*len(self.topo.ext)]*(len(loops)-len(remext)) + remext)
                        extshiti = loopsi * extshift

                        uv_structure['forest_to_cb_matrix'] = (loopsi.tolist(), shifti.tolist(), extshiti.tolist(), loops, shift.tolist(), extshift.tolist())
                    else:
                        uv_structure['forest_to_cb_matrix'] = ([[]], [[]], [[]], [[]], [[]], [[]])

                    if uv_structure['remaining_graph'].n_loops == 0 and len(uv_structure['uv_subgraphs']) > 0:
                        pure_forest_counter += 1

                
                for i, t in enumerate(diag_info['thresholds']):
                    t['id'] = i

                lmb_offset += amp_loops
            self.cuts.append(cut_info)

    def export(self, output_path, model=None, include_integration_channel_info=False, optimize_channels=False):
        
        out = {
            'name': self.name,
            'n_loops': self.topo.n_loops,
            'overall_numerator': self.overall_numerator,
            'n_incoming_momenta': len(self.incoming_momenta),
            'external_momenta': [self.external_momenta["q%d"%n] for n in sorted([int(qi.replace("q","")) for qi in self.external_momenta.keys()])],
            'default_fixed_cut_momenta': [[], []] if self.default_kinematics is None else self.default_kinematics,
            'topo': self.loop_topo.to_flat_format(),
            'topo_edges' : [ list(e)+[ (self.topo.powers[e[0]] if i not in self.topo.ext else 0), ]
                                for i, e in enumerate(self.topo.edge_map_lin) ],  # TODO: deprecate
            'edge_PDGs' : [[k,v] for k,v in self.particle_ids.items()],  # TODO: deprecate
            'edge_signatures' : self.topo.get_signature_map(), # TODO: deprecate
            'propagators': [
                {
                    'name': name,
                    'PDG': self.particle_ids[name] if name in self.particle_ids else 0,
                    'signature': self.topo.get_signature_map()[name],
                    'm_squared': self.masses[name]**2 if name in self.masses else 0.,
                    'vertices': (v1, v2),
                    'power': self.topo.powers[name] if name in self.topo.powers else 1,
                }
                for i, (name, v1, v2) in enumerate(self.topo.edge_map_lin) if i not in self.topo.ext
            ],
            'collinear_surfaces': self.sg_collinear,
            'FORM_integrand': self.FORM_integrand,
            # UNCOMMENT the entry below in order to output the information necessary for handling self-energies.
            #'subgraphs_info' : self.subgraphs_info,
            'loop_momentum_basis': [self.topo.edge_map_lin[e][0] for e in self.topo.loop_momenta],
            'e_cm_squared': sum(self.external_momenta[e][0] for e in self.incoming_momenta)**2 - sum(x*x for x in (sum(self.external_momenta[e][i] for e in self.incoming_momenta) for i in range(1, 4))),
            'cutkosky_cuts': [
                {
                    'cuts': [
                        {
                        'name': cut_edge['edge'],
                        'sign': cut_edge['sign'],
                        'power': cut_edge['power'],
                        'particle_id': cut_edge['particle_id'],
                        'signature': cut_edge['signature'],
                        'm_squared': self.masses[cut_edge['edge']]**2 if cut_edge['edge'] in self.masses else 0.,
                        }
                        for cut_edge in cuts['cuts']
                    ],
                    'id': cuts['id'],
                    'amplitudes': [{
                        'graph': diag['uv'][0]['remaining_graph_loop_topo'].to_flat_format(),
                        'conjugate_deformation': diag['conjugate_deformation'],
                        'n_loops': diag['n_loops'],
                        'propagators': diag['propagators'],
                        'thresholds': diag['thresholds'],
                    } for diag in cuts['diagram_info']],
                    'lmb_to_cb': cuts['lmb_to_cb'],
                    'cb_to_lmb': cuts['cb_to_lmb']     
                }

                for cuts in self.cuts
            ]
        }

        multi_channeling_lmb_bases = []
        channel_id = 0
        for lmb_channel in self.topo.loop_momentum_bases():
            defining_edges_for_this_channel = [ self.topo.edge_map_lin[e_id][0] for e_id in lmb_channel ]
            defining_edge_masses = [ model.get_particle(self.particle_ids[e_name]).get('mass') for e_name in defining_edges_for_this_channel ]
            defining_edge_masses = [ 0.0 if mass_param.upper()=='ZERO' else model['parameter_dict'][mass_param].real for mass_param in defining_edge_masses ]
            signatures = [ out['edge_signatures'][edge_name] for edge_name in defining_edges_for_this_channel ]
            multi_channeling_lmb_bases.append(
                {
                    'channel_id' : channel_id,
                    'cutkosky_cut_id' : -1, 
                    'defining_propagator_masses' : defining_edge_masses,
                    'defining_propagators' : defining_edges_for_this_channel,
                    'signatures' : signatures
                }
            )
            channel_id += 1
        out['multi_channeling_lmb_bases'] = multi_channeling_lmb_bases

        optimal_channel_ids = []
        if include_integration_channel_info:

            sig_map=self.topo.get_signature_map()
            sg_edge_powers = { e_data[0]: list(sig_map.values()).count(sig_map[e_data[0]]) for e_data in self.topo.edge_map_lin } 

            from alpha_loop.run_interface import SuperGraph
            sg = SuperGraph(copy.deepcopy(out))

            #print("Now handling SG %s"%sg['name'])
            sg.set_integration_channels()
            multi_channeling_bases = []
            channel_id = 0
            all_channels_added_so_far = []
            for channel in sg['SG_multichannel_info']:
                # Ignore left-side vs right-side cut and only take the left-side representative
                if channel['side'] != 'left':
                    continue

                edge_names_in_basis = []

                # First add the momenta from the Cutkosky cut. Prefer to select vector bosons as independent momenta
                identified_dependent_edge = False
                for cut_edge in sg['cutkosky_cuts'][channel['cutkosky_cut_id']]['cuts']:
                    if not identified_dependent_edge and abs(cut_edge['particle_id']) not in [21,22,23,24]:
                        identified_dependent_edge = True
                        continue
                    edge_names_in_basis.append(cut_edge['name'])
                if not identified_dependent_edge:
                    del edge_names_in_basis[-1]

                # Now complement the basis with all possible ones from the remaining loops to the left and right of the cut
                for LMB in channel['loop_LMBs']:
                    defining_edges_for_this_channel = edge_names_in_basis+LMB['loop_edges']
                    defining_edge_masses = [ model.get_particle(self.particle_ids[e_name]).get('mass') for e_name in defining_edges_for_this_channel ]
                    defining_edge_masses = [ 0.0 if mass_param.upper()=='ZERO' else model['parameter_dict'][mass_param].real for mass_param in defining_edge_masses ]
                    signatures = [ sg['edge_signatures'][edge_name] for edge_name in defining_edges_for_this_channel ]
                    canonical_signatures = set([ (tuple(sig[0]),tuple(sig[1])) for sig in signatures] )
                    if set(canonical_signatures) in all_channels_added_so_far:
                        continue
                    all_channels_added_so_far.append(set(canonical_signatures))
                    multi_channeling_bases.append(
                        {
                            'channel_id' : channel_id,
                            'cutkosky_cut_id' : channel['cutkosky_cut_id'],
                            'defining_propagator_masses' : defining_edge_masses,
                            'defining_propagators' : defining_edges_for_this_channel,
                            'signatures' : signatures,
                            'defining_propagators_powers' : [sg_edge_powers[de_name] for de_name in defining_edges_for_this_channel],
                        }
                    )
                    channel_id += 1

            if optimize_channels:

                all_cc_cuts_edges = [ set([cut_edge['edge'] for cut_edge in cuts['cuts']]) for cuts in self.cuts ]

                channels_score = []
                for i_channel, channel in enumerate(multi_channeling_bases):
                    
                    pdg_score = 0
                    for edge in channel['defining_propagators']:
                        if model is None:
                            if self.particle_ids[edge] in [21,22]:
                                pdg_score +=1
                        else:
                            particle = model.get_particle(self.particle_ids[edge])
                            if particle.get('spin') == 3 and particle.get('mass').upper()=='ZERO':
                                pdg_score += 1
                    
                    power_score = 0
                    for edge in channel['defining_propagators']:
                        power_score += (sg_edge_powers.get(edge,1)-1)
                    
                    cc_score = 0
                    for cut_edges in all_cc_cuts_edges:
                        cc_score -= len( cut_edges.difference(set(channel['defining_propagators'])) )

                    channels_score.append( (pdg_score, power_score, cc_score, i_channel) )
                    
                channels_score.sort( key=lambda chan:(chan[0],chan[1],chan[2]), reverse=True )
                optimal_channel_ids = [ multi_channeling_bases[channels_score[0][-1]]['channel_id'], ]

                # selected_channel_id = channels_score[0][-1]
                # multi_channeling_bases[selected_channel_id]['channel_id'] = 0
                # multi_channeling_bases = [ multi_channeling_bases[selected_channel_id], ]

        else:

            multi_channeling_bases = [
                { 
                    'channel_id' : 0, 
                    'cutkosky_cut_id' : -1, 
                    'defining_propagators' : out['loop_momentum_basis'],
                    'signatures' : [ out['edge_signatures'][edge_name] for edge_name in out['loop_momentum_basis'] ] 
                }
            ]

        out['multi_channeling_bases'] = multi_channeling_bases
        out['optimal_channel_ids'] = optimal_channel_ids

        try:
            import yaml
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        class NoAliasDumper(yaml.SafeDumper):
            def ignore_aliases(self, data):
                return True

        open(output_path,'w').write(yaml.dump(out, Dumper=NoAliasDumper, default_flow_style=None))

if __name__ == "__main__":
    pdgs = {
        'd': 1,
        'u': 2,
        'c': 3,
        's': 4,
        'b': 5,
        't': 6,
        'd~': -1,
        'u~': -2,
        'c~': -3,
        's~': -4,
        'b~': -5,
        't~': -6,
        'g': 21,
        'a': 22,
        'e-': 11,
        'e+': -11,
        'H': 25,
    }

    #ee_to_dd_3l_nested = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 8, 103), ('q4', 8, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 5),
    #('p5', 5, 6), ('p6', 6, 3), ('p7', 6, 7), ('p8', 7, 8), ('p9', 7, 2), ('p10', 5, 4)],
    #    "ee_to_dd_3l_bubble", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
    #    loop_momenta_names=('p2', 'p3', 'p4'),
    #    particle_ids={'p%s' % i: i for i in range(11)},
    #    overall_numerator=1.0)
    #ee_to_dd_3l_nested.export('ee_to_dd_3l_nested.yaml')

    #ee_to_dd_3l_two_bubble_one_line = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 8, 103), ('q4', 8, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 3),
    #('p5', 4, 5), ('p6', 5, 6), ('p7', 6, 5), ('p8', 6, 7), ('p9', 7, 8), ('p10', 7, 2)],
    #    "ee_to_dd_3l_bubble", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
    #    loop_momenta_names=('p2', 'p3', 'p6'),
    #    particle_ids={'p%s' % i: i for i in range(11)},
    #    overall_numerator=1.0)

    # Construct a cross section
    # result is -2 Zeta[3] 3 Pi/(16 Pi^2)^3 = -5.75396*10^-6
    mercedes = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 6),
                                        ('p4', 6, 5), ('p5', 5, 1), ('p6', 2, 4), ('p7', 3, 4), ('p8', 4, 5), ('q2', 6, 7)], "M", ['q1'], 2,
                                        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
                                        loop_momenta_names=('p1', 'p2', 'p3'))
    mercedes.export('mercedes_squared.yaml')

    # result is -5 Zeta[5] 4 Pi/(16 Pi^2)^4 = -1.04773*10^-7
    doublemercedes = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 7), ('p3', 7, 3), ('p4', 3, 6),
                                        ('p5', 6, 5), ('p6', 5, 1), ('p7', 2, 4), ('p8', 3, 4), ('p9', 4, 5), ('p10', 7, 4), ('q2', 6, 8)], "DM", ['q1'], 2,
                                        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]},
                                        loop_momenta_names=('p1', 'p2', 'p3', 'p4'))
    doublemercedes.export('doublemercedes_squared.yaml')

    bubble = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 2), ('q2', 2, 3)], "B", ['q1'], 0,
    {'q1': [2., 0., 0., 0.], 'q2': [2., 0., 0., 0.]},
    masses={'p1': 0.24, 'p2': 0.24})
    bubble.export('bubble_squared.yaml')

    t1 = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 4, 3), ('p4', 4, 1), ('p5', 2, 4), ('q2', 3, 5)], "T", ['q1'], 2,
        {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]}
        #masses={'p1': 0.1, 'p2': 0.1, 'p3': 0.1, 'p4': 0.1,})
        )
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

    # NOTE: for 2 -> N, the first two entries need to be the two incoming momenta
    # the outgoing momenta will be set to the input momenta in the same order, i.e., q3=q1, q4=q2.
    tth = SquaredTopologyGenerator([('q1', 0, 1), ('q2', 6, 7), ('q3', 4, 5), ('q4', 10, 11), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 10),
        ('p5', 10, 9), ('p6', 9, 8), ('p7', 8, 7), ('p8', 1, 7), ('p9', 2, 8), ('p10', 3, 9), ], "TTH", ['q1', 'q2'], 0,
        {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        final_state_particle_ids=(pdgs['t'], pdgs['t~'], pdgs['H']), particle_ids={'p1': pdgs['t'], 'p2': pdgs['t'], 'p3': pdgs['t'], 'p4': pdgs['t~'], 'p5': pdgs['t~'], 'p6': pdgs['t~'], 'p7': pdgs['t~'],
            'p8': pdgs['t'], 'p9': pdgs['g'], 'p10': pdgs['H']})
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
        masses={'p1': 1.0, 'p2': 1.0, 'p3': 1.0, 'p4': 1.0, 'p6': 1.0})
    two_to_three.export('two_to_three_squared.yaml')
