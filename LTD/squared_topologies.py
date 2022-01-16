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

class SquaredTopologyGenerator:
    def __init__(self, edges, name, incoming_momentum_names, n_jets, external_momenta, final_state_particle_ids=(),
        loop_momenta_names=None, loop_momenta_signs=None, masses={}, powers=None, particle_ids={}, jet_ids=None,
        MG_numerator={}, subgraphs_info={},overall_numerator=1., numerator_structure={},
        cut_filter=set(), FORM_numerator={}, FORM_integrand={},
        vertex_weights={}, edge_weights={}, generation_options={},analytic_result=None,
        default_kinematics=None):

        self.name = name
        self.topo = TopologyGenerator(edges, powers)
        self.topo.generate_momentum_flow(loop_momenta_names)
        self.external_momenta = external_momenta
        self.MG_numerator = MG_numerator
        self.FORM_numerator = FORM_numerator
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
            numerator_tensor_coefficients=[[0., 0.]],
            shift_map=None,
            analytic_result=analytic_result
            )

        cutkosky_cuts = self.topo.find_cutkosky_cuts(n_jets, incoming_momentum_names, final_state_particle_ids, 
                                    particle_ids, PDGs_in_jet=jet_ids, fuse_repeated_edges=True)

        if len(cut_filter) > 0:
            self.cuts = [c for c in self.cuts if tuple(n['edge'] for n in c['cuts']) in cut_filter]

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
                        'power': cut_powers[c[0]] }
                    for c in cutkosky_cut],
                'diagram_sets': [{
                        'diagram_info': [{
                            'graph':  g,
                            'conjugate_deformation': g.conjugate,
                        } for g in graphs]
                    }
                ]
            }

            cut_infos.append(cut_info)

        has_2l_bubble = any(len(cut_info['diagram_sets'][0]['diagram_info']) > 2 and
                            any(len(g['graph'].ext) == 2 and len(g['graph'].edge_map_lin) > 4 for g in cut_info['diagram_sets'][0]['diagram_info'])
                        for cut_info in cut_infos)

        for cutkosky_cut, cut_info in zip(cutkosky_cuts, cut_infos):
            c = cut_info['cuts']

            # determine the signature of the cuts
            for cut_edge in c:
                cut_edge['signature'] = copy.deepcopy(edge_map[cut_edge['edge']])
                cut_edge['particle_id'] = self.particle_ids[cut_edge['edge']] if cut_edge['edge'] in self.particle_ids else 0

            cut_name = tuple(a['edge'] for a in c)

            for diag_set in cut_info['diagram_sets']:
                # add the uv structure to the diagram
                for i, diag_info in enumerate(diag_set['diagram_info']):
                    diag_info['graph'].inherit_loop_momentum_basis(self.topo)

                    vw = {v: vertex_weights[v] if v in vertex_weights else 0 for e in diag_info['graph'].edge_map_lin for v in e[1:]}
                    ew = {e: edge_weights[e] if e in edge_weights else -2 for e, _, _ in diag_info['graph'].edge_map_lin}
                    
                    uv_limits = diag_info['graph'].construct_uv_limits(vw, ew, 
                                UV_min_dod_to_subtract=self.generation_options.get('UV_min_dod_to_subtract',0) )

                    internal_bubbles = []
                    for uv_limit in uv_limits:
                        # generate bubble derivative graphs out of remaining graphs
                        bubble_derivatives = []

                        uv_limit['bubble'] = []

                        # find all internal bubble momenta, we use that in gauge theories an internal bubble will be UV divergent
                        for uv_sg in uv_limit['uv_subgraphs']:
                            subgraph_internal_edges = [e[0] for i, e in enumerate(uv_sg['graph'].edge_map_lin) if i not in uv_sg['graph'].ext]

                            uv_sg['internal_bubble'] = None
                            uv_sg['mass_ct'] = None

                            if len(uv_sg['graph'].edge_map_lin) - len(subgraph_internal_edges) == 2:
                                if masses[uv_sg['graph'].edge_map_lin[uv_sg['graph'].ext[0]][0]] > 0.:
                                    uv_sg['mass_ct'] = True

                    diag_info['internal_bubbles'] = internal_bubbles
            
                    # give every subdiagram a globally unique id
                    for uv_limit in uv_limits:
                        for uv_sg in uv_limit['uv_subgraphs']:
                            uv_sg['id'] = graph_counter
                            graph_counter += 1

                            if uv_sg['internal_bubble'] is not None:
                                uv_sg['internal_bubble_id'] = graph_counter
                                graph_counter += 1

                            for dgi, dg in enumerate(uv_sg['derived_graphs']):
                                dg['id'] = graph_counter
                                graph_counter += 1

                                if dg['derivatives'] < uv_sg['taylor_order']:
                                    dg['soft_ct_id'] = graph_counter
                                    graph_counter += 1

                            uv_sg['integrated_ct_id'] = graph_counter
                            graph_counter += 1
                        uv_limit['remaining_graph_id'] = graph_counter
                        graph_counter += 1

                        for b in uv_limit['bubble']:
                            b['id'] = graph_counter
                            graph_counter += 1

                    diag_info.pop('graph')
                    diag_info['uv'] = [{
                        'uv_subgraphs': uv_limit['uv_subgraphs'],
                        'uv_spinney': [[list(g), dod] for g, dod in uv_limit['spinney']],
                        'uv_vertices': [x for x in uv_limit['uv_vertices']],
                        'uv_propagators': [m for g, _ in uv_limit['spinney'] for m in g],
                        'remaining_graph': uv_limit['remaining_graph'],
                        'remaining_graph_id' : uv_limit['remaining_graph_id'],
                        'bubble': uv_limit['bubble']
                    } for uv_limit in uv_limits]

                diag_set['id'] = diagram_set_counter
                diagram_set_counter += 1

                # construct a matrix from the cut basis to the loop momentum basis
                # this is useful if the numerator is specified in the loop momentum basis
                # the matrix will be padded with the loop momentum maps
                cut_to_lmb = [ cut_edge['signature'][0] for cut_edge in c[:-1]]

                for i, diag_info in enumerate(diag_set['diagram_info']):
                    for uv_structure in diag_info['uv']:
                        # create the loop topo of the remaing graph
                        (loop_mom_map, shift_map) = self.topo.build_proto_topology(uv_structure['remaining_graph'], c, skip_shift=False)

                        if uv_structure['uv_spinney'] == []:
                            cut_to_lmb.extend([x[0] for x in loop_mom_map])

                        for bubble in uv_structure['bubble']:
                            # create the loop topo of the derivative bubble graphs
                            (loop_mom_map, shift_map) = self.topo.build_proto_topology(bubble['graph'], c, skip_shift=False)

                            bubble['remaining_graph_loop_topo'] = bubble['graph'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                                # provide dummy external momenta
                                ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                                fixed_deformation=False,
                                mass_map=masses,
                                loop_momentum_map=loop_mom_map,
                                numerator_tensor_coefficients=[[0., 0.,]],
                                shift_map=shift_map,
                                check_external_momenta_names=False,
                                analytic_result=0)
                            bubble['remaining_graph_loop_topo'].external_kinematics = []

                        uv_structure['remaining_graph_loop_topo'] = uv_structure['remaining_graph'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                            # provide dummy external momenta
                            ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                            fixed_deformation=False,
                            mass_map=masses,
                            loop_momentum_map=loop_mom_map,
                            numerator_tensor_coefficients=[[0., 0.,]],
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
                lmb_to_cb_matrix = lmb_to_cb_matrix**-1

                diag_set['cb_to_lmb'] = [int(x) for x in lmb_to_cb_matrix]

                # we know that the loop momenta in the cmb are direcly one of the loop momenta in the lmb
                used_loop_momenta = [i for i in range(self.topo.n_loops) if
                    len([cc for cc in diag_set['cb_to_lmb'][i * self.topo.n_loops:(i+1) * self.topo.n_loops] if cc != 0]) == 1 and
                    any(cc != 0 for cc in diag_set['cb_to_lmb'][i * self.topo.n_loops + len(c) - 1:(i+1) * self.topo.n_loops])]
                assert(len(used_loop_momenta) == self.topo.n_loops - len(c) + 1)

                # construct the forest matrix that maps the amplitude momenta in the cmb
                # to ones suitable for the spinney
                for i, diag_info in enumerate(diag_set['diagram_info']):
                    for uv_structure in diag_info['uv']:
                        # create the LTD representation of the derived UV graph
                        forest_to_cb = []
                        for uv_subgraph in uv_structure['uv_subgraphs']:
                            for di, d in enumerate(uv_subgraph['derived_graphs']):
                                (loop_mom_map, shift_map) = self.topo.build_proto_topology(d['graph'], c, skip_shift=True)

                                # Construct a loop momentum map for the loop momenta that remain after removing
                                # all dependence on external momenta in the cmb.
                                # ie for fmb defining edge c3-c4+c1+p1+p2, where c1 is cut we construct
                                # f1 = c3 - c4, f1_shift = c1+p1+p2
                                # The shift should be added later to the parametric shifts of each propagator containing f1.
                                basis_shift_map = []
                                new_lm_map = []
                                for lmp, shp in loop_mom_map:
                                    # get the dependence of the loop momenta
                                    dep = numpy.array([0]*self.topo.n_loops, dtype=int)
                                    for ii, x in enumerate(lmp):
                                        dep += x* numpy.array(diag_set['cb_to_lmb'][ii * self.topo.n_loops:(ii+1) * self.topo.n_loops], dtype=int)
                                    deps = [0 if u not in used_loop_momenta else
                                            dep.dot(numpy.array(diag_set['cb_to_lmb'][u * self.topo.n_loops:(u+1) * self.topo.n_loops], dtype=int))
                                            for u in range(self.topo.n_loops)]

                                    basis_shift_map.append([[a-b for a,b in zip(lmp,deps)],shp])
                                    new_lm_map.append([deps, [0]*len(shp)])

                                loop_mom_map = new_lm_map

                                if di == 0:
                                    forest_to_cb.extend([x[0] for x in loop_mom_map])

                                # note: the shift map signs may get swapped when edges switch orientation
                                # the basis_shift_map should be unaffected since defining edges are not swapped
                                loop_topo = d['graph'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                                    # provide dummy external momenta
                                    ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                                    fixed_deformation=False,
                                    mass_map=masses,
                                    loop_momentum_map=loop_mom_map,
                                    numerator_tensor_coefficients=[[0., 0.,]],
                                    shift_map=shift_map,
                                    check_external_momenta_names=False,
                                    analytic_result=0)
                                loop_topo.external_kinematics = []

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

                                    # note: we assume that every propagator has power 1 and that only merging of identical edges raises it
                                    uv_loop_lines.append((ll.signature, [(p.name, p.parametric_shift,
                                        1 + len([1 for _,mn in loop_topo.propagator_map.items() if mn == p.name])) for p in ll.propagators], derived_ll_power - orig_ll_power))
                                    prop = ll.propagators[0]

                                    # keep the original mass for the soft CT part of the Taylor expansion
                                    if d['derivatives'] >= uv_subgraph['taylor_order']:
                                        prop.uv = True
                                        prop.m_squared = mu_uv**2
                                    prop.power = derived_ll_power
                                    prop.parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]
                                    ll.propagators = [prop]

                                loop_topo.uv_loop_lines = (uv_loop_lines, basis_shift_map)

                                if d['derivatives'] < uv_subgraph['taylor_order']:
                                    d['loop_topo_orig_mass'] = copy.deepcopy(loop_topo)
                                    for ll in loop_topo.loop_lines:
                                        ll.propagators[0].uv = True
                                        ll.propagators[0].m_squared = mu_uv**2

                                d['loop_topo'] = loop_topo

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
                                numerator_tensor_coefficients=[[0., 0.,]],
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

                            if uv_subgraph['internal_bubble']:
                                # keep the shift in terms of the local topology
                                (loop_mom_map, shift_map) = self.topo.build_proto_topology(uv_subgraph['internal_bubble'], c, skip_shift=True)

                                uv_subgraph['internal_bubble_loop_topo'] = uv_subgraph['internal_bubble'].create_loop_topology(name + '_' + ''.join(cut_name) + '_' + str(i),
                                    # provide dummy external momenta
                                    ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                                    fixed_deformation=False,
                                    mass_map=masses,
                                    loop_momentum_map=loop_mom_map,
                                    numerator_tensor_coefficients=[[0., 0.,]],
                                    shift_map=shift_map,
                                    check_external_momenta_names=False,
                                    analytic_result=0)
                                uv_subgraph['internal_bubble_loop_topo'].external_kinematics = []

                        forest_to_cb.extend([x[0] for x in uv_structure['remaining_graph_loop_topo'].loop_momentum_map])

                        if forest_to_cb != []:
                            fmb_in_cb = (Matrix(forest_to_cb) * lmb_to_cb_matrix).tolist()
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
                                for i, e in enumerate(self.topo.edge_map_lin) ],
            'edge_PDGs' : [[k,v] for k,v in self.particle_ids.items()],
            'edge_signatures' : self.topo.get_signature_map(),
            'MG_numerator': self.MG_numerator,
            'FORM_numerator': self.FORM_numerator,
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
                    'diagram_sets': [
                        {
                            'id': diag_set['id'],
                            'diagram_info': [{
                                'graph': diag['uv'][0]['remaining_graph_loop_topo'].to_flat_format(),
                                'conjugate_deformation': diag['conjugate_deformation'],
                            } for diag in diag_set['diagram_info']],
                            'cb_to_lmb': diag_set['cb_to_lmb']
                        }
                    for diag_set in cuts['diagram_sets']]
                }

                for cuts in self.cuts
            ]
        }

        multi_channeling_lmb_bases = []
        channel_id = 0
        for lmb_channel in self.topo.loop_momentum_bases():
            defining_edges_for_this_channel = [ self.topo.edge_map_lin[e_id][0] for e_id in lmb_channel ]
            signatures = [ out['edge_signatures'][edge_name] for edge_name in defining_edges_for_this_channel ]
            multi_channeling_lmb_bases.append(
                {
                    'channel_id' : channel_id,
                    'cutkosky_cut_id' : -1, 
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
                    signatures = [ sg['edge_signatures'][edge_name] for edge_name in defining_edges_for_this_channel ]
                    canonical_signatures = set([ (tuple(sig[0]),tuple(sig[1])) for sig in signatures] )
                    if set(canonical_signatures) in all_channels_added_so_far:
                        continue
                    all_channels_added_so_far.append(set(canonical_signatures))
                    multi_channeling_bases.append(
                        {
                            'channel_id' : channel_id,
                            'cutkosky_cut_id' : channel['cutkosky_cut_id'], 
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
