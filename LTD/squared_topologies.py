import os
import sys
root_path = os.path.dirname(os.path.realpath( __file__ ))
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
                                    particle_ids, PDGs_in_jet=jet_ids)

        self.cuts = [self.topo.bubble_cuts(c, incoming_momentum_names) for c in cutkosky_cuts]

        if len(cut_filter) > 0:
            self.cuts = [c for c in self.cuts if tuple(n['edge'] for n in c['cuts']) in cut_filter]

        self.masses = copy.deepcopy(masses)
        self.overall_numerator = overall_numerator
        self.incoming_momenta = incoming_momentum_names

        edge_map = self.topo.get_signature_map()
        mu_uv = 2. * math.sqrt(sum(self.external_momenta[e][0] for e in self.incoming_momenta)**2 - 
                        sum(x*x for x in (sum(self.external_momenta[e][i] for e in self.incoming_momenta) for i in range(1, 4))))

        #self.cut_diagrams = []
        diagram_set_counter = 0
        graph_counter = 0
        diagram_sets = []
        for cut_info in self.cuts:
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
                    # correct the edge weight of the bubble derivative
                    # the derivative wrt the numerator does not change the UV scaling
                    if diag_info['derivative'] and diag_info['derivative'][0] != diag_info['derivative'][1]:
                        ew[diag_info['derivative'][1]] -= 1
                    
                    uv_limits = diag_info['graph'].construct_uv_limits(vw, ew, 
                                UV_min_dod_to_subtract=self.generation_options.get('UV_min_dod_to_subtract',0) )
            
                    # give every subdiagram a globally unique id
                    for uv_limit in uv_limits:
                        for uv_sg in uv_limit['uv_subgraphs']:
                            uv_sg['id'] = graph_counter
                            graph_counter += 1
                            uv_sg['graph_index'] += i * 100
                            uv_sg['subgraph_indices'] = [j + i * 100 for j in uv_sg['subgraph_indices']]

                            for dg in uv_sg['derived_graphs']:
                                dg['id'] = graph_counter
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

                """
                uv_diag_sets_with_integrated_ct = []

                # unpack the factorized UV subgraphs
                for uv_diag_set in uv_diag_sets:
                    unfolded_diag_info = []
                    for di in uv_diag_set['diagram_info']:
                        # only the remaining graph gets the derivative flag to prevent
                        # it being applied more than once per diagram set
                        unfolded_diag_info.append({
                            'uv_info': None,
                            'uv_vertices': di['uv_vertices'],
                            'graph': di['remaining_graph'],
                            'derivative': di['derivative'],
                            'derivative_edge': None,
                            'bubble_momenta': di['bubble_momenta'],
                            'conjugate_deformation': di['conjugate_deformation']
                        })

                        for uv_lim in di['uv_subgraphs']:
                            unfolded_diag_info.append({
                                'uv_info': uv_lim,
                                'uv_vertices': None,
                                'graph': uv_lim['graph'],
                                'derivative': None,
                                'derivative_edge': di['derivative'][1] if di['derivative'] and di['derivative'][0] != di['derivative'][1] else None,
                                'bubble_momenta': [],
                                'conjugate_deformation': di['conjugate_deformation']
                            })

                    # take the cartesian product over all local+integrated CT
                    ct_opts = [[False] if (di['uv_info'] is None or not self.generation_options.get('generate_integrated_UV_CTs',True)) else [False, True] for di in unfolded_diag_info]
                    #print(ct_opts)
                    for o in product(*ct_opts):
                        new_diag_info = []
                        for ictflag, uv_diag_info in zip(o, unfolded_diag_info):
                            integrated_diag_info = copy.deepcopy(uv_diag_info)
                            integrated_diag_info['integrated_ct'] = ictflag
                            new_diag_info.append(integrated_diag_info)

                        uv_diag_sets_with_integrated_ct.append(
                            {
                                'diagram_info': new_diag_info,
                                'uv_spinney': copy.deepcopy(uv_diag_set['uv_spinney']),
                                'uv_propagators': copy.deepcopy(uv_diag_set['uv_propagators']),
                            }
                        )
                """

                    #uv_diag_set['diagram_info'] = unfolded_diag_info

                # add integrated counterterms to the diagram set
                #uv_diag_sets_with_integrated_ct = []
                #for uv_diag_set in uv_diag_sets:
                #    uv_diag_set['integrated_ct'] = False
                #    uv_diag_sets_with_integrated_ct.append(uv_diag_set)
                #    if self.generation_options.get('generate_integrated_UV_CTs',True):
                #        if any(di['uv_info'] is not None for di in uv_diag_set['diagram_info']):
                #            integrated_diag_set = copy.deepcopy(uv_diag_set)
                #            integrated_diag_set['integrated_ct'] = True
                #            uv_diag_sets_with_integrated_ct.append(integrated_diag_set)

                #uv_diagram_sets.extend(uv_diag_sets_with_integrated_ct)

            #cut_info['diagram_sets'] = uv_diagram_sets

                diag_set['id'] = diagram_set_counter
                diagram_set_counter += 1

                # construct a matrix from the cut basis to the loop momentum basis
                # this is useful if the numerator is specified in the loop momentum basis
                # the matrix will be padded with the loop momentum maps
                cut_to_lmb = [ cut_edge['signature'][0] for cut_edge in c[:-1]]

                loop_topos = []
                for i, diag_info in enumerate(diag_set['diagram_info']):
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
                            numerator_tensor_coefficients=[[0., 0.,]],
                            shift_map=shift_map,
                            check_external_momenta_names=False,
                            analytic_result=0)
                        uv_structure['remaining_graph_loop_topo'].external_kinematics = []

                    """
                    pprint(diag_info)
                    s = diag_info['graph']



                    #print('uvl', uv_loop_lines)

                    if diag_info['integrated_ct']:
                        # replace the graphs by finite vacuum bubbles
                        # the numerator will be the integrated CT
                        lm = [s.edge_map_lin[i][0] for i in s.loop_momenta]
                        g_int = TopologyGenerator([(lm, i, i) for i, lm in enumerate(lm)],
                            powers={lm: 3 for lm in lm})
                        (loop_mom_map, shift_map) = self.topo.build_proto_topology(g_int, c, skip_shift=diag_info['uv_info'] is not None)
                        loop_topo = g_int.create_loop_topology(self.name + '_' + ''.join(cut_name) + uv_name + '_' + str(i),
                            # provide dummy external momenta
                            ext_mom={edge_name: vectors.LorentzVector([0, 0, 0, 0]) for (edge_name, _, _) in self.topo.edge_map_lin},
                            fixed_deformation=False,
                            mass_map=masses,
                            loop_momentum_map=loop_mom_map,
                            numerator_tensor_coefficients=[[0., 0.,]],#[[0., 0.] for _ in range(numerator_entries)],
                            shift_map=shift_map,
                            check_external_momenta_names=False,
                            analytic_result=0)
                        for ll in loop_topo.loop_lines:
                            ll.propagators[0].uv = True
                            ll.propagators[0].m_squared = mu_uv**2
                            ll.propagators[0].power = 3
                            ll.propagators[0].parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]

                    loop_topo.external_kinematics = []

                    loop_topos.append(
                        {
                            'graph': loop_topo,
                            'loop_momentum_map': loop_mom_map,
                            'uv_loop_lines': uv_loop_lines,
                            'conjugate_deformation': diag_info['conjugate_deformation']
                        })
                     """

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
                                    orig_ll_power = sum(uv_subgraph['graph'].powers[pp.name] for pp in ll.propagators)

                                    # FIXME: repeated propagators will only appear once in ll.propagators and will therefore not be
                                    # getting the correct power
                                    uv_loop_lines.append((ll.signature, [(p.name, p.parametric_shift) for p in ll.propagators], derived_ll_power - orig_ll_power))
                                    prop = ll.propagators[0]
                                    prop.uv = True
                                    prop.m_squared = mu_uv**2
                                    prop.power = derived_ll_power
                                    prop.parametric_shift = [[0 for _ in c], [0 for _ in range(len(incoming_momentum_names) * 2)]]
                                    ll.propagators = [prop]

                                loop_topo.uv_loop_lines = (uv_loop_lines, basis_shift_map)

                                d['loop_topo'] = loop_topo

                        forest_to_cb.extend([x[0] for x in uv_structure['remaining_graph_loop_topo'].loop_momentum_map])

                        if forest_to_cb != []:
                            fmb_in_cb = (Matrix(forest_to_cb) * lmb_to_cb_matrix).tolist()
                            shift = Matrix([r[:len(c) - 1] for r in fmb_in_cb])
                            loops = [r[len(c) - 1:] for r in fmb_in_cb]
                            # filter out all columns with zeroes as these are cmb momenta belonging to another amplitude
                            loops = [[a for i, a in enumerate(r) if any(b[i] != 0 for b in loops) ] for r in loops]
                            loopsi = Matrix(loops)**-1
                            shifti = loopsi * shift

                            uv_structure['forest_to_cb_matrix'] = (loopsi.tolist() , shifti.tolist(), loops, shift.tolist())
                        else:
                            uv_structure['forest_to_cb_matrix'] = ([[]], [[]], [[]], [[]])

    def export(self, output_path):
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

        try:
            import yaml
            from yaml import Loader, Dumper
        except ImportError:
            raise BaseException("Install yaml python module in order to import topologies from yaml.")

        open(output_path,'w').write(yaml.dump(out, Dumper=Dumper, default_flow_style=None))

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

    ee_to_dd_2l_bubble = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 6, 103), ('q4', 6, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 4), ('p4', 4, 3),
    ('p5', 4, 5), ('p6', 5, 6), ('p7', 5, 2)],
        "ee_to_dd_2l_bubble", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        loop_momenta_names=('p2', 'p3'),
        particle_ids={ 'p1': pdgs['a'], 'p2': pdgs['d'], 'p3': pdgs['d'], 'p4': pdgs['g'], 'p5': pdgs['d'], 'p6': pdgs['a'], 'p7': pdgs['d~'] },
        overall_numerator=1.0,
        #cut_filter={('p2', 'p3', 'p4', 'p7')},
        numerator_structure={('p3', 'p4', 'p7'):
            {():
            [[[0,4,4],[0.,-1.477080880741393e8]],
            [[0,5,5],[0.,-1.477080880741393e8]],
            [[0,6,6],[0.,-1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[0,4,4,4],[0.,-7.385404403706966e7]],
            [[0,4,5,5],[0.,-7.385404403706966e7]],
            [[0,4,6,6],[0.,-7.385404403706966e7]],
            [[0,4,7,7],[0.,+7.385404403706966e7]],
            [[1,4,5],[0.,+2.954161761482786e8]],
            [[1,4,4,5],[0.,+1.477080880741393e8]],
            [[1,5,7,7],[0.,-1.477080880741393e8]],
            [[2,4,6],[0.,+2.954161761482786e8]],
            [[2,4,4,6],[0.,+1.477080880741393e8]],
            [[2,6,7,7],[0.,-1.477080880741393e8]],
            [[3,4,7],[0.,+2.954161761482786e8]],
            [[3,4,4,7],[0.,+7.385404403706966e7]],
            [[3,5,5,7],[0.,+7.385404403706966e7]],
            [[3,6,6,7],[0.,+7.385404403706966e7]],
            [[3,7,7,7],[0.,-7.385404403706966e7]],
            [[0,0,4],[0.,+2.954161761482786e8]],
            [[0,0,4,4],[0.,+2.215621321112089e8]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,-7.385404403706966e7]],
            [[0,1,4,5],[0.,-1.477080880741393e8]],
            [[0,2,4,6],[0.,-1.477080880741393e8]],
            [[0,3,4,7],[0.,-2.954161761482786e8]],
            [[1,1,4],[0.,-2.954161761482786e8]],
            [[1,1,4,4],[0.,-1.477080880741393e8]],
            [[1,1,7,7],[0.,+1.477080880741393e8]],
            [[1,3,5,7],[0.,+1.477080880741393e8]],
            [[2,2,4],[0.,-2.954161761482786e8]],
            [[2,2,4,4],[0.,-1.477080880741393e8]],
            [[2,2,7,7],[0.,+1.477080880741393e8]],
            [[2,3,6,7],[0.,+1.477080880741393e8]],
            [[3,3,4],[0.,-2.954161761482786e8]],
            [[3,3,4,4],[0.,-7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,+2.215621321112089e8]],
            [[0,0,0],[0.,-1.477080880741393e8]],
            [[0,0,0,4],[0.,-2.215621321112089e8]],
            [[0,0,3,7],[0.,+2.215621321112089e8]],
            [[0,1,1],[0.,+1.477080880741393e8]],
            [[0,1,1,4],[0.,+2.215621321112089e8]],
            [[0,2,2],[0.,+1.477080880741393e8]],
            [[0,2,2,4],[0.,+2.215621321112089e8]],
            [[0,3,3],[0.,+1.477080880741393e8]],
            [[0,3,3,4],[0.,+2.215621321112089e8]],
            [[1,1,3,7],[0.,-2.215621321112089e8]],
            [[2,2,3,7],[0.,-2.215621321112089e8]],
            [[3,3,3,7],[0.,-2.215621321112089e8]],
            [[0,0,0,0],[0.,+7.385404403706966e7]],
            [[0,0,1,1],[0.,-7.385404403706966e7]],
            [[0,0,2,2],[0.,-7.385404403706966e7]],
            [[0,0,3,3],[0.,-1.477080880741393e8]],
            [[1,1,3,3],[0.,+7.385404403706966e7]],
            [[2,2,3,3],[0.,+7.385404403706966e7]],
            [[3,3,3,3],[0.,+7.385404403706966e7]]]
            },
            ('p5','p3','p4','p7'):
            {():
            [[[4,8,8,8],[2.350847235165231e7,0.]],
            [[4,8,9,9],[2.350847235165231e7,0.]],
            [[4,8,10,10],[2.350847235165231e7,0.]],
            [[4,8,11,11],[2.350847235165231e7,0.]],
            [[5,8,8,9],[-4.701694470330462e7,0.]],
            [[6,8,8,10],[-4.701694470330462e7,0.]],
            [[7,8,8,11],[-4.701694470330462e7,0.]],
            [[4,4,8,8],[-7.052541705495694e7,0.]],
            [[4,4,9,9],[-2.350847235165231e7,0.]],
            [[4,4,10,10],[-2.350847235165231e7,0.]],
            [[4,4,11,11],[-2.350847235165231e7,0.]],
            [[4,5,8,9],[4.701694470330462e7,0.]],
            [[4,6,8,10],[4.701694470330462e7,0.]],
            [[4,7,8,11],[4.701694470330462e7,0.]],
            [[5,5,8,8],[4.701694470330462e7,0.]],
            [[6,6,8,8],[4.701694470330462e7,0.]],
            [[7,7,8,8],[4.701694470330462e7,0.]],
            [[4,4,4,8],[7.052541705495694e7,0.]],
            [[4,5,5,8],[-7.052541705495694e7,0.]],
            [[4,6,6,8],[-7.052541705495694e7,0.]],
            [[4,7,7,8],[-7.052541705495694e7,0.]],
            [[4,4,4,4],[-2.350847235165231e7,0.]],
            [[4,4,5,5],[2.350847235165231e7,0.]],
            [[4,4,6,6],[2.350847235165231e7,0.]],
            [[4,4,7,7],[2.350847235165231e7,0.]],
            [[0,4,8,8],[-2.350847235165231e7,0.]],
            [[0,4,9,9],[-2.350847235165231e7,0.]],
            [[0,4,10,10],[-2.350847235165231e7,0.]],
            [[0,4,11,11],[-2.350847235165231e7,0.]],
            [[0,4,8,8,8],[-1.175423617582616e7,0.]],
            [[0,4,8,9,9],[-1.175423617582616e7,0.]],
            [[0,4,8,10,10],[-1.175423617582616e7,0.]],
            [[0,4,8,11,11],[-1.175423617582616e7,0.]],
            [[0,5,8,9],[4.701694470330462e7,0.]],
            [[0,5,8,8,9],[2.350847235165232e7,0.]],
            [[0,6,8,10],[4.701694470330462e7,0.]],
            [[0,6,8,8,10],[2.350847235165232e7,0.]],
            [[0,7,8,11],[4.701694470330462e7,0.]],
            [[0,7,8,8,11],[2.350847235165232e7,0.]],
            [[0,4,4,8],[4.701694470330462e7,0.]],
            [[0,4,4,8,8],[3.526270852747848e7,0.]],
            [[0,4,4,9,9],[1.175423617582616e7,0.]],
            [[0,4,4,10,10],[1.175423617582616e7,0.]],
            [[0,4,4,11,11],[1.175423617582616e7,0.]],
            [[0,4,5,8,9],[-2.350847235165232e7,0.]],
            [[0,4,6,8,10],[-2.350847235165232e7,0.]],
            [[0,4,7,8,11],[-2.350847235165232e7,0.]],
            [[0,5,5,8],[-4.701694470330462e7,0.]],
            [[0,5,5,8,8],[-2.350847235165232e7,0.]],
            [[0,6,6,8],[-4.701694470330462e7,0.]],
            [[0,6,6,8,8],[-2.350847235165232e7,0.]],
            [[0,7,7,8],[-4.701694470330462e7,0.]],
            [[0,7,7,8,8],[-2.350847235165232e7,0.]],
            [[0,4,4,4],[-2.350847235165231e7,0.]],
            [[0,4,4,4,8],[-3.526270852747848e7,0.]],
            [[0,4,5,5],[2.350847235165231e7,0.]],
            [[0,4,5,5,8],[3.526270852747848e7,0.]],
            [[0,4,6,6],[2.350847235165231e7,0.]],
            [[0,4,6,6,8],[3.526270852747848e7,0.]],
            [[0,4,7,7],[2.350847235165231e7,0.]],
            [[0,4,7,7,8],[3.526270852747848e7,0.]],
            [[0,4,4,4,4],[1.175423617582616e7,0.]],
            [[0,4,4,5,5],[-1.175423617582616e7,0.]],
            [[0,4,4,6,6],[-1.175423617582616e7,0.]],
            [[0,4,4,7,7],[-1.175423617582616e7,0.]],
            [[3,4,8,8,11],[2.350847235165232e7,0.]],
            [[3,5,8,9,11],[-2.350847235165232e7,0.]],
            [[3,6,8,10,11],[-2.350847235165232e7,0.]],
            [[3,7,8,8,8],[-1.175423617582616e7,0.]],
            [[3,7,8,9,9],[1.175423617582616e7,0.]],
            [[3,7,8,10,10],[1.175423617582616e7,0.]],
            [[3,7,8,11,11],[-1.175423617582616e7,0.]],
            [[3,4,4,8,11],[-4.701694470330463e7,0.]],
            [[3,4,5,9,11],[2.350847235165232e7,0.]],
            [[3,4,6,10,11],[2.350847235165232e7,0.]],
            [[3,4,7,8,8],[1.175423617582616e7,0.]],
            [[3,4,7,9,9],[-1.175423617582616e7,0.]],
            [[3,4,7,10,10],[-1.175423617582616e7,0.]],
            [[3,4,7,11,11],[1.175423617582616e7,0.]],
            [[3,5,5,8,11],[2.350847235165232e7,0.]],
            [[3,6,6,8,11],[2.350847235165232e7,0.]],
            [[3,7,7,8,11],[2.350847235165232e7,0.]],
            [[3,4,4,4,11],[2.350847235165232e7,0.]],
            [[3,4,4,7,8],[1.175423617582616e7,0.]],
            [[3,4,5,5,11],[-2.350847235165232e7,0.]],
            [[3,4,6,6,11],[-2.350847235165232e7,0.]],
            [[3,4,7,7,11],[-2.350847235165232e7,0.]],
            [[3,5,5,7,8],[-1.175423617582616e7,0.]],
            [[3,6,6,7,8],[-1.175423617582616e7,0.]],
            [[3,7,7,7,8],[-1.175423617582616e7,0.]],
            [[3,4,4,4,7],[-1.175423617582616e7,0.]],
            [[3,4,5,5,7],[1.175423617582616e7,0.]],
            [[3,4,6,6,7],[1.175423617582616e7,0.]],
            [[3,4,7,7,7],[1.175423617582616e7,0.]],
            [[0,0,4,8,8],[1.175423617582616e7,0.]],
            [[0,0,4,9,9],[1.175423617582616e7,0.]],
            [[0,0,4,10,10],[1.175423617582616e7,0.]],
            [[0,0,4,11,11],[1.175423617582616e7,0.]],
            [[0,0,5,8,9],[-2.350847235165232e7,0.]],
            [[0,0,6,8,10],[-2.350847235165232e7,0.]],
            [[0,0,7,8,11],[-2.350847235165232e7,0.]],
            [[0,0,4,4,8],[-2.350847235165232e7,0.]],
            [[0,0,5,5,8],[2.350847235165232e7,0.]],
            [[0,0,6,6,8],[2.350847235165232e7,0.]],
            [[0,0,7,7,8],[2.350847235165232e7,0.]],
            [[0,0,4,4,4],[1.175423617582616e7,0.]],
            [[0,0,4,5,5],[-1.175423617582616e7,0.]],
            [[0,0,4,6,6],[-1.175423617582616e7,0.]],
            [[0,0,4,7,7],[-1.175423617582616e7,0.]],
            [[0,3,4,8,11],[-2.350847235165232e7,0.]],
            [[0,3,5,9,11],[2.350847235165232e7,0.]],
            [[0,3,6,10,11],[2.350847235165232e7,0.]],
            [[0,3,7,8,8],[1.175423617582616e7,0.]],
            [[0,3,7,9,9],[-1.175423617582616e7,0.]],
            [[0,3,7,10,10],[-1.175423617582616e7,0.]],
            [[0,3,7,11,11],[1.175423617582616e7,0.]],
            [[0,3,4,4,11],[2.350847235165232e7,0.]],
            [[0,3,5,5,11],[-2.350847235165232e7,0.]],
            [[0,3,6,6,11],[-2.350847235165232e7,0.]],
            [[0,3,7,7,11],[-2.350847235165232e7,0.]],
            [[0,3,4,4,7],[-1.175423617582616e7,0.]],
            [[0,3,5,5,7],[1.175423617582616e7,0.]],
            [[0,3,6,6,7],[1.175423617582616e7,0.]],
            [[0,3,7,7,7],[1.175423617582616e7,0.]]],
            # the UV approximation cannot be generated automatically
            #('p3',):
            #[[[0],[0.,-7.562654109395933e10]],
            #[[0,4,4],[0.,2.363329409186229e9]],
            #[[3,4,4,7],[0.,-1.181664704593115e9]],
            #[[0,0],[0.,3.781327054697966e10]],
            #[[0,0,4,4],[0.,-1.181664704593115e9]],
            #[[3,3],[0.,-3.781327054697966e10]],
            #[[3,3,4,4],[0.,1.181664704593115e9]]]
            },
            ('p2','p3','p4','p7'):
            {():
            [[[4,8,8,8],[-2.350847235165231e7,0.]],
            [[4,8,9,9],[-2.350847235165231e7,0.]],
            [[4,8,10,10],[-2.350847235165231e7,0.]],
            [[4,8,11,11],[-2.350847235165231e7,0.]],
            [[5,8,8,9],[4.701694470330462e7,0.]],
            [[6,8,8,10],[4.701694470330462e7,0.]],
            [[7,8,8,11],[4.701694470330462e7,0.]],
            [[4,4,8,8],[7.052541705495694e7,0.]],
            [[4,4,9,9],[2.350847235165231e7,0.]],
            [[4,4,10,10],[2.350847235165231e7,0.]],
            [[4,4,11,11],[2.350847235165231e7,0.]],
            [[4,5,8,9],[-4.701694470330462e7,0.]],
            [[4,6,8,10],[-4.701694470330462e7,0.]],
            [[4,7,8,11],[-4.701694470330462e7,0.]],
            [[5,5,8,8],[-4.701694470330462e7,0.]],
            [[6,6,8,8],[-4.701694470330462e7,0.]],
            [[7,7,8,8],[-4.701694470330462e7,0.]],
            [[4,4,4,8],[-7.052541705495694e7,0.]],
            [[4,5,5,8],[7.052541705495694e7,0.]],
            [[4,6,6,8],[7.052541705495694e7,0.]],
            [[4,7,7,8],[7.052541705495694e7,0.]],
            [[4,4,4,4],[2.350847235165231e7,0.]],
            [[4,4,5,5],[-2.350847235165231e7,0.]],
            [[4,4,6,6],[-2.350847235165231e7,0.]],
            [[4,4,7,7],[-2.350847235165231e7,0.]],
            [[0,4,8,8],[2.350847235165231e7,0.]],
            [[0,4,9,9],[2.350847235165231e7,0.]],
            [[0,4,10,10],[2.350847235165231e7,0.]],
            [[0,4,11,11],[2.350847235165231e7,0.]],
            [[0,4,8,8,8],[1.175423617582616e7,0.]],
            [[0,4,8,9,9],[1.175423617582616e7,0.]],
            [[0,4,8,10,10],[1.175423617582616e7,0.]],
            [[0,4,8,11,11],[1.175423617582616e7,0.]],
            [[0,5,8,9],[-4.701694470330462e7,0.]],
            [[0,5,8,8,9],[-2.350847235165232e7,0.]],
            [[0,6,8,10],[-4.701694470330462e7,0.]],
            [[0,6,8,8,10],[-2.350847235165232e7,0.]],
            [[0,7,8,11],[-4.701694470330462e7,0.]],
            [[0,7,8,8,11],[-2.350847235165232e7,0.]],
            [[0,4,4,8],[-4.701694470330462e7,0.]],
            [[0,4,4,8,8],[-3.526270852747848e7,0.]],
            [[0,4,4,9,9],[-1.175423617582616e7,0.]],
            [[0,4,4,10,10],[-1.175423617582616e7,0.]],
            [[0,4,4,11,11],[-1.175423617582616e7,0.]],
            [[0,4,5,8,9],[2.350847235165232e7,0.]],
            [[0,4,6,8,10],[2.350847235165232e7,0.]],
            [[0,4,7,8,11],[2.350847235165232e7,0.]],
            [[0,5,5,8],[4.701694470330462e7,0.]],
            [[0,5,5,8,8],[2.350847235165232e7,0.]],
            [[0,6,6,8],[4.701694470330462e7,0.]],
            [[0,6,6,8,8],[2.350847235165232e7,0.]],
            [[0,7,7,8],[4.701694470330462e7,0.]],
            [[0,7,7,8,8],[2.350847235165232e7,0.]],
            [[0,4,4,4],[2.350847235165231e7,0.]],
            [[0,4,4,4,8],[3.526270852747848e7,0.]],
            [[0,4,5,5],[-2.350847235165231e7,0.]],
            [[0,4,5,5,8],[-3.526270852747848e7,0.]],
            [[0,4,6,6],[-2.350847235165231e7,0.]],
            [[0,4,6,6,8],[-3.526270852747848e7,0.]],
            [[0,4,7,7],[-2.350847235165231e7,0.]],
            [[0,4,7,7,8],[-3.526270852747848e7,0.]],
            [[0,4,4,4,4],[-1.175423617582616e7,0.]],
            [[0,4,4,5,5],[1.175423617582616e7,0.]],
            [[0,4,4,6,6],[1.175423617582616e7,0.]],
            [[0,4,4,7,7],[1.175423617582616e7,0.]],
            [[3,4,8,8,11],[-2.350847235165232e7,0.]],
            [[3,5,8,9,11],[2.350847235165232e7,0.]],
            [[3,6,8,10,11],[2.350847235165232e7,0.]],
            [[3,7,8,8,8],[1.175423617582616e7,0.]],
            [[3,7,8,9,9],[-1.175423617582616e7,0.]],
            [[3,7,8,10,10],[-1.175423617582616e7,0.]],
            [[3,7,8,11,11],[1.175423617582616e7,0.]],
            [[3,4,4,8,11],[4.701694470330463e7,0.]],
            [[3,4,5,9,11],[-2.350847235165232e7,0.]],
            [[3,4,6,10,11],[-2.350847235165232e7,0.]],
            [[3,4,7,8,8],[-1.175423617582616e7,0.]],
            [[3,4,7,9,9],[1.175423617582616e7,0.]],
            [[3,4,7,10,10],[1.175423617582616e7,0.]],
            [[3,4,7,11,11],[-1.175423617582616e7,0.]],
            [[3,5,5,8,11],[-2.350847235165232e7,0.]],
            [[3,6,6,8,11],[-2.350847235165232e7,0.]],
            [[3,7,7,8,11],[-2.350847235165232e7,0.]],
            [[3,4,4,4,11],[-2.350847235165232e7,0.]],
            [[3,4,4,7,8],[-1.175423617582616e7,0.]],
            [[3,4,5,5,11],[2.350847235165232e7,0.]],
            [[3,4,6,6,11],[2.350847235165232e7,0.]],
            [[3,4,7,7,11],[2.350847235165232e7,0.]],
            [[3,5,5,7,8],[1.175423617582616e7,0.]],
            [[3,6,6,7,8],[1.175423617582616e7,0.]],
            [[3,7,7,7,8],[1.175423617582616e7,0.]],
            [[3,4,4,4,7],[1.175423617582616e7,0.]],
            [[3,4,5,5,7],[-1.175423617582616e7,0.]],
            [[3,4,6,6,7],[-1.175423617582616e7,0.]],
            [[3,4,7,7,7],[-1.175423617582616e7,0.]],
            [[0,0,4,8,8],[-1.175423617582616e7,0.]],
            [[0,0,4,9,9],[-1.175423617582616e7,0.]],
            [[0,0,4,10,10],[-1.175423617582616e7,0.]],
            [[0,0,4,11,11],[-1.175423617582616e7,0.]],
            [[0,0,5,8,9],[2.350847235165232e7,0.]],
            [[0,0,6,8,10],[2.350847235165232e7,0.]],
            [[0,0,7,8,11],[2.350847235165232e7,0.]],
            [[0,0,4,4,8],[2.350847235165232e7,0.]],
            [[0,0,5,5,8],[-2.350847235165232e7,0.]],
            [[0,0,6,6,8],[-2.350847235165232e7,0.]],
            [[0,0,7,7,8],[-2.350847235165232e7,0.]],
            [[0,0,4,4,4],[-1.175423617582616e7,0.]],
            [[0,0,4,5,5],[1.175423617582616e7,0.]],
            [[0,0,4,6,6],[1.175423617582616e7,0.]],
            [[0,0,4,7,7],[1.175423617582616e7,0.]],
            [[0,3,4,8,11],[2.350847235165232e7,0.]],
            [[0,3,5,9,11],[-2.350847235165232e7,0.]],
            [[0,3,6,10,11],[-2.350847235165232e7,0.]],
            [[0,3,7,8,8],[-1.175423617582616e7,0.]],
            [[0,3,7,9,9],[1.175423617582616e7,0.]],
            [[0,3,7,10,10],[1.175423617582616e7,0.]],
            [[0,3,7,11,11],[-1.175423617582616e7,0.]],
            [[0,3,4,4,11],[-2.350847235165232e7,0.]],
            [[0,3,5,5,11],[2.350847235165232e7,0.]],
            [[0,3,6,6,11],[2.350847235165232e7,0.]],
            [[0,3,7,7,11],[2.350847235165232e7,0.]],
            [[0,3,4,4,7],[1.175423617582616e7,0.]],
            [[0,3,5,5,7],[-1.175423617582616e7,0.]],
            [[0,3,6,6,7],[-1.175423617582616e7,0.]],
            [[0,3,7,7,7],[-1.175423617582616e7,0.]]]
        }
        })
    ee_to_dd_2l_bubble.export('ee_to_dd_2l_bubble.yaml')

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

    gamma_to_dd_1l = SquaredTopologyGenerator([('q1', 0, 1), ('p1', 1, 2), ('p2', 1, 2), ('q2', 2, 3)], "gamma_to_dd_1l", ['q1'], 0,
        {'q1': [2., 0., 0., 0.], 'q2': [2., 0., 0., 0.]},
        particle_ids={ 'p1': pdgs['d'], 'p2': pdgs['d~'] },
        overall_numerator=0.25,
        numerator_structure={('p1', 'p2'):
            { ():
                [[[0], (-1.49418e8, 0.)],
                [[0, 0], (7.47091e7, 0.)],
                [[3, 3], (-7.47091e7, 0.)]]
            }
        }
        )
    gamma_to_dd_1l.export('gamma_to_dd_1l.yaml')

    ee_to_dd_1l = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 4, 103), ('q4', 4, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 2, 3), ('p4', 3, 4)], 
        "ee_to_dd_1l", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        particle_ids={ 'p1': pdgs['a'], 'p2': pdgs['d'], 'p3': pdgs['d~'], 'p4': pdgs['a'] },
        overall_numerator=0.25,
        numerator_structure={('p2', 'p3'):
            { ():
                [[[0], (-1.49418e8, 0.)],
                [[0, 0], (7.47091e7, 0.)],
                [[3, 3], (-7.47091e7, 0.)]]
            }
        }
        )
    ee_to_dd_1l.export('ee_to_dd_1l.yaml')

    ee_to_dd_2l_doubletriangle = SquaredTopologyGenerator([('q1', 101, 1), ('q2', 102, 1), ('q3', 6, 103), ('q4', 6, 104), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 5), ('p4', 5, 6),
    ('p5',5, 4), ('p6', 4, 2), ('p7', 4, 3)],
        "ee_to_dd_2l_doubletriangle", ['q1', 'q2'], 2, {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]},
        loop_momenta_names=('p2', 'p3'),
        particle_ids={ 'p1': pdgs['a'], 'p2': pdgs['d'], 'p3': pdgs['d'], 'p4': pdgs['a'], 'p5': pdgs['d~'], 'p6': pdgs['d~'], 'p7': pdgs['g'] },
        overall_numerator=1.0,
#        cut_filter={('p3', 'p5')},
        numerator_structure={('p2', 'p5', 'p7'):
            { (): # uv structure
            [[[0,4],[0.,+2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,5],[0.,-2.954161761482786e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,6],[0.,-2.954161761482786e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,-1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,+1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,+1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,+1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]]
            },
            ('p3', 'p6', 'p7'):
            { (): # uv structure
            [[[0,4],[0.,+2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,5],[0.,-2.954161761482786e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,6],[0.,-2.954161761482786e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,-1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,+1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,+1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,+1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]]
            },
            ('p2', 'p6'):
            {():
            [[[0,4],[0.,-2.954161761482786e8]],
            [[0,4,4],[0.,+1.477080880741393e8]],
            [[0,7,7],[0.,-1.477080880741393e8]],
            [[1,4,5],[0.,-1.477080880741393e8]],
            [[2,4,6],[0.,-1.477080880741393e8]],
            [[3,7],[0.,-2.954161761482786e8]],
            [[0,0,4],[0.,+1.477080880741393e8]],
            [[0,0,4,4],[0.,-7.385404403706966e7]],
            [[0,0,5,5],[0.,+7.385404403706966e7]],
            [[0,0,6,6],[0.,+7.385404403706966e7]],
            [[0,0,7,7],[0.,+7.385404403706966e7]],
            [[0,1,5],[0.,-1.477080880741393e8]],
            [[0,1,4,5],[0.,+1.477080880741393e8]],
            [[0,2,6],[0.,-1.477080880741393e8]],
            [[0,2,4,6],[0.,+1.477080880741393e8]],
            [[1,1,5,5],[0.,-1.477080880741393e8]],
            [[1,2,5,6],[0.,-2.954161761482786e8]],
            [[1,3,5,7],[0.,-1.477080880741393e8]],
            [[2,2,6,6],[0.,-1.477080880741393e8]],
            [[2,3,6,7],[0.,-1.477080880741393e8]],
            [[3,3,4],[0.,-1.477080880741393e8]],
            [[3,3,4,4],[0.,+7.385404403706966e7]],
            [[3,3,5,5],[0.,-7.385404403706966e7]],
            [[3,3,6,6],[0.,-7.385404403706966e7]],
            [[3,3,7,7],[0.,-7.385404403706966e7]]],
            ('p3',):
            [[[0],[0.,+4.726658818372458e9]],
            [[0,4,4],[0.,-1.477080880741393e8]],
            [[0,7,7],[0.,+1.477080880741393e8]],
            [[1,4,5],[0.,+1.477080880741393e8]],
            [[2,4,6],[0.,+1.477080880741393e8]],
            [[0,0],[0.,-2.363329409186229e9]],
            [[0,0,4,4],[0.,+7.385404403706966e7]],
            [[0,0,5,5],[0.,-7.385404403706966e7]],
            [[0,0,6,6],[0.,-7.385404403706966e7]],
            [[0,0,7,7],[0.,-7.385404403706966e7]],
            [[0,1,4,5],[0.,-1.477080880741393e8]],
            [[0,2,4,6],[0.,-1.477080880741393e8]],
            [[1,1,5,5],[0.,+1.477080880741393e8]],
            [[1,2,5,6],[0.,+2.954161761482786e8]],
            [[1,3,5,7],[0.,+1.477080880741393e8]],
            [[2,2,6,6],[0.,+1.477080880741393e8]],
            [[2,3,6,7],[0.,+1.477080880741393e8]],
            [[3,3],[0.,+2.363329409186229e9]],
            [[3,3,4,4],[0.,-7.385404403706966e7]],
            [[3,3,5,5],[0.,+7.385404403706966e7]],
            [[3,3,6,6],[0.,+7.385404403706966e7]],
            [[3,3,7,7],[0.,+7.385404403706966e7]]]
            },
            ('p3', 'p5'):
            {():
                [[[0,4],[0.,-2.954161761482786e8]],
                [[0,4,4],[0.,+1.477080880741393e8]],
                [[0,7,7],[0.,-1.477080880741393e8]],
                [[1,4,5],[0.,-1.477080880741393e8]],
                [[2,4,6],[0.,-1.477080880741393e8]],
                [[3,7],[0.,-2.954161761482786e8]],
                [[0,0,4],[0.,+1.477080880741393e8]],
                [[0,0,4,4],[0.,-7.385404403706966e7]],
                [[0,0,5,5],[0.,+7.385404403706966e7]],
                [[0,0,6,6],[0.,+7.385404403706966e7]],
                [[0,0,7,7],[0.,+7.385404403706966e7]],
                [[0,1,5],[0.,-1.477080880741393e8]],
                [[0,1,4,5],[0.,+1.477080880741393e8]],
                [[0,2,6],[0.,-1.477080880741393e8]],
                [[0,2,4,6],[0.,+1.477080880741393e8]],
                [[1,1,5,5],[0.,-1.477080880741393e8]],
                [[1,2,5,6],[0.,-2.954161761482786e8]],
                [[1,3,5,7],[0.,-1.477080880741393e8]],
                [[2,2,6,6],[0.,-1.477080880741393e8]],
                [[2,3,6,7],[0.,-1.477080880741393e8]],
                [[3,3,4],[0.,-1.477080880741393e8]],
                [[3,3,4,4],[0.,+7.385404403706966e7]],
                [[3,3,5,5],[0.,-7.385404403706966e7]],
                [[3,3,6,6],[0.,-7.385404403706966e7]],
                [[3,3,7,7],[0.,-7.385404403706966e7]]],
            ('p2',):
                [[[0],[0.,4.726658818372458e9]],
                [[0,4,4],[0.,-1.477080880741393e8]],
                [[0,7,7],[0.,1.477080880741393e8]],
                [[1,4,5],[0.,1.477080880741393e8]],
                [[2,4,6],[0.,1.477080880741393e8]],
                [[0,0],[0.,-2.363329409186229e9]],
                [[0,0,4,4],[0.,7.385404403706966e7]],
                [[0,0,5,5],[0.,-7.385404403706966e7]],
                [[0,0,6,6],[0.,-7.385404403706966e7]],
                [[0,0,7,7],[0.,-7.385404403706966e7]],
                [[0,1,4,5],[0.,-1.477080880741393e8]],
                [[0,2,4,6],[0.,-1.477080880741393e8]],
                [[1,1,5,5],[0.,1.477080880741393e8]],
                [[1,2,5,6],[0.,2.954161761482786e8]],
                [[1,3,5,7],[0.,1.477080880741393e8]],
                [[2,2,6,6],[0.,1.477080880741393e8]],
                [[2,3,6,7],[0.,1.477080880741393e8]],
                [[3,3],[0.,2.363329409186229e9]],
                [[3,3,4,4],[0.,-7.385404403706966e7]],
                [[3,3,5,5],[0.,7.385404403706966e7]],
                [[3,3,6,6],[0.,7.385404403706966e7]],
                [[3,3,7,7],[0.,7.385404403706966e7]]]
            },
            }
        )
    ee_to_dd_2l_doubletriangle.export('ee_to_dd_2l_doubletriangle.yaml')

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
