#!/usr/bin/env python3

import sys
import os
import stat
pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
template_dir = pjoin(root_path,os.path.pardir,os.path.pardir,'alpha_loop','Templates')
sys.path.insert(0, pjoin(root_path,os.path.pardir,os.path.pardir) )
import LTD.ltd_utils
import sys
import sympy as sp
import itertools
import re
import subprocess
import uuid
import argparse
import multiprocessing
import progressbar
import shutil
from pprint import pprint, pformat

import logging
logger = logging.getLogger("IntegratedUVTester")
logger.setLevel(logging.INFO)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)

regexp_times = re.compile(r'(\w+\d+)\*(\w+\d+)')
def repl_times(mach_obj):
    return 'DOT_%s_%s_'%(mach_obj.group(1),mach_obj.group(2))

regexp_power = re.compile(r'(\w+\d+)\*\*(\d+)')
def repl_power(mach_obj):
    if mach_obj.group(2)!='2':
        logger.critical("Powers above two are not safe when converting back to dot products.")
        sys.exit(1)
    return 'DOT_%s_%s_'%(mach_obj.group(1),mach_obj.group(1))

regexp_dot = re.compile(r'DOT\_([\w|\d]+)\_([\w|\d]+)\_')
def repl_dot(mach_obj):
    return 'g(%s,%s)'%(mach_obj.group(1),mach_obj.group(2))
def repl_dot_to_raw(mach_obj):
    return "(%s.%s)"%(mach_obj.group(1).replace('fmb','k'),mach_obj.group(2).replace('fmb','k'))

#a='fmb2*fm4+fm5*fmb2+fmb4**2'
#a=re.sub(regexp_times,repl_times,a)
#a=re.sub(regexp_power,repl_power,a)
#print(a)

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

class IntegratedUVTester(object):
    
    def __init__(self, *args, **opts):
        pass
    
    def get_three_loop_topo_generators(self):
        
        # mercedes_edges = [
        #     ('q1',1,2),
        #     ('q2',1,3),
        #     ('q3',1,4),
        #     ('q4',2,3),
        #     ('q5',3,4),
        #     ('q6',4,2),
        # ]
        # Use the same convention as in rust, for convenience
        mercedes_edges = [
            ('q1',1,2),
            ('q2',2,3),
            ('q3',3,1),
            ('q4',2,4),
            ('q5',3,4),
            ('q6',4,1),
        ]
        all_topos = []
        tmp_topo = LTD.ltd_utils.TopologyGenerator(mercedes_edges)
        lmbs_to_consider = [ tuple(sorted([tmp_topo.edge_map_lin[e_lmb][0] for e_lmb in lmb])) for lmb in tmp_topo.loop_momentum_bases() ]
        # Place the natural basis first
        lmbs_to_consider.pop(lmbs_to_consider.index(('q1','q2','q3')))
        lmbs_to_consider = [('q1','q2','q3'),]+lmbs_to_consider
        #lmbs_to_consider = [('q1','q2','q5'),('q1','q2','q4'),('q1','q4','q5'),('q4','q5','q6')]
        #print(lmbs_to_consider)
        for lmb in lmbs_to_consider:
            topo_generator = LTD.ltd_utils.TopologyGenerator(mercedes_edges)
            topo_generator.generate_momentum_flow( loop_momenta = lmb )
            all_topos.append(topo_generator)

        return all_topos

    @staticmethod
    def to_dot_products(sp_expr):
        sp_expr = str(sp_expr.expand())
        sp_expr = re.sub(regexp_times,repl_times,sp_expr)
        sp_expr = re.sub(regexp_power,repl_power,sp_expr)
        return sp.parse_expr(sp_expr)

    @staticmethod
    def format_numerator(sp_expr):
        sp_expr = sp_expr.expand()
        res = str(sp_expr).replace('**','^')
        res = re.sub(regexp_dot,repl_dot,res)
        return res

    @staticmethod
    def parse_FORM_output(FORM_output):
        
        res = []
        for line in FORM_output.split('\n'):
            if len(res)==0 and not line.strip().startswith("F ="):
                continue
            tmp_line = line.strip()
            if line.strip().startswith("F ="):
                tmp_line = tmp_line[3:]
            if tmp_line.endswith(';'):
                res.append(tmp_line[:-1])
                break
            else:
                if tmp_line.endswith('\\'):
                    res.append(tmp_line[:-1])
                else:
                    res.append(tmp_line)

        return ''.join(l.strip() for l in res)

    @staticmethod
    def read_FORM_output(FORM_str):
        
        if FORM_str.strip()=='0':
            return {}
        try:
            ep = sp.Symbol('ep')
            all_coeffs = sp.Poly((sp.parse_expr(FORM_str.strip().replace('rat','').replace('^','**'))*ep**20).expand(),ep).all_coeffs()
            return {(len(all_coeffs)-(20+1))-i_term: float(val) for i_term, val in enumerate(all_coeffs) if float(val)!=0.}
        except Exception as e:
            logger.critical("Could not parse FORM string output evaluation '%s'. Exception: %s"%(FORM_str, str(e)))
            sys.exit(1)

    @staticmethod
    def test_num_denom_cancellation(args):
        
        n_loops, sig_map, edges, topo_identifier, uvdenoms_placeholder, base_num_rank, min_power_per_edge, test_template, pc = args

        fmbs = [sp.symbols("fmb%d"%(i+1)) for i in range(n_loops)]
        mUV = sp.symbols("mUV2")
        if isinstance(base_num_rank,int):
            basenum = (IntegratedUVTester.to_dot_products(sum(fmbs)**2)+mUV**2)**(base_num_rank//2)
        else:
            basenum = sp.parse_expr(base_num_rank)

        target_pc = [ min(p,min_power_per_edge) for p in pc ]

        #logger.info("Now testing (%s) -> (%s)"%(','.join('%d'%p for p in pc),','.join('%d'%p for p in target_pc)))

        num_A = IntegratedUVTester.format_numerator(basenum)
        num_B = basenum
        for i_edge, (start_p, target_p) in enumerate(zip(pc,target_pc)):
            if start_p == target_p:
                continue
            if start_p<target_p:
                logger.critical("Logical mistake in test.")
                sys.exit(1)
            num_B *= ( IntegratedUVTester.to_dot_products(sum( fmbs[i_mom]*wgt for i_mom, wgt in enumerate(sig_map[edges[i_edge]][0]) if wgt!=0 )**2) - mUV**2)**(start_p-target_p)
        num_B = IntegratedUVTester.format_numerator(num_B)
        F = '%s*(\n\t%s*(%s)\n\t-%s*(%s))'%(
            topo_identifier,
            uvdenoms_placeholder%tuple(target_pc),
            num_A,
            uvdenoms_placeholder%tuple(pc),
            num_B
        )
        #logger.info("About to process:\n%s"%F)
        random_file_name = pjoin(root_path,'form_test_%s.frm'%str(uuid.uuid4()))
        with open(random_file_name,'w') as f:
            f.write(test_template%F)
        FORM_output = subprocess.check_output(['form',random_file_name],cwd=root_path).decode()
        res = IntegratedUVTester.parse_FORM_output(FORM_output)
        if res == '0':
            try:
                os.remove(random_file_name)
            except Exception as e:
                pass
        
        return (res, target_pc, random_file_name)

    def test_num_cancel_denom(self, max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=2, power_combs=None, lmbs=None, n_cores=None):
        topos = self.get_three_loop_topo_generators()
        logger.info("A total of %d topologies/LMB will be considered for this test."%len(topos))

        test_template = open(pjoin(root_path,'test_template.frm'),'r').read()

        for i_topo, topo in enumerate(topos):
            if lmbs is not None and i_topo not in lmbs:
                continue
            logger.info("Now testing routing #%d (loop momenta: %s) ..."%(i_topo+1, ','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta) ))

            sig_map = topo.get_signature_map()

            # Build a list of edges
            edges = sorted(list(sig_map.keys()))
            # Build a list of nodes
            nodes = {}
            for edge_name, start_node, end_node in topo.edge_map_lin:
                if start_node==end_node:
                    logger.critical("ERROR : tadpoles not supported")
                    sys.exit(1)
                if end_node in nodes:
                    nodes[end_node].append((edge_name,+1))
                else:
                    nodes[end_node] = [(edge_name,+1),]
                if start_node in nodes:
                    nodes[start_node].append((edge_name,-1))
                else:
                    nodes[start_node] = [(edge_name,-1),]            
            topo_identifier = '*'.join(
                'vxs(%s)'%(
                    ','.join(
                        ''.join(
                            '%sfmb%d'%('+' if sgn*wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_map[e][0]) if wgt != 0
                        ) for e, sgn in es
                    )
                ) for n, es in nodes.items()
            )
            uvdenoms_placeholder = '*'.join( 'uvprop(%s,%s)'%(
                    ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_map[e][0]) if wgt != 0),
                    '%d'
                ) for e in edges )

            # Build a list of edge powers and momenta combinations
            if power_combs is None:
                power_combinations = list(
                    el for el in itertools.product(*([list(range(max_total_denom_powers+1)),]*len(edges))) if 
                    all(e<=max_power_per_edge for e in el) and all(e>=min_power_per_edge for e in el) and sum(el)<=max_total_denom_powers and any(e>min_power_per_edge for e in el)
                )
            else:
                power_combinations = power_combs

            logger.info("Will consider the following %d combinations of denominator powers:\n%s"%(len(power_combinations),str(power_combinations)))

            with progressbar.ProgressBar(
                    prefix = 'Integrated UV num cancel tests for routing #%d (loop momenta: %s): '%(i_topo+1, ','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta)),
                    max_value=len(power_combinations)) as bar:
                pool = multiprocessing.Pool(processes=n_cores if n_cores is not None else multiprocessing.cpu_count())
                tests_it = pool.imap(IntegratedUVTester.test_num_denom_cancellation, [
                    (topo.n_loops, sig_map, edges, topo_identifier, uvdenoms_placeholder, base_num_rank, min_power_per_edge, test_template, pc) for pc in power_combinations
                ])
                n_tests_done = 0
                for res, target_pc, form_test_filename in tests_it:
                    n_tests_done += 1
                    bar.update(n_tests_done)
                    if res != '0':
                        logger.critical("Test failed for routing #%d (loop momenta: %s) when testing (%s) -> (%s). Test file name: %s"%(
                            i_topo, ','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta),
                            ','.join('%d'%p for p in pc),','.join('%d'%p for p in target_pc),
                            form_test_filename))
                        sys.exit(1)

    @staticmethod
    def test_one_lmb_change(args):
        
        n_loops, sig_map, edges, topo_identifiers, uvdenoms_placeholders, base_num_rank, test_template, test_combination = args
        pc, basis_change, i_lmb = test_combination

        fmbs = [sp.symbols("fmb%d"%(i+1)) for i in range(n_loops)]
        mUV = sp.symbols("mUV2")
        if isinstance(base_num_rank,int):
            basenum = (IntegratedUVTester.to_dot_products(sum(fmbs)**2)+mUV**2)**(base_num_rank//2)
        else:
            basenum = sp.parse_expr(base_num_rank)
        # Example of simple handcrafted numerator
        #basenum = IntegratedUVTester.to_dot_products(fmbs[0]*fmbs[0])

        # Inefficient, but who cares...
        all_dot_products_substitutions = {}
        for i in range(len(basis_change)):
            for j in range(len(basis_change)):
                all_dot_products_substitutions['DOT_fmb%d_fmb%d_'%(i+1,j+1)] = \
                    '(%s)'%str(IntegratedUVTester.to_dot_products(sp.parse_expr("(%s)*(%s)"%(basis_change[i],basis_change[j]))))

        num_B = str(basenum)
        # use these three lines to do the replacement
        rep = dict((re.escape(k), v) for k, v in all_dot_products_substitutions.items()) 
        pattern = re.compile("|".join(rep.keys()))
        num_B = pattern.sub(lambda m: rep[re.escape(m.group(0))], num_B)
        num_B = sp.parse_expr(num_B)

        num_A = IntegratedUVTester.format_numerator(basenum)
        num_B = IntegratedUVTester.format_numerator(num_B)

        # pprint(basis_change)
        # pprint(num_A)
        # pprint(num_B)

        # print(i_lmb)
        # pprint(topo_identifiers[0])
        # pprint(topo_identifiers[i_lmb])
        # pprint(uvdenoms_placeholders[0])
        # pprint(uvdenoms_placeholders[i_lmb])
        F = '(%s*(%s)*(%s)\n\t-%s*(%s)*(%s))'%(
            topo_identifiers[0],
            uvdenoms_placeholders[0]%tuple(pc),
            num_A,
            topo_identifiers[i_lmb],
            uvdenoms_placeholders[i_lmb]%tuple(pc),
            num_B
        )
        #logger.info("About to process:\n%s"%F)
        random_file_name = pjoin(root_path,'form_test_%s.frm'%str(uuid.uuid4()))
        with open(random_file_name,'w') as f:
            f.write(test_template%F)

        # logger.info("About to test LMB for powers (%s) and LMB transition #0->#%d, i.e. (%s)->(%s). Test file name: %s"%(
        #     ','.join('%d'%p for p in pc), i_lmb, ','.join(expr for expr in ['fmb1','fmb2','fmb3']),
        #     ','.join(expr for expr in basis_change),
        #     os.path.basename(random_file_name)))

        FORM_output = subprocess.check_output(['form',random_file_name],cwd=root_path).decode()
        res = IntegratedUVTester.parse_FORM_output(FORM_output)
        # logger.info("Got res=%s"%res)
        # stop

        if res == '0':
            try:
                os.remove(random_file_name)
            except Exception as e:
                pass
        
        return (res, basis_change, pc, i_lmb, random_file_name)

    def test_lmb_change(self, max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=6, power_combs=None, lmbs=None, n_cores = None):

        test_template = open(pjoin(root_path,'test_template_lmb_change.frm'),'r').read()

        topos = self.get_three_loop_topo_generators()
        logger.info("A total of %d topologies/LMB will be considered for this test."%len(topos))

        # Build a list of edges
        edges = sorted([e[0] for e in topos[0].edge_map_lin])

        sig_maps = []
        topo_identifiers = []
        uvdenoms_placeholders = []
        basis_changes = []
        for i_topo, topo in enumerate(topos):
            sig_maps.append(topo.get_signature_map())

            if i_topo == 0:
                basis_changes.append(tuple(['fmb%d'%i for i in range(1,topo.n_loops+1)]))
            else:
                basis_changes.append(tuple(
                    [   
                        ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_maps[i_topo][topos[0].edge_map_lin[topos[0].loop_momenta[i]][0]][0]) if wgt!=0)
                        for i in range(topo.n_loops)
                    ]
                ))
                # print([topos[0].loop_momenta[i] for i in range(topo.n_loops)])
                # print([topos[0].edge_map_lin[topos[0].loop_momenta[i]] for i in range(topo.n_loops)])
                # print([topo.loop_momenta[i] for i in range(topo.n_loops)])
                # print([topo.edge_map_lin[topo.loop_momenta[i]] for i in range(topo.n_loops)])
                # print("--")
                # print(sig_maps[0][ topos[0].edge_map_lin[topos[0].loop_momenta[0]][0] ])
                # print(sig_maps[i_topo][ topos[0].edge_map_lin[topos[0].loop_momenta[0]][0] ])
                # print(
                #     ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in 
                #         enumerate(sig_maps[i_topo][ topos[0].edge_map_lin[topos[0].loop_momenta[0]][0] ][0])  if wgt!=0
                #     )
                # )
                # print(basis_changes[-1])

            # Build a list of nodes
            nodes = {}
            for edge_name, start_node, end_node in topo.edge_map_lin:
                if start_node==end_node:
                    logger.critical("ERROR : tadpoles not supported")
                    sys.exit(1)
                if end_node in nodes:
                    nodes[end_node].append((edge_name,+1))
                else:
                    nodes[end_node] = [(edge_name,+1),]
                if start_node in nodes:
                    nodes[start_node].append((edge_name,-1))
                else:
                    nodes[start_node] = [(edge_name,-1),]            
            topo_identifiers.append('*'.join(
                'vxs(%s)'%(
                    ','.join(
                        ''.join(
                            '%sfmb%d'%('+' if sgn*wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_maps[i_topo][e][0]) if wgt != 0
                        ) for e, sgn in es
                    )
                ) for n, es in nodes.items()
            ))
            uvdenoms_placeholders.append( '*'.join( 'uvprop(%s,%s)'%(
                    ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_maps[i_topo][e][0]) if wgt != 0),
                    '%d'
                ) for e in edges ) )

        # Build a list of edge powers and momenta combinations
        if power_combs is None:
            power_combinations = list(
                el for el in itertools.product(*([list(range(max_total_denom_powers+1)),]*len(edges))) if 
                all(e<=max_power_per_edge for e in el) and all(e>=min_power_per_edge for e in el) and sum(el)<=max_total_denom_powers and any(e>min_power_per_edge for e in el)
            )
        else:
            power_combinations = power_combs

        logger.info("Will consider the following %d combinations of denominator powers:\n%s"%(len(power_combinations),str(power_combinations)))
        tests_combination = [(pc, basis_changes[i_lmb], i_lmb) for pc in power_combinations for i_lmb in ( range(1,len(topos)) if lmbs is None else lmbs ) ]
        with progressbar.ProgressBar(
                prefix = 'Integrated UV basis change tests: ',
                max_value=len(tests_combination)) as bar:
            pool = multiprocessing.Pool(processes=n_cores if n_cores is not None else multiprocessing.cpu_count())
            tests_it = pool.imap(IntegratedUVTester.test_one_lmb_change, [
                (topo.n_loops, sig_maps, edges, topo_identifiers, uvdenoms_placeholders, base_num_rank, test_template, tc) for tc in tests_combination
            ])
            n_tests_done = 0
            for res, basis_change, pc, i_lmb, form_test_filename in tests_it:
                n_tests_done += 1
                bar.update(n_tests_done)
                if res != '0':
                    logger.critical("\nLMB change failed when testing powers (%s) and LMB transition #0->#%d, i.e. (%s)->(%s), i.e. (LMB=%s) vs (LMB=%s).\nTest file name: %s"%(
                        ','.join('%d'%p for p in pc), i_lmb, ','.join(expr for expr in basis_changes[0]),
                        ','.join(expr for expr in basis_change),
                        ','.join(topos[0].edge_map_lin[l][0] for l in topos[0].loop_momenta),
                        ','.join(topos[i_lmb].edge_map_lin[l][0] for l in topos[i_lmb].loop_momenta),
                        os.path.basename(form_test_filename)))
                    sys.exit(1)

    @staticmethod
    def write_pySecDec_run(output_dir, numerator_input_path, topology_generator_graph, input_power_list, graph_name, 
                            graph_id, couplings_prefactor_str, additional_overall_factor, params, mUV, template_run, max_epsilon_order_for_3L=0):
        """ Writes out a standalone script that can run this supergraph in pySecDec."""

        g = topology_generator_graph

        couplings_prefactor_str = re.sub(r'rat\((\s*-?\s*\d+),(\s*-?\s*\d+)\)',r'(\1/float(\2))',couplings_prefactor_str)
        couplings_prefactor_str = couplings_prefactor_str.replace('^','**').replace('i_','(1j)')

        processed_params = {}
        for p, v in params.items():
            if isinstance(v, float): 
                processed_params[p] = v
            elif p=='pi':
                processed_params[p] = float(math.pi)

        output_path = pjoin(output_dir,'run_graph_{}.py'.format(graph_id))
        output_numerator_path = pjoin(output_dir,'numerator_{}.txt'.format(graph_id))

        lorentx_indices_count = {'count': 0}
        regexp_lorentz = re.compile(r'([p|k]\d+\.[p|k]\d+)')
        def repl_lorentz(matchobj):
            lorentx_indices_count['count'] += 1
            dot_args = matchobj.group(0).split('.')
            return '(%s(mu%d)*%s(mu%d))'%(dot_args[0],lorentx_indices_count['count'],dot_args[1],lorentx_indices_count['count'])
        regexp_rat = re.compile(r'rat\(([-|+|\s|\d|\w|\*]+),([-|+|\s|\d|\w|\*]+)\)')
        def repl_rat(mach_obj):
            if mach_obj.group(2)!='1':
                raise FormProcessingError("The numerator to feed pySecDec with has a rational coefficient whose denominator is not 1.")
            return '(%s)'%(mach_obj.group(1).replace(' ',''))
        regexp_power = re.compile(r'([\w|\d]+)\.([\w|\d]+)\^([\d]+)')
        def repl_power(mach_obj):
            return '*'.join(['(%s.%s)'%(mach_obj.group(1),mach_obj.group(2)),]*int(mach_obj.group(3)))
        with open(output_numerator_path,'w') as num_out:
            numerator = []
            numerator_lines = []
            with open(numerator_input_path,'r') as num_in:
                for line in num_in.readlines():
                    numerator_lines.append(line.strip())
            # Sadly FORM breaks down line across rat parenthesis, so we are forced to process the whole numerator string at once.
            for line in [''.join(numerator_lines),]:#num_in.readlines():
                orig_line = line
                line = line.strip().replace('ep','eps').replace(' ','')
                line = re.sub(regexp_power,repl_power,line)
                line = line.replace('^','**')
                line = re.sub(regexp_lorentz,repl_lorentz,line)
                line = re.sub(regexp_rat,repl_rat,line)
                if 'rat' in line:
                    raise FormProcessingError("Bug when doing rat substitution in the following FORM expression: Original line:\n%s\nProcessed line:\n%s"%(orig_line,line))
                numerator.append(line)
            num_out.write('\n'.join(numerator))
        lorentz_indices = ['mu%d'%(i+1) for i in range(lorentx_indices_count['count'])]

        internal_edges = []
        external_edges = []
        power_list = []
        propagators = []
        real_parameters = []
        real_parameters_input = []
        default_externals = []
        sig_map=g.get_signature_map()
        all_internal_nodes = set(sum([[e[1],e[2]] for i_e, e in enumerate(g.edge_map_lin) if i_e not in g.ext],[]))
        for i_e, e in enumerate(g.edge_map_lin):
            if i_e in g.ext:
                external_edges.append([e[0],(e[1] if e[1] in all_internal_nodes else e[2])])
                if not all(s==0 for s in list(sig_map[e[0]][1])[:len(sig_map[e[0]][1])//2]):
                    mass_param = 'mUV2' # model.get_particle(particle_PDGs[e[0]]).get('mass')
                    if mass_param.upper()!='ZERO':
                        #default_externals.append( (model['parameter_dict'][mass_param].real, 0., 0., 0.) )
                        default_externals.append( (mUV, 0., 0., 0.) )
                    else:
                        # If there is a single incoming, do not put it onshell
                        if len(g.ext)==2:
                            default_externals.append( (1., 0., 0., 0.) )
                        elif len(g.ext)==1:
                            default_externals.append( (1., 0., 0., ((-1)**len(default_externals))*1.) )
                        else:
                            pass
            else:
                internal_edges.append([e[0],[e[1],e[2]]])
                power_list.append(g.powers[e[0]])
                loop_momentum = None
                for i_l, lsig in enumerate(list(sig_map[e[0]][0])):
                    if lsig == 0:
                        continue
                    if loop_momentum is None:
                        loop_momentum = 'k%d'%(i_l+1) if lsig > 0 else '-k%d'%(i_l+1)
                    else:
                        loop_momentum += '+k%d'%(i_l+1) if lsig > 0 else '-k%d'%(i_l+1)
    
                for externals_list in [
                    list(sig_map[e[0]][1])[:len(sig_map[e[0]][1])//2],
                    list(sig_map[e[0]][1])[len(sig_map[e[0]][1])//2:]
                ]:
                    for i_p, psig in enumerate(externals_list):
                        if psig == 0:
                            continue
                        # Propagator from a purely external tree
                        if loop_momentum is None:
                            loop_momentum = 'p%d'%(i_p+1) if psig > 0 else '-p%d'%(i_p+1)
                        else:
                            loop_momentum += '+p%d'%(i_p+1) if psig > 0 else '-p%d'%(i_p+1)
                
                mass_str = ''
                mass_param = 'mUV2' #model.get_particle(particle_PDGs[e[0]]).get('mass')
                if mass_param.upper()=='ZERO':
                    mass_str = ''
                else:
                    if isinstance(mUV,float):
                        if mass_param not in real_parameters:
                            real_parameters.append(mass_param)
                            #real_parameters_input.append(model['parameter_dict'][mass_param].real)
                            real_parameters_input.append(mUV)
                        mass_str = '-%s**2'%mass_param
                    else:
                        mass_str = '-%s'%mUV

                propagators.append('(%s)**2%s'%(loop_momentum, mass_str))
    
        loop_momenta_str = None
        external_momenta_str = None
        replacement_rules = []
        if loop_momenta_str is None:
            n_loop_momenta = g.n_loops
            n_externals = len(g.ext)//2
            loop_momenta_str = ['k%d'%(i_l+1) for i_l in range(n_loop_momenta)]
            external_momenta_str = ['p%d'%(i_l+1) for i_l in range(n_externals)]
            for i in range(n_externals):
                for j in range(i, n_externals):
                    replacement_rules.append(('p%d*p%d'%(i+1,j+1),'p%d%d'%(i+1,j+1)))
                    real_parameters.append('p%d%d'%(i+1,j+1))

        repl_dict = {}
        repl_dict['graph_name'] = str(graph_name)
        repl_dict['default_externals'] = str(default_externals)

        sorted_edge_names = sorted([(e[0],i_e) for i_e, e in enumerate(g.edge_map_lin)], key=lambda el: el[0])

        # Drawing replacement
        repl_dict['drawing_input_internal_lines'] = str([internal_edges[i_e] for (_, i_e) in sorted_edge_names])
        if input_power_list is None:
            repl_dict['drawing_input_power_list'] = str([power_list[i_e] for (_, i_e) in sorted_edge_names])
        else:
            repl_dict['drawing_input_power_list'] = str(input_power_list)
        repl_dict['drawing_input_external_lines'] = str(external_edges)

        # Loop package replacement
        # Specify propagators according to their sorted names
        repl_dict['propagators'] = str([propagators[i_e] for (_, i_e) in sorted_edge_names])
        repl_dict['loop_momenta'] = str(loop_momenta_str)
        repl_dict['external_momenta'] = str(external_momenta_str)
        repl_dict['lorentz_indices'] = str(lorentz_indices)
        if input_power_list is None:
            repl_dict['power_list'] = str([power_list[i_e] for (_, i_e) in sorted_edge_names])
        else:
            repl_dict['power_list'] = str(input_power_list)
        repl_dict['numerator_path'] = './%s'%(os.path.basename(output_numerator_path))
        repl_dict['replacement_rules'] = str(replacement_rules)
        repl_dict['real_parameters'] = str(real_parameters)
        repl_dict['loop_additional_prefactor'] = '1'
        repl_dict['contour_deformation'] = 'False'
        repl_dict['max_epsilon_order'] = max(3-g.n_loops,0)+max_epsilon_order_for_3L

        # pySecDec integration replacement
        repl_dict['n_loops'] = str(g.n_loops)
        repl_dict['additional_overall_factor'] = str(additional_overall_factor)
        repl_dict['real_parameters_input'] = str(real_parameters_input)
        repl_dict['complex_parameters_input'] = str([])
        repl_dict['couplings_prefactor'] = str(couplings_prefactor_str)
        repl_dict['couplings_values'] = str(processed_params)
        with open(output_path,'w') as f:
            f.write(template_run.format(**repl_dict))

        # Make the file executable
        os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IEXEC)

    @staticmethod
    def compare_FORM_vs_pySecDec(FORM_result, pySecDec_result, n_sigma_threshold):

        all_poles = sorted(list(set(list(pySecDec_result.keys())+list(FORM_result.keys()))))

        res_msg = [
            '\n%sFORM result:%s\n%s%s%s'%(Colour.GREEN,Colour.END,Colour.BLUE, pformat(FORM_result), Colour.END),
            '\n%sPySecDec result:%s\n%s%s%s\n'%(Colour.GREEN,Colour.END,Colour.BLUE, pformat(pySecDec_result), Colour.END),
            ]
        test_result = True

        for pole in all_poles:
            if pole in FORM_result and pole not in pySecDec_result:
                res_msg.append('%-80s : %s'%("Pole %-2d"%pole,"%sPresent in FORM but not pySecDec!%s"%(Colour.RED,Colour.END)))
                test_result = False
            elif pole in pySecDec_result and pole not in FORM_result:
                discr = abs(pySecDec_result[pole][0])/abs(pySecDec_result[pole][1]) if abs(pySecDec_result[pole][1])!=0. else 100
                if discr>n_sigma_threshold:
                    test_result = False
                res_msg.append('%-80s : %s'%("Pole %-2d : 0 (absent) vs %.16e +/- %.2e"%(pole, abs(pySecDec_result[pole][0]), abs(pySecDec_result[pole][1])),
                    '%s%.3g %s %.1g sigmas%s'%(Colour.RED if discr>n_sigma_threshold else Colour.GREEN,discr,'>' if discr>n_sigma_threshold else '<',n_sigma_threshold,Colour.END)))
            else:
                # Never ask for more than 1.e-12 in relative precision as we then become sensitive to the fact that we're doing I/O on double precision
                discr = abs(FORM_result[pole]-pySecDec_result[pole][0])/max( abs(pySecDec_result[pole][1]), abs(FORM_result[pole])*1.0e-12 ) if abs(pySecDec_result[pole][1])!=0. else 100
                if discr>n_sigma_threshold:
                    test_result = False
                res_msg.append('%-80s : %s'%("Pole %-2d : %.16e vs %.16e +/- %.2e"%(pole, FORM_result[pole], pySecDec_result[pole][0], pySecDec_result[pole][1]),
                    '%s%.3g %s %.1g sigmas%s.'%(Colour.RED if discr>n_sigma_threshold else Colour.GREEN,discr,'>' if discr>n_sigma_threshold else '<',n_sigma_threshold,Colour.END)))
        
        return test_result, '\n'.join(res_msg)

    @staticmethod
    def compare_vs_pySecDec(args):
        
        topo, topo_identifier, uvdenoms_placeholder, base_num_rank, test_template, pySecDecTemplate, target_acc, n_sigma_threshold, verbosity, pc = args

        n_loops = topo.n_loops

        #comp_dir_name = pjoin(root_path,'comparison_vs_pySecDec','comparison_%s'%str(uuid.uuid4()))
        comp_dir_name = pjoin(root_path,'comparison_vs_pySecDec','comparison_%s'%('_'.join('%d'%p for p in pc)))
        if not os.path.isdir(comp_dir_name):
            os.mkdir(comp_dir_name)
        
        fmbs = [sp.symbols("fmb%d"%(i+1)) for i in range(n_loops)]
        mUV = sp.symbols("mUV2")
        if isinstance(base_num_rank,int):
            basenum = (IntegratedUVTester.to_dot_products(sum(fmbs)**2)+mUV**2)**(base_num_rank//2)
        else:
            basenum = sp.parse_expr(base_num_rank)
        numerator = IntegratedUVTester.format_numerator(basenum)

        F = '%s*(%s)*(%s)'%(
            topo_identifier,
            uvdenoms_placeholder%tuple(pc),
            numerator
        )

        form_test_path = pjoin(comp_dir_name,'run.frm')
        with open(form_test_path,'w') as f:
            f.write(test_template%F)

        if verbosity>1:
            logger.info("Running form with 'form %s' in '%s'..."%(form_test_path, comp_dir_name))
        FORM_output = subprocess.check_output(['form',form_test_path],cwd=comp_dir_name).decode()
        res = IntegratedUVTester.parse_FORM_output(FORM_output)
        FORM_result = IntegratedUVTester.read_FORM_output(res)
        if verbosity>1:
            logger.info("%sFORM result for power combination (%s):%s\n%s%s%s"%(Colour.BLUE,','.join('%d'%p for p in pc),Colour.END, Colour.GREEN,pformat(FORM_result),Colour.END))

        # If mUV_value set to float it will be kept parameteric, but if hardcoded to 1 it will directly be written as such
        #mUV_value = 1.0
        mUV_value = '1'

        pySecDec_repl_dict = {'numerator_path' : './FORM_numerator.txt'}
        if isinstance(mUV_value, float):
            raw_numerator = str(basenum.expand())
        else:
            raw_numerator = str(sp.parse_expr(str(basenum.expand()).replace('mUV2',mUV_value)).expand())
        raw_numerator = re.sub(regexp_dot,repl_dot_to_raw,raw_numerator)
        raw_numerator = raw_numerator.replace('**','^')
        with open(pjoin(comp_dir_name,'raw_numerator.txt'),'w') as f:
            f.write(raw_numerator)
        
        IntegratedUVTester.write_pySecDec_run(comp_dir_name, pjoin(comp_dir_name,'raw_numerator.txt'), topo, pc, 'test_integrated_UV_%s'%('_'.join('%d'%p for p in pc)), 
                            graph_id=0, couplings_prefactor_str='1', additional_overall_factor=1.0, params={}, mUV=mUV_value, template_run=pySecDecTemplate)

        if os.path.isfile(pjoin(comp_dir_name,'pySecDec_res.txt')):
            os.remove(pjoin(comp_dir_name,'pySecDec_res.txt'))
        cmd = ['./run_graph_0.py','-i','qmc','-r','%g'%target_acc,'-d','./pySecDec_res.txt']
        if verbosity<=2:
            cmd.append('-q')
        if verbosity>1:
            logger.info("Running pySecDec with '%s' in '%s' ..."%(' '.join(cmd),comp_dir_name))
        if verbosity<=2:
            subprocess.run(cmd,cwd=comp_dir_name,capture_output=True)
        else:
            subprocess.call(cmd,cwd=comp_dir_name)

        if not os.path.isfile(pjoin(comp_dir_name,'pySecDec_res.txt')):
            logger.critical("Running pySecDec with '%s' in '%s' failed to produce an output."%(' '.join(cmd),comp_dir_name))
            sys.exit(1)
        
        pySecDec_result = eval(open(pjoin(comp_dir_name,'pySecDec_res.txt'),'r').read())
        pySecDec_result = {k: (v[0].real,v[1].real) for k, v in pySecDec_result.items()}

        if verbosity>1:
            logger.info("%spySecDec result for power combination (%s):%s\n%s%s%s"%(Colour.BLUE,','.join('%d'%p for p in pc),Colour.END, Colour.GREEN,pformat(pySecDec_result),Colour.END))
        
        test_passed, msg = IntegratedUVTester.compare_FORM_vs_pySecDec(FORM_result, pySecDec_result, n_sigma_threshold)
        with open(pjoin(comp_dir_name,'comparison.txt'),'w') as f:
            f.write(msg)

        if verbosity>1:
            logger.info("%sTest result for pc=(%s):%s %s%s%s\n%s"%(Colour.BLUE, ','.join('%d'%p for p in pc), Colour.END, Colour.GREEN if test_passed else Colour.RED, 'PASSED' if test_passed else 'FAILED', Colour.END, msg))

        return (test_passed, comp_dir_name, msg, pc)

    def test_vs_pySecDec(self, max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=6, 
                            power_combs=None, n_cores = None, target_acc=1.0e-06, n_sigma_threshold=3, verbosity=0):


        test_template = open(pjoin(root_path,'test_template_vs_pySecDec.frm'),'r').read()
        pySecDecTemplate = open(pjoin(template_dir,'run_pySecDec_template.py'),'r').read()

        topo = self.get_three_loop_topo_generators()[0]
        logger.info("The topology with LMB edges (%s) will be considered for comparisons against pySecDec."%(','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta)))

        # Build a list of edges
        edges = sorted([e[0] for e in topo.edge_map_lin])

        sig_map = topo.get_signature_map()
        
        # Build a list of nodes
        nodes = {}
        for edge_name, start_node, end_node in topo.edge_map_lin:
            if start_node==end_node:
                logger.critical("ERROR : tadpoles not supported")
                sys.exit(1)
            if end_node in nodes:
                nodes[end_node].append((edge_name,+1))
            else:
                nodes[end_node] = [(edge_name,+1),]
            if start_node in nodes:
                nodes[start_node].append((edge_name,-1))
            else:
                nodes[start_node] = [(edge_name,-1),]

        topo_identifier = '*'.join(
            'vxs(%s)'%(
                ','.join(
                    ''.join(
                        '%sfmb%d'%('+' if sgn*wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_map[e][0]) if wgt != 0
                    ) for e, sgn in es
                )
            ) for n, es in nodes.items()
        )
        uvdenoms_placeholder = '*'.join( 'uvprop(%s,%s)'%(
                ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_map[e][0]) if wgt != 0),
                '%d'
            ) for e in edges )

        # Build a list of edge powers and momenta combinations
        if power_combs is None:
            power_combinations = list(
                el for el in itertools.product(*([list(range(max_total_denom_powers+1)),]*len(edges))) if 
                all(e<=max_power_per_edge for e in el) and all(e>=min_power_per_edge for e in el) and sum(el)<=max_total_denom_powers and any(e>min_power_per_edge for e in el)
            )
        else:
            power_combinations = power_combs

        logger.info("Will consider the following %d combinations of denominator powers:\n%s"%(len(power_combinations),str(power_combinations)))

        if not os.path.isdir(pjoin(root_path,'comparison_vs_pySecDec')):
            os.mkdir(pjoin(root_path,'comparison_vs_pySecDec'))
        else:
            logger.critical("%sWARNING: The test will reuse content of directory 'comparison_vs_pySecDec', remove it or rename it if you want to start from scratch (for example when changin numerator).%s"%(Colour.RED, Colour.END))

        tests_combination = power_combinations
        with progressbar.ProgressBar(
                prefix = 'Comparison vs pySecDec tests: ',
                max_value=len(tests_combination)) as bar:
            pool = multiprocessing.Pool(processes=n_cores if n_cores is not None else multiprocessing.cpu_count())
            tests_it = pool.imap(IntegratedUVTester.compare_vs_pySecDec, [
                (topo, topo_identifier, uvdenoms_placeholder, base_num_rank, test_template, pySecDecTemplate, target_acc, n_sigma_threshold, verbosity, tc) for tc in tests_combination
            ])
            n_tests_done = 0
            bar.update(n_tests_done)
            for test_passed, comp_dir_name, msg, pc in tests_it:
                n_tests_done += 1
                bar.update(n_tests_done)
                if not test_passed:
                    logger.critical("\n%sComparison with pySecDec in '%s' failed for power combination (%s):%s\n%s"%(
                        Colour.RED, pjoin('.','comparison_vs_pySecDec',os.path.dirname(comp_dir_name)), ','.join('%d'%p for p in pc), Colour.END, msg
                    ))
                    #sys.exit(1)
                elif verbosity>0:
                    logger.info("\n%sComparison with pySecDec in '%s' passed for power combination (%s):%s\n%s"%(
                        Colour.GREEN, pjoin('.','comparison_vs_pySecDec',os.path.dirname(comp_dir_name)), ','.join('%d'%p for p in pc), Colour.END, msg
                    ))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Testing script for integrated UV CTs""")
    requiredNamed = parser.add_argument_group('required named arguments')
    # Required options
    #requiredNamed.add_argument(...)
    # Optional options

    parser.add_argument('--verbosity', '-v', dest='verbosity', type=int, default=0,
                        help='Specify test verbosity (default: %(default)s).')
    parser.add_argument('--threshold', '-t', dest='threshold', type=float, default=3.0,
                        help='Threshold in number of standard deviation to classify as a fail (default: %(default)s).')
    parser.add_argument('--required_accuracy', '-r', dest='required_accuracy', type=float, default=1.0e-6,
                        help='Required accuracy for the comparison (default: %(default)s).')

    args = parser.parse_args()

    tester = IntegratedUVTester()
    #tester.test_num_cancel_denom(max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=2, power_combs=None, lmbs=None, n_cores=None)
    #
    # You can also specify the numerator directly using the base_num_rank option as a string; for instance: base_num_rank="(DOT_fmb3_fmb3_-mUV2**2)**3"
    #
    #tester.test_lmb_change(max_total_denom_powers = 9, max_power_per_edge = 3, min_power_per_edge = 0, base_num_rank=8, power_combs=None, lmbs=None, n_cores = None)
    #tester.test_lmb_change(max_total_denom_powers = 6, max_power_per_edge = 4, min_power_per_edge = 0, base_num_rank=0, power_combs=[(3,3,3,3,3,3),], lmbs=None, n_cores = None)
    #
    # Here's a failure
    #
    #tester.test_lmb_change(power_combs=[(1,1,1,1,1,2),],lmbs=[7,],base_num_rank=6, n_cores=1)
    #
    # Which we can further simplify
    #
    #tester.test_lmb_change(power_combs=[(0,-3,1,1,1,1),],lmbs=[1,],base_num_rank='1', n_cores=1)
    #tester.test_lmb_change(power_combs=[(3,3,3,3,3,3),],lmbs=[1,],base_num_rank='DOT_fmb2_fmb2_**3', n_cores=1)

    # tester.test_vs_pySecDec(max_total_denom_powers = 9, max_power_per_edge = 3, min_power_per_edge = 0, 
    #     base_num_rank=0, power_combs=[(1,1,0,1,2,1),(1,1,1,1,1,1),(1,0,0,1,0,1),], n_cores = 3, 
    #     target_acc=args.required_accuracy, n_sigma_threshold=args.threshold, verbosity=args.verbosity)
    
    tester.test_vs_pySecDec(max_total_denom_powers = 7, max_power_per_edge = 2, min_power_per_edge = 1, 
        base_num_rank=2, power_combs=None, n_cores = 8, 
        target_acc=args.required_accuracy, n_sigma_threshold=args.threshold, verbosity=args.verbosity)
