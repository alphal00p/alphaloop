#!/usr/bin/env python3

import sys
import os
pjoin = os.path.join
root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.path.pardir,os.path.pardir) )
import LTD.ltd_utils
import sys
import sympy as sp
import itertools
import re
import subprocess
import uuid
import multiprocessing
import progressbar
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

#a='fmb2*fm4+fm5*fmb2+fmb4**2'
#a=re.sub(regexp_times,repl_times,a)
#a=re.sub(regexp_power,repl_power,a)
#print(a)

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
        
        for line in FORM_output.split('\n'):
            if line.strip().startswith("F = "):
                return line.strip()[4:-1]
        return None

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

if __name__ == '__main__':
    tester = IntegratedUVTester()
    #tester.test_num_cancel_denom(max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=2, power_combs=None, lmbs=None, n_cores=None)
    #
    # You can also specify the numerator directly using the base_num_rank option as a string; for instance: base_num_rank="(DOT_fmb3_fmb3_-mUV2**2)**3"
    #
    tester.test_lmb_change(max_total_denom_powers = 9, max_power_per_edge = 3, min_power_per_edge = 0, base_num_rank=8, power_combs=None, lmbs=None, n_cores = None)
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