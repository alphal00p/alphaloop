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
        
        mercedes_edges = [
            ('q1',1,2),
            ('q2',1,3),
            ('q3',1,4),
            ('q4',2,3),
            ('q5',3,4),
            ('q6',4,2),
        ]
        all_topos = []
        for lmb in [('q1','q2','q5'),('q1','q2','q4'),('q1','q4','q5'),('q4','q5','q6')]:
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
    def test_combination(args):
        
        n_loops, sig_map, edges, topo_identifier, uvdenoms_placeholder, base_num_rank, min_power_per_edge, test_template, pc = args

        fmbs = [sp.symbols("fmb%d"%(i+1)) for i in range(n_loops)]
        mUV = sp.symbols("mUV2")
        basenum = IntegratedUVTester.to_dot_products(sum(fmbs)**2)**(base_num_rank//2)
        
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

    def test(self, max_total_denom_powers = 8, max_power_per_edge = 2, min_power_per_edge = 1, base_num_rank=8):
        topos = self.get_three_loop_topo_generators()

        test_template = open(pjoin(root_path,'test_template.frm'),'r').read()

        for i_topo, topo in enumerate(topos):
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
                        ) for e, sgn in edges
                    )
                ) for n, edges in nodes.items()
            )
            uvdenoms_placeholder = '*'.join( 'uvprop(%s,%s)'%(
                    ''.join('%sfmb%d'%('+' if wgt>0 else '-', i_mom+1) for i_mom, wgt in enumerate(sig_map[e][0]) if wgt != 0),
                    '%d'
                ) for e in edges )

            # Build a list of edge powers and momenta combinations
            power_combinations = list(
                el for el in itertools.product(*([list(range(max_total_denom_powers+1)),]*len(edges))) if 
                all(e<=max_power_per_edge for e in el) and all(e>=min_power_per_edge for e in el) and sum(el)<=max_total_denom_powers and any(e>min_power_per_edge for e in el)
            )
            logger.info("Will consider the following %d combinations of denominator powers:\n%s"%(len(power_combinations),str(power_combinations)))

            with progressbar.ProgressBar(
                    prefix = 'Integrated UV tests for routing #%d (loop momenta: %s): '%(i_topo+1, ','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta)),
                    max_value=len(power_combinations)+1) as bar:
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
                tests_it = pool.imap(IntegratedUVTester.test_combination, [
                    (topo.n_loops, sig_map, edges, topo_identifier, uvdenoms_placeholder, base_num_rank, min_power_per_edge, test_template, pc) for pc in power_combinations
                ])
                n_tests_done = 0
                for res, target_pc, random_file_name in tests_it:
                    n_tests_done += 1
                    bar.update(n_tests_done)
                    if res != '0':
                        logger.critical("Test failed for routing #%d (loop momenta: %s) when testing (%s) -> (%s). Test file name: %s"%(
                            i_topo, ','.join(topo.edge_map_lin[l][0] for l in topo.loop_momenta),
                            ','.join('%d'%p for p in pc),','.join('%d'%p for p in target_pc),
                            form_test_filename))
                        sys.exit(1)

if __name__ == '__main__':
    IntegratedUVTester().test()
