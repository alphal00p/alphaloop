#!/usr/bin/env python3

import copy
import logging
import os
from pathlib import Path
from pprint import pprint, pformat
import math
import igraph
import time
import numpy as np
from itertools import combinations_with_replacement
from collections import OrderedDict
import glob as glob_module

import progressbar
from itertools import chain, product
import sys
import functools
import subprocess
import argparse
import shutil
import py_compile
from warnings import catch_warnings
import os
import stat

pjoin = os.path.join

if __name__ == "__main__":
    root_path = os.path.dirname(os.path.realpath( __file__ ))
    sys.path.insert(0, pjoin(root_path,os.path.pardir))
    sys.path.insert(0, pjoin(root_path,os.path.pardir,os.path.pardir,os.path.pardir))
    if 'MG5DIR' in os.environ:
        sys.path.insert(0, os.environ['MG5DIR'])
    else:
        print("\033[91mYou are using ./FORM_processing.py in standalone, it is recommended then "+
              "that you define the environment variable 'MG5DIR' pointing to the root directory of MG5aMC in your system.\033[0m")
        sys.exit(1)

import alpha_loop.utils as utils
import re

import madgraph.core.base_objects as base_objects
import madgraph.various.misc as misc
import madgraph.iolibs.file_writers as writers
from madgraph import MadGraph5Error, InvalidCmd, MG5DIR
import models.model_reader as model_reader
import multiprocessing
import platform

import LTD.squared_topologies
import LTD.ltd_utils
import LTD.partial_fractioning

import alpha_loop.formset as formset

logger = logging.getLogger('alphaLoop.FORM_processing')

try:
    import yaml
    from yaml import Loader, Dumper
except ImportError:
    raise BaseException("Install yaml python module in order to import/export topologies from/to yaml.")

if __name__ == "__main__":
    logging.basicConfig()
    logger.setLevel(logging.INFO)

plugin_path = os.path.dirname(os.path.realpath( __file__ ))

FORM_processing_options = {
    'FORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'form').resolve()),
    'tFORM_path': str(Path(plugin_path).parent.joinpath('libraries', 'form', 'sources', 'tform').resolve()),
    # Define the extra aguments for the compilation
    'compilation-options': [],
    'FORM_parallel_cores': 1,
    'cores': 2, #multiprocessing.cpu_count(),
    'extra-options': {'OPTIMITERATIONS': 1000, 'NUMERATOR': 0, 'SUMDIAGRAMSETS': 'nosum', 'MAXVARSFOROPTIM': 3500},
    # If None, only consider the LMB originally chosen.
    # If positive and equal to N, consider the first N LMB from the list of LMB automatically generated
    # If negative consider all possible LMBs.
    'number_of_lmbs' : None,
    # If None, the reference LMB will be the one originally chosen.
    # If positive and equal to N, the Nth LMB will be used for the reference implementation of the supergraph.
    'reference_lmb' : None,
    'FORM_call_sig_id_offset_for_additional_lmb' : 1000000,
    'generate_arb_prec_output' : False,
    'generate_integrated_UV_CTs' : True,
    'on_shell_renormalisation' : True,
    'perform_msbar_subtraction' : True,
    'uv_test': None,
    'generate_renormalisation_graphs' : False,
    # Set the option below to a positive integer so as to enable the splitting of large source files into many
    # with at most around the number of lines specified here. Note that this disables static and inline optimisations!
    'max_n_lines_in_C_source' : 20000,
    'max_n_lines_in_C_function' : 1000,
    'max_n_files_in_library' : 100,
    'include_integration_channel_info' : True,
    'UV_min_dod_to_subtract' : 0,
    'selected_epsilon_UV_order' : 0,
    # Select how to include the finite part of the renormalisation. Possible values are:
    # a) 'together' : all contributions are ketp.
    # b) 'only' : only the contribution from the finite part of the renormalisation is kept.
    # c) 'removed' : the finitie part of the renormalisation is removed. 
    'renormalisation_finite_terms' : 'together',
    'optimisation_strategy' : 'CSEgreedy',
    'optimise_generation_lmb' : True,
    'optimise_integration_channels' : True,
    'FORM_use_max_mem_fraction': 0.75,
    'FORM_setup': {
    #   'MaxTermSize':'100K',
    #   'Workspace':'1G'
    },
    'optimize_c_output_mode_per_pair_of_factors' : 'flat'
}

dummy_scalar_PDGs = {
    3370 : 'ZERO',
    3371 : 'massdummya',
    3372 : 'massdummyb',
    3373 : 'massdummyc',
    3374 : 'massdummyd',
    3375 : 'massdummye',
    3376 : 'massdummyf',
    3377 : 'massdummyg',
    3378 : 'massdummyh',
    3379 : 'massdummyi'
}

forced_edge_scalings = {}
forced_node_scalings = {}

# Can switch to tmpdir() if necessary at some point
FORM_workspace = pjoin(plugin_path,'FORM_workspace')
Path(FORM_workspace).mkdir(parents=True, exist_ok=True)
resources_to_link = ['diacolor.h']
for resource in resources_to_link:
    if not os.path.exists(pjoin(FORM_workspace,resource)):
        utils.ln(pjoin(plugin_path,resource),starting_dir=FORM_workspace)

# Check if there is a multiline return (i.e. by checking for balanced parenthesis), in which case we reformat it to make it 
# multiple lines so as to help compiler.
def balanced_parenthesis(s):
    """ Returns None if balanced, and otherwise index of first offending parenthesis """
    #pairs = {"{": "}", "(": ")", "[": "]"}
    pairs = {"(": ")"}
    stack = []
    for i_c, c in enumerate(s):
        #if c in "{[(":
        if c=="(":
            stack.append((i_c,c))
        elif stack and c == pairs[stack[-1][1]]:
            stack.pop()
        elif c==")":
            return i_c
    
    return (stack[0][0] if len(stack)>0 else None)

# This is only necessary when running with older versions of FORM that could have terms on multiple lines.
def temporary_preprocess_multiline_blocks(FORM_output):

    def combine_lines(lines, split_terms=False):

        if not split_terms:
            combined_line = ''
            for i_line, line in enumerate(lines):
                l = line if not line.endswith('\\') else line[:-1]
                if i_line == 0:
                    combined_line += l
                else:
                    combined_line += l.strip()

            return [combined_line,]
        else:
            combined_lines = [lines[0],]

            current_running_line = ''
            for i_line, line in enumerate(lines[1:]):
                if not (line.startswith('       - ') or line.startswith('       + ')):
                    if current_running_line == '':
                        combined_lines[-1] = combine_lines([combined_lines[-1],line],split_terms=False)[0]
                        continue
                    else:
                        current_running_line = combine_lines([current_running_line,line],split_terms=False)[0]
                else:
                    current_running_line = combine_lines([current_running_line,line],split_terms=False)[0]

                if len(combined_lines)==1:
                   if balanced_parenthesis(combined_lines[0]) is None:
                       raise FormProcessingError("No imbalanced parenthesis in first line to combine within an output C context.")

                b_par = balanced_parenthesis(current_running_line)
                if b_par is None:
                    combined_lines.append('      _ += %s'%current_running_line)
                    current_running_line = ''
                else:
                    if current_running_line[b_par]==')':
                        if i_line != len(lines)-2:
                            raise FormProcessingError('The closing parenthesis of the an output C block did not happen on the last line of the block.')
                        combined_lines.append('      _ += %s'%current_running_line)
                        current_running_line = ''
                    else:
                        continue

            if current_running_line!='':
                raise FormProcessingError('Accumulated line of a C output block was not added:\n%s'%current_running_line)

            return combined_lines

    new_lines = []
    
    prev_line = None
    lines_in_current_block = []
    is_in_multiline_block = False
    was_last_line_a_return = False
    for i_line, line in enumerate(FORM_output.split('\n')):
        stripped_line = line.strip()
        if '//CMODE' in stripped_line:
            is_in_multiline_block =True
            continue
        # Old detection
        #if stripped_line.startswith('return') and not stripped_line[-1]==';':
        #    lines_in_current_block.append(line)
        #    is_in_multiline_block =True
        #    continue
        if stripped_line==';':
            if not is_in_multiline_block:
                if len(lines_in_current_block)!=0:
                    raise FormProcessingError('ERROR A in preprocessing of multiline blocks of FORM output at line %d, line:\n%s'%(i_line, line))
                new_lines[-1] += ';'
                continue
            if len(lines_in_current_block)==0:
                raise FormProcessingError('ERROR B in preprocessing of multiline blocks of FORM output at line %d, line:\n%s'%(i_line, line))
            new_lines.extend(combine_lines(lines_in_current_block, split_terms=False))
            lines_in_current_block = []
            continue

        if stripped_line.startswith('_ +='):
            if not is_in_multiline_block:
                raise FormProcessingError('ERROR C in preprocessing of multiline blocks of FORM output at line %d, line:\n%s'%(i_line, line))
            if len(lines_in_current_block)!=0:
                raise FormProcessingError('ERROR D in preprocessing of multiline blocks of FORM output at line %d, line:\n%s'%(i_line, line))

        if not is_in_multiline_block:
            new_lines.append(line)
            continue

        lines_in_current_block.append(line)

        if stripped_line[-1]==';':
            new_lines.extend(combine_lines(lines_in_current_block, split_terms=True))
            lines_in_current_block = [] 
            is_in_multiline_block = False

    if is_in_multiline_block or len(lines_in_current_block)>0:
        raise FormProcessingError('ERROR E in preprocessing of multiline blocks of FORM output at line %d, line:\n%s'%(i_line, line))

    return new_lines

def temporary_fix_FORM_output(FORM_output):

    # new_output = []
    # previous_line = None
    # for line in FORM_output.split('\n'):
    #     if line.startswith('      _ +=  '):
    #         line = '      %s'%line[12:]
    #         if previous_line is not None:
    #             if previous_line[:-1].strip()!='':
    #                 new_output.append(previous_line[:-1])
    #         previous_line = line
    #     else:
    #         if previous_line is not None:
    #             if previous_line.strip()!='':
    #                 new_output.append(previous_line)
    #         previous_line = line
    # if previous_line is not None:
    #     new_output.append(previous_line)

    # return '\n'.join(new_output)

    # FORM can still break lines into multiple ones that need to be merged,
    # however within a block of "C output" style, each term within the big overall parenthesis should not be split.
    new_output = temporary_preprocess_multiline_blocks(FORM_output)
    # With the new FORM version we can directly use the current output.
    #new_output = FORM_output.split('\n')
    
    currently_in_unbalanced_context = None
    forest_re = re.compile(r'forestid\((\d+)\)')
    power_re = re.compile(r'pow\((\w*)\,(\d+)\)')
#    power_re = re.compile(r'pow\((mUV)\,(\d+)\)')

    Z_counter = {'count' : 2}
    structures_seen = {}
    known_structures = {
        'logmUVmu' : None,
    }
    def reset_function_context_vars():
        Z_counter['count'] = 2
        structures_seen.clear()
        for k in known_structures:
            known_structures[k] = None

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield tuple(lst[i:i + n])

    def process_line(l, first=False):
        new_lines = []
        #processed_line = l.replace(' ','')
        processed_line = l

        # Optimize forests
        for forestID in list(re.findall(forest_re,l)):
            forest_structure = 'forestid(%s)'%forestID
            if forest_structure not in structures_seen:
                Z_counter['count'] += 1
                structures_seen[forest_structure] = 'Z%d_'%Z_counter['count']
                new_lines.append('%s = %s;'%(structures_seen[forest_structure], forest_structure))
            processed_line = processed_line.replace(forest_structure,structures_seen[forest_structure])

        # Optimize know structures
        for structure, var in list(known_structures.items()):
            if var is None and structure in processed_line:
                Z_counter['count'] += 1
                known_structures[structure] = 'Z%d_'%Z_counter['count']
                new_lines.append('%s = %s;'%(known_structures[structure], structure))
            if structure in processed_line:
                processed_line = processed_line.replace(structure,known_structures[structure])

        # Optimize powers
        for power_arg, power_int in list(re.findall(power_re,l)):
            power_structure = 'pow(%s,%s)'%(power_arg, power_int)
            if power_structure not in structures_seen:
                Z_counter['count'] += 1
                structures_seen[power_structure] = 'Z%d_'%Z_counter['count']
                new_lines.append('%s = %s;'%(structures_seen[power_structure], power_structure))
            processed_line = processed_line.replace(power_structure,structures_seen[power_structure])

        # Iteratively optimize doublets
        if FORM_processing_options['optimize_c_output_mode_per_pair_of_factors'] is not None:
            processed_line=processed_line.replace(' ','')
            sign_in_front = ''
            if processed_line[0] in ['+','-']:
                sign_in_front = processed_line[0]
                processed_line = processed_line[1:]
            factors = processed_line.split('*')
            while len(factors)>1:
                factors.sort()
                new_factors = []
                for factors_pair in chunks(factors,2):
                    if len(factors_pair)==2:
                        if factors_pair not in structures_seen:
                            Z_counter['count'] += 1
                            structures_seen[factors_pair] = 'Z%d_'%Z_counter['count']
                            new_lines.append('%s = %s;'%(structures_seen[factors_pair], '*'.join(factors_pair)))
                        new_factors.append(structures_seen[factors_pair])
                    else:
                        new_factors.append(factors_pair[0])
                factors = new_factors
                if FORM_processing_options['optimize_c_output_mode_per_pair_of_factors']=='flat':
                    break
            processed_line = '%s%s'%(sign_in_front, '*'.join(factors))

        if first:
            new_lines.append('Z2_ = %s'%processed_line+';')
        else:
            new_lines.append('Z2_ += %s'%processed_line+';')
        return new_lines

    processed_output = []
    for i_line, line in enumerate(new_output):
        if currently_in_unbalanced_context is None:
            unbalanced_index = balanced_parenthesis(line)
            if unbalanced_index is None:
                processed_output.append(line)
            else:
                if not line.strip().startswith('return'):
                    raise FormProcessingError("Unbalanced parenthesis for a line that is not a return line, this is not expected. Line #%d: '%s'"%(i_line+1, line))
                currently_in_unbalanced_context = unbalanced_index
                # Reset all global replacement rules
                reset_function_context_vars()

                # Make sure the symbol before the unbalenced parenthesis is indeed '*'
                if line[unbalanced_index]!='(' or line[unbalanced_index-1]!='*':
                    raise FormProcessingError("Unexpected strucuture for line with unbalanced parenthesis. Line: '%s'"%line)
                processed_output.append(line[:(unbalanced_index-1)].replace('return','Z1_ = ')+';')
                processed_output.extend(process_line(line[(unbalanced_index+1):], first=True))
        else:
            if not line.startswith('      _ +='):
            #if not (line.startswith('       +') or line.startswith('       -')):
                raise FormProcessingError("Unexpected format of line in unbalanced context. Line: '%s'"%line)
            unbalanced_index = balanced_parenthesis(line)
            if unbalanced_index is None:
                processed_output.extend(process_line(line[10:]))
                #processed_output.extend(process_line(line[7:]))
            else:
                if line[-1]!=';' or unbalanced_index!=(len(line)-2):
                    raise FormProcessingError("Unexpected end of balanced context. Line: '%s'"%line)
                currently_in_unbalanced_context = None
                processed_output.extend(process_line(line[10:-2]))
                #processed_output.extend(process_line(line[7:-2]))
                processed_output.append('return Z1_*Z2_;')

    return '\n'.join(processed_output)

class FormProcessingError(MadGraph5Error):
    """ Error for the FORM processing phase."""
    pass

class FORMSuperGraph(object):
    """ Simplified SuperGraph object with only the necessary information to have it processed by FORM."""

    _FORM_Feynman_rules_conventions = {
        # SSS
        ( 1, 1, 1 ): (0, 1, 2),
        # SSSS
        ( 1, 1, 1 ): (0, 1, 2),
        # GGG
        ( 3, 3, 3 ): (0, 1, 2),
        # GGGG
        ( 3, 3, 3, 3 ): (0, 1, 2, 3),
        # FxFV
        ( -2, 2, 3 ): (0, 2, 1),
        # GHxGHV
        ( -1, 1, 3 ): (0, 2, 1),
        # FxFS
        ( -2, 1, 2 ): (0, 1, 2),
        # SSS
        ( 1, 1, 1 ): (0, 1, 2),
        # SSSS
        ( 1, 1, 1, 1 ): (0, 1, 2),
        # VxVV (e.g. W+ W- a )
        ( -3, 3, 3 ): (1, 0, 2),
        # Fx F
        (-2, 2): (0, 1),
        # V V
        (-3,3): (0,1),
        (3,3): (0,1),
        # S S
        (1,1): (0,1),
        # S V V
        (1,3,3): (0,1,2),
        # S V V V
        (1,3,3,3): (0,1,2,3),
        # S V V V V
        (1,3,3,3,3): (0,1,2,3,4),
        # S S V V
        (1,1,3,3): (0,1,2,3),
        # S S V V V
        (1,1,3,3,3): (0,1,2,3,4),
        # S S S V V
        (1,1,1,3,3): (0,1,2,3,4),
        # S S S V V V
        (1,1,1,3,3,3): (0,1,2,3,4,5)
    }

    _include_momentum_routing_in_rendering=False
    _include_edge_name_in_rendering=True
    _rendering_size = (1.0*(11.0*60),1.0*(8.5*60)) # 1.0 prefactor should be about 1 landscape A4 format per graph
    # Choose graph layout strategy. Interesting options are in comment.
    _graph_layout_strategy = '{"PackingLayout"->"ClosestPacking"}' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringElectricalEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> "SpringEmbedding"' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredEmbedding", "Orientation" -> Left, "RootVertex" -> "I1"}' 
    #_graph_layout_strategy = 'GraphLayout -> {"LayeredDigraphEmbedding", "Orientation" -> Left}'
    
    # '{"SpringEmbedding"}' gives interesting results too.

    def __init__(self, *args,
        call_identifier=None,
        name=None,
        edges=None,
        nodes=None,
        overall_factor="1",
        multiplicity=1,
        benchmark_result=0.0,
        default_kinematics=None,
        effective_vertex_id=None
    ):
        """ initialize a FORM SuperGraph from several options."""

        self.is_zero = False
        self.edges = edges
        self.nodes = nodes
        self.overall_factor = overall_factor
        self.multiplicity = multiplicity
        # Give the possibility of specifying a benchmark result to the output yaml file for future ref.
        self.benchmark_result = benchmark_result
        # A hashable call signature
        self.call_identifier = call_identifier
        if name is None:
            self.name = str(self.call_identifier)
        else:
            self.name = name

        self.default_kinematics = default_kinematics
        self.effective_vertex_id = effective_vertex_id
        self.squared_topology = None
        self.configurations = None

        # Will be filled in during FORM generation
        self.code_generation_statistics = None

        # Store copies of self for different choices of LMB to be used for cross-check.
        self.additional_lmbs = []

    def filter_valid_cuts(self, cuts):
        """ 
        filter graph base on a list of allowed cuts whose entires are defined as:
            ([allowd pdg], n) : to ensure to cut "n" times particles contained 
                                in the pdg list
            ('any', n)        : cut "n" edges that could be anything 
                                (useful for extra real radiations)
        example: e+e- > ggh @NLO => [([21],2), ([25],1), ('any',1)]
        """
        g = igraph.Graph()
        g.add_vertices(len(self.nodes))
        undirected_edges = set()

        for e in self.edges.values():
            undirected_edges.add(tuple(sorted(e['vertices'])))

        cut_edges = [[] for _ in range(len(cuts))]
        incoming_vertices = []
        outgoing_vertices = []
        edge_colors = {tuple(sorted([x - 1 for x in ue])) : [] for ue in undirected_edges}
        for ue in undirected_edges:
            multiple = 0
            for e in self.edges.values():
                if tuple(sorted(e['vertices'])) == ue:
                    multiple += 1
                    if e['type'] == 'virtual':
                        for ci, c in enumerate(cuts):
                            if c[0] == 'any':
                                cut_edges[ci] += [tuple(sorted([x - 1 for x in ue]))]
                            else:
                                if abs(e['PDG']) in c[0]:
                                    cut_edges[ci] += [tuple(sorted([x - 1 for x in ue]))]
                        edge_colors[tuple(sorted([x - 1 for x in ue]))] += [abs(e['PDG'])]
                    elif e['type'] == 'in':
                        incoming_vertices += [min(e['vertices'])-1]
                    elif e['type'] == 'out':
                        outgoing_vertices += [min(e['vertices'])-1]
            g.add_edges([tuple(sorted([x - 1 for x in ue]))]*multiple)

        take_cuts = []
        cut_colors = []
        n_optional = 0
        valid_cut = False
        for ci, cut in enumerate(cuts):
            cut_colors += [cut[0]] * cut[1]
            if cut[0] == 'any' or None in cut[0]:
                n_optional += cut[1]
                take_cuts += [cut_edges[ci]+[None]] * cut[1]
            else:
                take_cuts += [cut_edges[ci]] * cut[1]
        #print(take_cuts)
        #print(cut_colors)
        #print(incoming_vertices)
        #print(outgoing_vertices)
        invalid_cuts = []
        count_checks = 0
        #print(cuts)
        for cut_edges in product(*take_cuts):
            # When a valid is cut we know we need to keep this graph
            if valid_cut:
                break
            
            # Check if cut has to be dropped based on previous failed attempts
            if any(all(cut_edges.count(c) >= veto_c.count(c) for c in set(veto_c)) for veto_c in invalid_cuts):
                continue
            count_checks += 1
            
            # Check valid color cuts
            ec = copy.deepcopy(edge_colors)
            for ce, cc in zip(cut_edges,cut_colors):
                if cc == 'any' or ce == None:
                    continue
                for color in cc:
                    if color in ec[ce]:
                        ec[ce].remove(color)
                        break
                else:
                    break
            else:
                # Allow for Pure Virtual corrections
                virtual_loops = cut_edges.count(None)
                for _ in range(virtual_loops):
                    cut_edges = list(cut_edges)
                    cut_edges.remove(None)
                    cut_edges = tuple(cut_edges)

                # Apply set of cuts
                gtmp = g.copy()
                for ci in range(len(cut_edges)):
                    if not gtmp.are_connected(*cut_edges[ci]):
                        # update invalid cuts
                        new_veto = cut_edges[:ci+1]
                        invalid_cuts = list(filter(lambda veto_c: not all(veto_c.count(c) >= new_veto.count(c) for c in set(new_veto)),invalid_cuts))
                        invalid_cuts += [new_veto]
                        break
                    gtmp.delete_edges(gtmp.get_eid(*cut_edges[ci]))
                    if not gtmp.is_connected():
                        if ci+1 < len(cut_edges):
                            if ci+1 < len(cut_edges) - n_optional + virtual_loops:
                                # update invalid cuts
                                new_veto = cut_edges[:ci+1]
                                invalid_cuts = list(filter(lambda veto_c: not all(veto_c.count(c) >= new_veto.count(c) for c in set(new_veto)),invalid_cuts))
                                invalid_cuts += [new_veto]
                            break
                        else:
                            # check that the vertices are correctly connected 
                            # to the left and right paths
                            # Meaning two vertices involved in a cut should belong do
                            # the opposite disconnected graphs
                            with catch_warnings(record=True) as caught_warnings:
                                if any(len(gtmp.get_shortest_paths(c[0],to=c[1])[0])>0 for c in cut_edges):
                                    break
                                # check that incoming and outgoing are on separate disconnected graphs
                                if any(len(gtmp.get_shortest_paths(incoming_vertices[0],to=c)[0])==0 for c in incoming_vertices[1:]) or\
                                    any(len(gtmp.get_shortest_paths(incoming_vertices[0],to=c)[0])>0 for c in outgoing_vertices):
                                    break
                                valid_cut = True
        return valid_cut
    
    def filter_valid_cuts_helper(args):
        (iso_graph, cuts) = args
        return (iso_graph, iso_graph[1][0].filter_valid_cuts(cuts))


    def get_mathematica_rendering_code(self, model, FORM_id=None, lmb_id=None):
        """ Generate mathematica expression for drawing this graph."""

        repl_dict = {
            'edge_font_size'     : 10,
            'vetex_font_size'    : 10,
            'width'              : self._rendering_size[0],
            'height'             : self._rendering_size[1],
            'graph_layout_strategy' : self._graph_layout_strategy
        }
        if self.call_identifier and all(k in self.call_identifier for k in ['proc_id','left_diagram_id','right_diagram_id']):
            graph_name='MG: %s'%('P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%self.call_identifier)
        else:
            graph_name='MG: %s'%self.name
        if FORM_id is not None:
            graph_name +=' | FORM: #%d'%FORM_id
        if lmb_id is not None:
            graph_name +=' | LMB: #%d'%lmb_id
        repl_dict['graph_name'] = graph_name

        # Special name rendering rules
        def get_part_name(pdg):
            quark_names = {1:"d",2:"u",3:"s",4:"c",5:"b",6:"t"}
            if pdg==11:
                return r"\!\(\*SuperscriptBox[\(e\), \(+\)]\)"
            elif pdg==-11:
                return r"\!\(\*SuperscriptBox[\(e\), \(-\)]\)"
            elif pdg==22:
                return r"\[Gamma]"
            elif pdg in [-1,-2,-3,-4,-5,-6]:
                return r"\!\(\*OverscriptBox[\(%s\), \(_\)]\)"%quark_names[abs(pdg)]
            else:
                return model.get_particle(pdg).get_name()

        node_key_to_node_name = {}
        def get_node_label(node_key):
            if node_key in node_key_to_node_name:
                return node_key_to_node_name[node_key]
            node = self.nodes[node_key]
            if 'renormalisation_vertex_n_loops' in node:
                label = 'UV%dL'%node['renormalisation_vertex_n_loops']
                label += '%dP'%node['renormalisation_vertex_n_shrunk_edges']
                # If that label is not unique, then prefix it the actual key.
                if label in node_key_to_node_name.values():
                    label = '%s_%s'%(str(node_key),label)
            else:
                label = str(node_key)

            node_key_to_node_name[node_key] = label
            return label

        # Generate edge list
        edge_template = """Labeled[Style[CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]%(edge_style)s,%(edge_color)s,Thickness[%(thickness)f]],"%(edge_label)s"]"""
        all_edge_definitions = []
        all_edge_shape_functions = []
        for edge_key, edge_data in self.edges.items():
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key)
            edge_repl_dict = {}
            # is_LMB = ('name' in edge_data and 'LMB' in str(edge_data['name']).upper())
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            is_LMB = (sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0)
            if is_LMB:
                edge_repl_dict['thickness'] = 0.005
            else:
                edge_repl_dict['thickness'] = 0.002
            edge_repl_dict['in_node'] = get_node_label(edge_key[0])
            edge_repl_dict['out_node'] = get_node_label(edge_key[1])
            edge_repl_dict['edge_key'] = edge_key[2]
            edge_repl_dict['arrow_style'] = 'Arrow' if not is_LMB else 'HalfFilledDoubleArrow'
            edge_repl_dict['arrow_size'] = 0.015 if not is_LMB else 0.025
            all_edge_shape_functions.append(
                'CreateEdge["%(in_node)s","%(out_node)s",%(edge_key)d]->GraphElementData["%(arrow_style)s", "ArrowSize" -> %(arrow_size)f]'%edge_repl_dict
            )
            if 'name' in edge_data and 'CUT' in str(edge_data['name']).upper():
                edge_repl_dict['edge_style'] = ",Dashed"
            else:
                edge_repl_dict['edge_style'] = ""
            if edge_data['PDG'] in [-1,-2,-3,-4,-5,1,2,3,4,5]:
                color = "Cyan"
            elif edge_data['PDG'] in [-6,6]:
                color = "Blue"
            elif edge_data['PDG'] in [21,]:
                color = "Red"
            elif edge_data['PDG'] in [82,-82]:
                color = "Pink"
            elif edge_data['PDG'] in [25,]:
                color = "Green"
            else:
                color = "Gray"
            edge_repl_dict['edge_color'] = color
            edge_label_pieces = ['psi' if edge_data['PDG'] in dummy_scalar_PDGs else get_part_name(edge_data['PDG']),]
            if 'name' in edge_data and self._include_edge_name_in_rendering:
                edge_label_pieces.append(edge_data['name'])
            if self._include_momentum_routing_in_rendering:
                edge_label_pieces.append(edge_data['momentum'])
            if is_LMB:
                edge_label_pieces.append('#%d'%(abs_sig[0].index(1)))
            edge_label = "|".join(edge_label_pieces)
            edge_repl_dict['edge_label'] = edge_label
            all_edge_definitions.append(edge_template%edge_repl_dict)

        repl_dict['edge_lists'] = ',\n'.join(all_edge_definitions)
        repl_dict['edge_shape_definitions'] = ',\n'.join(all_edge_shape_functions)
        return \
"""Labeled[GraphClass[{
%(edge_lists)s
},
EdgeShapeFunction -> {
%(edge_shape_definitions)s
},
EdgeLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(edge_font_size)d, Bold],
VertexLabelStyle -> Directive[FontFamily -> "CMU Typewriter Text", FontSize -> %(vetex_font_size)d, Bold],
VertexSize -> Large,
VertexLabels -> Placed[Automatic,Center],
GraphLayout -> %(graph_layout_strategy)s,
ImageSize -> {%(width)f, %(height)f}
],"%(graph_name)s"]"""%repl_dict
    
    def draw(self, model, output_dir,FORM_id=None, lmb_id=None):
        """ Outputs the mathematica code for rendering this FORMSuperGraph."""
        
        if FORM_id is not None:
            file_name = 'Graph_%04d'%FORM_id
        else:
            file_name = 'Graph_%s'%self.name
        
        if lmb_id is not None:
            file_name += '_LMB_%04d'%lmb_id

        MM_code = \
"""GraphClass = If[$VersionNumber > 12, EdgeTaggedGraph, Graph];
CreateEdge[u_,v_,t_]:=If[$VersionNumber > 12, DirectedEdge[u, v, t], DirectedEdge[u, v]];
aGraph=%s;
"""%self.get_mathematica_rendering_code(model,FORM_id=FORM_id, lmb_id=lmb_id)
        # Export to PDF in landscape format. One graph per page for now.
        # The 1.2 multiplier accounts for margins
        MM_code += 'Export["%s.pdf", GraphicsGrid[{{aGraph}}], ImageSize -> {%f, %f}];'%(
                    file_name,1.25*self._rendering_size[0],1.25*self._rendering_size[1])
        open(pjoin(output_dir,'%s.m'%file_name),'w').write(MM_code)

    def generate_numerator_form_input(self, additional_overall_factor='', only_algebra=False):
        # create the input file for FORM
        form_diag = self.overall_factor+additional_overall_factor
        for node in self.nodes.values():
            if node['vertex_id'] < 0:
                continue

            form_diag += '*\n vx({},{},{})'.format(
                ','.join(str(p) for p in node['PDGs']),
                ','.join(node['momenta']),
                ','.join(str(i) for i in node['indices'])
            )

        for edge in self.edges.values():
            form_diag += '*\n prop({},{},{},{})'.format(
                edge['PDG'],
                edge['type'],
                edge['momentum'],
                ','.join(str(i) for i in edge['indices'])
            )

        return form_diag

    def get_topo_generator(self, specified_LMB=None):
        """ Returns a topology generator for that FORMSuperGraph."""

        topo_edges = copy.deepcopy(self.edges)

        original_LMB = {}
        external_edges = []
        other_edges = []
        edge_name_to_key = {}
        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in topo_edges.items():
            # Fix for QGRAF pipeline
            if not isinstance(edge_key, tuple):
                edge_key = (*edge_data['vertices'], edge_key) 
          
            if edge_data['type'] == 'virtual':
                if not edge_data['name'].startswith('p'):
                    edge_data['name'] = 'p%s'%edge_data['name']
                other_edges.append((edge_data['name'],edge_key[0],edge_key[1]))
            else:
                if not edge_data['name'].startswith('q'):
                    edge_data['name'] = 'q%s'%edge_data['name'][1:]
                external_edges.append((edge_data['name'],edge_data['vertices'][0],edge_data['vertices'][1]))
            edge_name_to_key[edge_data['name']]=edge_key

            # Test if it is a defining edge of the lmb
            abs_sig = ( [abs(s) for s in edge_data['signature'][0]], [abs(s) for s in edge_data['signature'][1]])
            if sum(abs_sig[0]) == 1 and sum(abs_sig[1]) == 0:
                original_LMB[abs_sig[0].index(1)]=edge_data['name']

        topo_edges = external_edges+other_edges

        # Set the LMB to a sorted one
        original_LMB = sorted(list(original_LMB.items()),key=lambda e: e[0])
        assert(all(oLMBe[0]==i for i,oLMBe in enumerate(original_LMB)))
        original_LMB = [oLMBe[1] for oLMBe in original_LMB]

        topo_generator = LTD.ltd_utils.TopologyGenerator(topo_edges)

        topo_generator.generate_momentum_flow( loop_momenta = (original_LMB if specified_LMB is None else specified_LMB) )
        original_LMB = [edge_name_to_key[oLMBe] for oLMBe in original_LMB]

        return topo_generator, edge_name_to_key, original_LMB

    def sampling_score_lmb(self, lmb, model, edge_key_to_name, edge_name_to_power, edge_name_to_pdg):
        """ Returns a complexity (tuple)which is identical to that of the optimised sampling LMB. This however typically induces a significantly slower generation and runtime of the integrand."""

        lmb_edge_names = [ edge_key_to_name[e_key] for e_key in lmb ]
        lmb_pdg_score = 1
        for edge_name in lmb_edge_names:
            pdg = edge_name_to_pdg[edge_name]
            if model is None:
                if pdg in [21,22]:
                    lmb_pdg_score +=1
            else:
                particle = model.get_particle(pdg)
                if particle.get('spin') == 3 and particle.get('mass').upper()=='ZERO':
                    lmb_pdg_score += 1

        lmb_power_score = 0
        for edge_name in lmb_edge_names:
            lmb_power_score += (edge_name_to_power.get(edge_name,1)-1)

        return (lmb_pdg_score, lmb_power_score)

    def complexity_score_lmb(self, lmb, model, edge_key_to_name, edge_name_to_power, edge_name_to_pdg, weight_combination_rule='product'):
        """ Returns a complexity (tuple) score for this generation lmb. The idea being to pick the generation lmb with highest score so as to improve on generation time."""

        lmb_edge_names = [ edge_key_to_name[e_key] for e_key in lmb ]
        topo_generator, _, _ = self.get_topo_generator(specified_LMB=lmb_edge_names)

        # And adjust all signatures (incl. the string momenta assignment accordingly)
        lmb_edge_names = [ edge_name.replace('q','p') if edge_name.startswith('q') else edge_name[1:] for edge_name in lmb_edge_names]
        signatures = topo_generator.get_signature_map()
        signatures = { edge_name.replace('q','p') if edge_name.startswith('q') else edge_name[1:] : (tuple(sig[0]),
                tuple([ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ])
         ) for edge_name, sig in signatures.items() }
        edge_name_to_power = { edge_name.replace('q','p') if edge_name.startswith('q') else edge_name[1:] : v for edge_name,v in edge_name_to_power.items() }

        # from pprint import pprint, pformat
        # pprint(signatures)
        # print("EDGES")
        # for k, edge in self.edges.items():
        #     print('key=%s, edge=%s'%(str(k),pformat(edge)))
        #     print("New signature: %s"%pformat(signatures[edge['name']]))
        # print("NODES")
        # for k, node in self.nodes.items():
        #     print('key=%s, node=%s'%(str(k),pformat(node)))
        #     print("Edges connected names: %s"%pformat([ self.edges[e_id]['name'] for e_id in node['edge_ids']]))
        #     print("Edges momenta: %s"%pformat([ signatures[self.edges[e_id]['name']] for e_id in node['edge_ids']]))
        # stop

        # Count the powers beyond 1 inside the LMB
        raised_powers_in_LMB = 0

        edge_weights = []
        for edge_key in sorted(self.edges.keys()):
            particle = model.get_particle(abs(self.edges[edge_key]['PDG']))
            # ignore any edge that does not carry a loop momentum
            if all(s==0 for s in signatures[self.edges[edge_key]['name']][0]):
                edge_weights.append(1)
                continue
            sig = list(signatures[self.edges[edge_key]['name']][0])+list(signatures[self.edges[edge_key]['name']][1])
            # Fermion propagators carry momentum and thus have an edge weight (~ number of terms) proportional to the number of
            # terms in the linear decomposition of the momentum carried.            
            if particle.get('spin') == 2:
                edge_weights.append( sum([ 1 if abs(wgt)!=0 else 0 for wgt in sig ]) )
            else:
                edge_weights.append(1)
            
            if self.edges[edge_key]['name'] in lmb_edge_names:
                raised_powers_in_LMB += edge_name_to_power[self.edges[edge_key]['name']]-1

        node_weights = []
        for node_key in sorted(self.nodes.keys()):
            # ignore Pure tree vertices that do not involve any loop momentum (like for example when considering EW decays in QCD computations)
            found_a_loop_momentum = False
            for edge_key in self.nodes[node_key]['edge_ids']:
                if any(s!=0 for s in signatures[self.edges[edge_key]['name']][0]):
                    found_a_loop_momentum = True
                    break
            if not found_a_loop_momentum:
                node_weights.append(1)
                continue
            vertex_spins = tuple(sorted( [ model.get_particle(abs(pdg)).get('spin') for pdg in self.nodes[node_key]['PDGs'] ] ))
            if vertex_spins in [(3,3,3), (1,3,3),(1,1,3)]:
                # The Feynman rule of these vertices typically contain derivatives which can add as many terms as there are components
                # in the momenta of the edges attached
                node_weight = 0
                for edge_key in self.nodes[node_key]['edge_ids']:
                    sig = list(signatures[self.edges[edge_key]['name']][0])+list(signatures[self.edges[edge_key]['name']][1])
                    node_weight += sum([ 1 if abs(wgt)!=0 else 0 for wgt in sig ])
                node_weights.append(node_weight)
            elif vertex_spins in [(3,3,3,3),(2,2,3),(2,2,1)]:
                node_weights.append(1)
            elif all(s==1 for s in vertex_spins):
                node_weights.append(1)
            else:
                raise FormProcessingError("No rule for assigning complexity weight for vertex with spins %s when computing LMB complexitiy score."%(','.join('%d'%s for s in vertex_spins)))

        product_score = functools.reduce(lambda a, b: a*b, edge_weights)*functools.reduce(lambda a, b: a*b, node_weights)
        sum_squared_score = sum(wgt**2 for wgt in edge_weights) + sum(wgt**2 for wgt in node_weights)
        final_score = None
        # We use the negative value as we want to minimize the complexity, and entries beyond the first in the tuple are used for tie-breaker
        # We always place the "raised_powers_in_LMB" first because we always seek to maximise the number of raised propagator edges part of the generation LMB selected.
        if weight_combination_rule == 'product':
            final_score = (raised_powers_in_LMB, -product_score, -sum_squared_score)
        elif weight_combination_rule == 'sum_squared':
            final_score = (raised_powers_in_LMB, -sum_squared_score, -product_score)
        elif weight_combination_rule == 'combination':
            final_score = (raised_powers_in_LMB, -product_score-sum_squared_score,-product_score, -sum_squared_score)
        else:
            raise FormProcessingError("Function complexity_score_lmb only supports the combination rules 'product' and 'sum_squared' in the options.")

        return final_score

    def score_lmb(self, *args):

        # Choose which scoring algorithm to use
        #return self.sampling_score_lmb(*args)
        #return self.complexity_score_lmb(*args, weight_combination_rule='sum_squared')
        #return self.complexity_score_lmb(*args, weight_combination_rule='product')
        return self.complexity_score_lmb(*args, weight_combination_rule='combination')

    def adjust_LMBs(self, model):
        """ Depending on the FORM options 'number_of_lmbs' and 'reference_lmb', this function fills in the attribute 
        additional_lmbs of this class."""

        if (FORM_processing_options['reference_lmb'] is None) and (FORM_processing_options['number_of_lmbs'] is None) and not FORM_processing_options['optimise_generation_lmb']:
            return
    
        topo_generator, edge_name_to_key, original_LMB = self.get_topo_generator()
        edge_key_to_name = {v:k for k,v in edge_name_to_key.items()}
        
        # Get propagator powers
        sig_map=topo_generator.get_signature_map()
        sig_powers = {}
        for edge_name, sig in sig_map.items():
            a = (tuple(sig[0]),tuple(sig[1]))
            b =  (tuple([-s for s in sig[0]]),tuple([-s for s in sig[1]]))
            for s in [a,b]:
                if s in sig_powers:
                    sig_powers[s] += 1
                else:
                    sig_powers[s] = 1
        edge_name_to_power = { e_name: sig_powers[(tuple(sig_map[e_name][0]),tuple(sig_map[e_name][1]))] for e_name in edge_name_to_key } 

        all_lmbs = topo_generator.loop_momentum_bases()
        all_lmbs= [ tuple([edge_name_to_key[topo_generator.edge_map_lin[e][0]] for e in lmb]) for lmb in all_lmbs]

        edge_name_to_pdg = {edge['name']:edge['PDG'] for edge in self.edges.values()}
        # add the prefix p
        edge_name_to_pdg= {'p%s'%k: v for k,v in edge_name_to_pdg.items() if not k.startswith('p') }

        forced_LMB_index = None
        if FORM_processing_options['optimise_generation_lmb']:
            
            lmb_metric = []
            for (lmb_index, lmb) in enumerate(all_lmbs):
                score = self.score_lmb(lmb, model, edge_key_to_name, edge_name_to_power, edge_name_to_pdg)
                lmb_metric.append(tuple(list(score)+[lmb_index,]))

            lmb_metric.sort(key=lambda score:tuple(score[:-1]), reverse=True)
            # Convenient printout for debugging heuristics for choice of basis for particular supergraphs
            # if self.name == 'SG_QG0':
            #     print("All lmb scores for %s:\n%s"%(self.name,'\n'.join(
            #         '#%-2d: %s = %s'%(lmb_info[-1],','.join(edge_key_to_name[e] for e in all_lmbs[lmb_info[-1]]),','.join('%-5s'%s for s in lmb_info[:-1])) for lmb_info in lmb_metric
            #     )))
            # stop
            forced_LMB_index = lmb_metric[0][-1]

        reference_lmb = None
        if isinstance(FORM_processing_options['reference_lmb'], int):
            reference_lmb = FORM_processing_options['reference_lmb']
        elif isinstance(FORM_processing_options['reference_lmb'],dict) and self.name in FORM_processing_options['reference_lmb']:
            reference_lmb = 1
            all_lmbs = [ tuple([ edge_name_to_key[e] for e in FORM_processing_options['reference_lmb'][self.name] ]), ]

        # Then overwrite the reference LMB if the user requested it
        if (reference_lmb is not None) or (forced_LMB_index is not None):
            if reference_lmb is not None:
                original_LMB = all_lmbs[(reference_lmb-1)%len(all_lmbs)]
            else:
                original_LMB = all_lmbs[forced_LMB_index]

            # Sort the LMB according to PDGs. Only do it if the user did not specify the list of edges explicitly
            if not (isinstance(FORM_processing_options['reference_lmb'],dict) and self.name in FORM_processing_options['reference_lmb']):
                original_LMB=tuple(sorted(original_LMB,key=lambda e_key: abs(edge_name_to_pdg[edge_key_to_name[e_key]])))

            # Regenerate the topology with this new overwritten LMB
            topo_generator, _, _ = self.get_topo_generator(specified_LMB=[ edge_key_to_name[e_key] for e_key in original_LMB ])
            # And adjust all signatures (incl. the string momenta assignment accordingly)
            signatures = topo_generator.get_signature_map()
            signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in signatures.items() }
            for edge_key, edge_data in self.edges.items():
                # Fix for QGRAF pipeline
                if not isinstance(edge_key, tuple):
                    edge_key = (*edge_data['vertices'], edge_key) 
          
                edge_data['signature'] = signatures[edge_key]
                edge_data['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(edge_data['signature'], set_outgoing_equal_to_incoming=False)
            for node_key, node_data in self.nodes.items():
                node_data['momenta'] = tuple([                    
                        FORMSuperGraph.momenta_decomposition_to_string(
                            [ 
                                [ s*(1 if self.edges[e_key]['vertices'][1]==node_key else -1) for s in self.edges[e_key]['signature'][0] ],
                                [ s*(1 if self.edges[e_key]['vertices'][1]==node_key else -1) for s in self.edges[e_key]['signature'][1] ],
                            ],
                            set_outgoing_equal_to_incoming=False)
                    for e_key in node_data['edge_ids']])

        if FORM_processing_options['number_of_lmbs'] is None:
            return

        original_lmb_signatures = topo_generator.get_signature_map()
        original_lmb_signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in original_lmb_signatures.items() }

        # Now generate copies of this supergraph with different LMBs
        additional_lmbs_SGs = []
        for i_lmb, lmb in enumerate(all_lmbs):
            if lmb==original_LMB:
                continue
            if FORM_processing_options['number_of_lmbs']>=0 and len(additional_lmbs_SGs)==FORM_processing_options['number_of_lmbs']:
                break

            # Generate the topology with this additional LMB
            other_lmb_topo_generator, _, _ = self.get_topo_generator(specified_LMB=[ edge_key_to_name[e_key] for e_key in lmb ])
            # And adjust all signatures (incl. the string momenta assignment accordingly)
            other_lmb_signatures = other_lmb_topo_generator.get_signature_map()
            other_lmb_signatures = { edge_name_to_key[edge_name]: [sig[0],
                    [ i+o for i,o in zip(sig[1][:len(sig[1])//2],sig[1][len(sig[1])//2:]) ]
                ] for edge_name, sig in other_lmb_signatures.items() }

            # Compute the affine transformation to go from the original LMB to this additional LMB
            affine_transfo_other_lmb_to_original = []
            for defining_edge_key in lmb:
                other_sig = original_lmb_signatures[defining_edge_key]
                affine_transfo_other_lmb_to_original.append( [ list(other_sig[0]), list(other_sig[1]) ]  )
            affine_transfo_original_lmb_to_other_lmb = []
            for defining_edge_key in original_LMB:
                other_sig = other_lmb_signatures[defining_edge_key]
                affine_transfo_original_lmb_to_other_lmb.append( [ list(other_sig[0]), list(other_sig[1]) ]  )

            other_LMB_super_graph = copy.deepcopy(self)
            # Flag these additional supergraphs as *NOT* original representative by setting their 'additional_lmbs' to an integer instead
            # which identifies which additional LMB it corresponds to.
            other_LMB_super_graph.additional_lmbs = i_lmb+1
            for edge_key, edge_data in other_LMB_super_graph.edges.items():
                # Fix for QGRAF pipeline
                if not isinstance(edge_key, tuple):
                    edge_key = (*edge_data['vertices'], edge_key) 
          
                edge_data['signature'] = other_lmb_signatures[edge_key]
                edge_data['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(edge_data['signature'], set_outgoing_equal_to_incoming=False)
            for node_key, node_data in other_LMB_super_graph.nodes.items():
                node_data['momenta'] = tuple([
                    FORMSuperGraph.momenta_decomposition_to_string(
                        [ 
                            [ s*(1 if other_LMB_super_graph.edges[e_key]['vertices'][1]==node_key else -1) for s in other_LMB_super_graph.edges[e_key]['signature'][0] ],
                            [ s*(1 if other_LMB_super_graph.edges[e_key]['vertices'][1]==node_key else -1) for s in other_LMB_super_graph.edges[e_key]['signature'][1] ],
                        ],
                        set_outgoing_equal_to_incoming=False)
                     for e_key in node_data['edge_ids']])

            additional_lmbs_SGs.append( (i_lmb+1, affine_transfo_other_lmb_to_original, affine_transfo_original_lmb_to_other_lmb, other_LMB_super_graph ) )

        self.additional_lmbs = additional_lmbs_SGs

    def derive_signatures(self):
        n_incoming = sum([1 for edge in self.edges.values() if edge['type'] == 'in'])
        n_loops = len(self.edges) - len(self.nodes) + 1

        # parse momentum
        p = re.compile(r'(^|\+|-)(k|p)(\d*)')
        for edge in self.edges.values():
            parsed = [e.groups() for e in p.finditer(edge['momentum'])]
            signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
            for pa in parsed:
                if pa[1] == 'p':
                    signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                else:
                    signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1

            edge['signature'] = signature

    def impose_signatures(self):

        for eid, e in self.edges.items():
            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
            neg_mom =  FORMSuperGraph.momenta_decomposition_to_string([[-s for s in sp] for sp in e['signature']], set_outgoing_equal_to_incoming=False)
            for i, vi in enumerate(e['vertices']):
                e_index = self.nodes[vi]['edge_ids'].index(eid)
                mom = list(self.nodes[vi]['momenta'])
                mom[e_index] = neg_mom if i == 0 else e['momentum']
                self.nodes[vi]['momenta'] = mom

    @classmethod
    def sort_edges(cls, model, edges_to_sort):
        """ Sort the edges adjacent to a node so as to map the conventions
        from the Feynmna rules encoded in Form."""

        # The only relevant properties are the spin and wether an edge
        # is a particle or antiparticle.
        edge_particles = [model.get_particle(e['PDG']) for e in edges_to_sort]
        identities = [ ( p.get('spin') * (1 if p.get('is_part') else -1) , position ) for position, p in enumerate(edge_particles) ]
        identities.sort(key=lambda el: el[0])
        canonical_identifier = tuple([identity for identity, position in identities])
        
        # Always put in the trivial rules for multiscalar n-point vertex
        if all(i==1 for i in canonical_identifier):
            new_position = list(range(len(canonical_identifier)))
        elif canonical_identifier not in cls._FORM_Feynman_rules_conventions:
            raise FormProcessingError("Conventions for FORM Feynman rules of signature {} not specifed.".format(canonical_identifier))
        else:
            new_position = cls._FORM_Feynman_rules_conventions[canonical_identifier]

        return [ edges_to_sort[identities[position][1]] for position in new_position ]

    @classmethod
    def momenta_decomposition_to_string(cls, momenta_decomposition, set_outgoing_equal_to_incoming=True):
        """ Turns ((1,0,0,-1),(1,1)) into 'k1-k4+p1+p2'"""

        res = ""
        first=True
        # The outgoing momenta are set element-wise equal to the incoming ones.
        if set_outgoing_equal_to_incoming:
            momenta_decomposition = [
                momenta_decomposition[0],
                [inp+outp for inp, outp in zip(
                    momenta_decomposition[1][:len(momenta_decomposition[1])//2],
                    momenta_decomposition[1][len(momenta_decomposition[1])//2:]
                )]
            ]
        # Fuse outgoing and incoming
        for symbol, mom_decomposition in zip(('k','p'),momenta_decomposition):
            for i_k, wgt in enumerate(mom_decomposition):
                if wgt!=0:
                    if first:
                        if wgt<0:
                            res+="-"
                        first=False
                    else:
                        if wgt<0:
                            res+="-"
                        else:
                            res+="+"
                    if abs(wgt)!=1:
                        res+="%d*"%abs(wgt)
                    res+="%s%d"%(symbol,(i_k+1))

        return res

    @classmethod
    def from_LTD2SuperGraph(cls, LTD2_super_graph):
        """ Creates a FORMSuperGraph from an LTD2_super_graph."""

        # First make sure the momentum routing of the LTD2_super_graph
        # is set:
        LTD2_super_graph.set_momentum_routing()

        overall_factor = "1"

        if LTD2_super_graph.symmetry_factor is not None:
            multiplicity = LTD2_super_graph.symmetry_factor
        else:
            multiplicity = 1

        # Let us just match the overall phase picked for the MG num:
        overall_phase = complex(-1.0,0.0)**len(LTD2_super_graph.cuts)
        if overall_phase.imag != 0:
            raise FormProcessingError("No support for overall complex phase yet (Ben: how do we put a complex number in FORM? ^^)")
        else:
            overall_factor += '*%d'%int(overall_phase.real)

        fermion_factor = LTD2_super_graph.diag_left_of_cut.fermion_factor*LTD2_super_graph.diag_right_of_cut.fermion_factor
        overall_factor += '*%d'%fermion_factor

        model = LTD2_super_graph.model

        # Let us also include a factor -1 for each closed ghost loop.
        # This is automatically done in MadGraph already by considering the "wavefunction" of external scalar ghosts to be sqrt(i).
        # Said differently, we want a factor -1 for each pair of two ghosts in the final state
        n_ghosts = len([ 1 for c in LTD2_super_graph.cuts if model.get_particle(LTD2_super_graph.graph.edges[c[1]]['pdg']).get('ghost') ])
        assert(n_ghosts%2==0)
        if n_ghosts > 0:
            overall_factor += '*%d'%(-1**(n_ghosts//2))

        local_graph = copy.deepcopy(LTD2_super_graph.graph)
        
        # Collect the outer-most nodes
        outer_in_nodes = []
        outer_out_nodes = []
        for node_key, node_data in local_graph.nodes.items():
            # Assign the corresponding edge
            if node_key.startswith('I') or node_key.startswith('O'):
                for edge_key, edge_data in local_graph.edges.items():
                    if edge_key[0]==node_key or edge_key[1]==node_key:
                        node_data['edge_ids'] = tuple([edge_key,])
                        break
                if 'edge_ids' not in node_data:
                    node_data['edge_ids'] = tuple([])
            if node_key.startswith('I'):
                outer_in_nodes.append(node_key)
                node_data['momenta'] = ('p%d'%int(node_key[1:]),)
                node_data['vertex_id'] = -1 # Incoming node at assigned vertex ID = -1
            if node_key.startswith('O'):
                outer_out_nodes.append(node_key)
                node_data['momenta'] = ('p%d'%int(node_key[1:]),)
                node_data['vertex_id'] = -2 # Incoming node at assigned vertex ID = -2

        curr_index = 0
        # Assign indices to the outermost nodes
        for node_key in outer_in_nodes:
            this_index = int(node_key[1:])
            curr_index = max(curr_index,this_index)
            local_graph.nodes[node_key]['indices'] = (this_index,)

        n_incoming=int(curr_index)
        for node_key in outer_out_nodes:
            this_index = n_incoming+int(node_key[1:])
            curr_index = max(curr_index,this_index)
            local_graph.nodes[node_key]['indices'] = (this_index,)

        # assign indices to all edges now
        for edge_key, edge_data in local_graph.edges.items():
            edge_indices = []
            found_external_node = False
            for node_key in [edge_key[0],edge_key[1]]:
                if node_key in outer_in_nodes+outer_out_nodes:
                    edge_indices = [ local_graph.nodes[node_key]['indices'][0], ]
                    found_external_node = True
                    local_graph.nodes[node_key]['PDGs'] = (edge_data['pdg'],)
                else:
                    curr_index += 1
                    if not found_external_node:
                        edge_indices.append(curr_index)
            edge_data['indices'] = tuple(edge_indices)
            if edge_key[0] in outer_in_nodes or edge_key[1] in outer_in_nodes:
                edge_data['type'] = 'in'
#                loop_momenta = tuple([0 for _ in range(LTD2_super_graph.n_loops)])
#                external_momenta = [0 for _ in range(len(outer_in_nodes))]
#                if edge_key[0] in outer_in_nodes:
#                    external_momenta[int(edge_key[0][1:])-1] = 1
#                elif edge_key[1] in outer_in_nodes:
#                    external_momenta[int(edge_key[1][1:])-1] = -1
#                edge_data['momentum'] = (loop_momenta, tuple(external_momenta))
            elif edge_key[0] in outer_out_nodes or edge_key[1] in outer_out_nodes:
                edge_data['type'] = 'out'
#                loop_momenta = tuple([0 for _ in range(LTD2_super_graph.n_loops)])
#                external_momenta = [0 for _ in range(len(outer_in_nodes))]
#                if edge_key[1] in outer_out_nodes:
#                    external_momenta[int(edge_key[1][1:])-1] = 1
#                elif edge_key[0] in outer_out_nodes:
#                    external_momenta[int(edge_key[0][1:])-1] = -1
#                edge_data['momentum'] = (loop_momenta, tuple(external_momenta))
            else:
                edge_data['type'] = 'virtual'
            edge_data['PDG'] = edge_data['pdg']
            del edge_data['pdg']
            # Keep the name (used for rendering only)
#            del edge_data['name']

        # assign node extra information
        for node_key, node_data in local_graph.nodes.items():
            # Skip external nodes that have already been treated
            if node_key in outer_in_nodes+outer_out_nodes:
                continue
            # First collect all adjacent edges:
            adjacent_in_edges = [
                dict(list(in_edge_data.items())+[('key',(u,v,c))]) for u,v,c,in_edge_data in 
                local_graph.in_edges(node_key,data=True,keys=True)
            ]
            adjacent_out_edges = [
                dict(list(out_edge_data.items())+[('key',(u,v,c))]) for u,v,c,out_edge_data in 
                local_graph.out_edges(node_key,data=True,keys=True)
            ]
            # The direction only matters in so far as we must flip the
            # momentum carried by edges that are outgoing as well as
            # change their PDG from particle to antiparticle
            for data in adjacent_in_edges:
                if len(data['indices'])==2:
                    data['index_for_this_node'] =  data['indices'][1]
                else:
                    # This is an external edge with a single index.
                    data['index_for_this_node'] =  data['indices'][0]
            for data in adjacent_out_edges:
                data['index_for_this_node'] =  data['indices'][0]
                data['momentum'] = (
                    tuple([-p for p in data['momentum'][0]]),
                    tuple([-p for p in data['momentum'][1]])
                )
                data['PDG'] = model.get_particle(data['PDG']).get_anti_pdg_code()

            # We must now sort the edges according to the special rules for the vertex Feynman rules
            all_adjacent_edges = cls.sort_edges(model,adjacent_in_edges+adjacent_out_edges)
            node_data['PDGs'] = tuple([e['PDG'] for e in all_adjacent_edges])
            node_data['indices'] = tuple([e['index_for_this_node'] for e in all_adjacent_edges])
            node_data['momenta'] = tuple([cls.momenta_decomposition_to_string(e['momentum']) for e in all_adjacent_edges])
            node_data['edge_ids'] = tuple([e['key'] for e in all_adjacent_edges])

            # Example printout information about the vertex
            #misc.sprint("Vertex of node %s:\n%s"%(
            #        str(node_key),pformat(model.get_interaction(node_data['vertex_id']))))

        # Finally overwrite the edge momentum so as to be a string
        for edge_key, edge_data in local_graph.edges.items():
            edge_data['momentum'] =cls.momenta_decomposition_to_string(edge_data['momentum'])
            edge_data['vertices'] =(edge_key[0],edge_key[1])
            # In FORM conventions the fermion are going together with their flow, so we need
            # to flip the order of their fundamental/antifundamental indices so that FORM
            # builds the correct propagator. 
            # NO LONGER NEEDED: This now directly done right in FORM.
#            particle = model.get_particle(edge_data['PDG'])
#            if len(edge_data['indices'])>1 and (particle.get('spin')%2==0 or particle.get('ghost')):
#                if not particle.get('is_part'):
#                    edge_data['indices'] = tuple([edge_data['indices'][1],edge_data['indices'][0]])

        graph_name = 'P%(proc_id)dL%(left_diagram_id)dR%(right_diagram_id)d'%LTD2_super_graph.call_signature
        form_super_graph =  cls(
            name = graph_name,
            call_identifier=LTD2_super_graph.call_signature,
            edges = dict(local_graph.edges),
            nodes = dict(local_graph.nodes),
            overall_factor = overall_factor,
            multiplicity = multiplicity
        )
        #misc.sprint(graph_name)
        #misc.sprint(pformat(dict(LTD2_super_graph.graph.edges.items())))
        #misc.sprint(pformat(dict(local_graph.edges.items())))

        return form_super_graph

    @classmethod
    def from_dict(cls, file_path):
        """ Creates a FORMSuperGraph from a Python dict file path."""

        # TODO: Creates an instance from a Python dict dump.
        pass

    def to_dict(self, file_path=None):
        """ Outputs the FORMSuperGraph self to a Python dict file path."""

        dict_to_dump = {
            'edges' : self.edges,
            'nodes' : self.nodes,
            'overall_factor' : self.overall_factor,
            'multiplicity' : self.multiplicity
        }
        if file_path:
            open(file_path,'w').write(pformat(dict_to_dump))
        else:
            return dict_to_dump

    def generate_ltd_integrand(self, g, graph_id, signature_offset, unique_ltd, constants, unique_propagators, prop_id, prop_mom_in_lmb, energies, on_shell_condition):
        diagres = []
        propagators = []
        res = []

        # enumerate all cut options
        for co in [x for css in g.ltd_cut_structure for x in 
                product(*[[(cs, (li, i)) for i in range(len(l.propagators))] for li, (cs, l) in enumerate(zip(css, g.loop_lines)) if cs != 0])]:
            ltd_closure_factor = int(np.prod([s for (s, _) in co]))

            # construct the cut basis to LTD loop momentum basis mapping, used to substitute the numerator
            mat = [g.loop_lines[li].signature for (_, (li, _)) in co]
            if mat == []:
                nmi = []
                cb_to_lmb = []
            else:
                nmi = np.linalg.inv(np.array(mat).transpose()) # the tranpose matrix is used for signatures
                cb_to_lmb = np.linalg.inv(np.array(mat))

            m = []
            ltdenergy = ['ltd{0},{1}E{0}'.format(prop_id[(li, pi)], '+' if cut_sign == 1 else '-') for (cut_sign, (li, pi)) in co]
            for i, r in enumerate(cb_to_lmb):
                mm = []
                for (c, (_, (li, pi))) in zip(r, co):
                    if c != 0:
                        ext = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(li, pi)]][1], prop_mom_in_lmb[prop_id[(li, pi)]][2]), True)
                        if ext == '':
                            ext = '0'
                        
                        mm.append('{}ltd{}{}energies({})'.format('' if c == 1 else '-', prop_id[(li, pi)], '-' if c == 1 else '+', ext))
                m += ['fmb{}'.format(i + signature_offset + 1), '+'.join(mm)]

            r = []
            der = []
            for li, l in enumerate(g.loop_lines):
                if all(s == 0 for s in l.signature):
                    continue

                energy = []
                energy_full = ''
                sig_map = nmi.dot(l.signature)
                for (sig_sign, (cut_sign, (lci, pci))) in zip(sig_map, co):
                    if sig_sign != 0:
                        sig = tuple(list(prop_mom_in_lmb[prop_id[(lci, pci)]][1]) + list(prop_mom_in_lmb[prop_id[(lci, pci)]][2]))
                        minsig = tuple(list(-prop_mom_in_lmb[prop_id[(lci, pci)]][1]) + list(-prop_mom_in_lmb[prop_id[(lci, pci)]][2]))
                        
                        if sig in on_shell_condition:
                            momp = '+{}'.format(on_shell_condition[sig])
                        elif minsig in on_shell_condition:
                            momp = '-{}'.format(on_shell_condition[minsig])
                        else:
                            momp = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(lci, pci)]][1], prop_mom_in_lmb[prop_id[(lci, pci)]][2]), True)
                            momp = '0' if momp == '' else '+energies({})'.format(momp)

                        energy.append('{},ltd{},{}{}'.format(int(sig_sign), prop_id[(lci, pci)], '+' if -sig_sign == 1 else '-', momp))
                        # the full energy including the cut sign
                        energy_full += '{}E{}{}{}'.format('+' if int(cut_sign * sig_sign) == 1 else '-', prop_id[(lci, pci)], '+' if -int(sig_sign) == 1 else '-', momp)

                for pi, p in enumerate(l.propagators):
                    powmod = '' if p.power == 1 else '^' +  str(p.power)
                    if (1, (li, pi)) in co:
                        prop_plus = '2*E{0}'.format(prop_id[(li, pi)])
                        if prop_plus not in unique_propagators:
                            unique_propagators[prop_plus] = len(unique_propagators)
                        propagators.append('invd{},{}'.format(unique_propagators[prop_plus], prop_plus))
                        r.append('prop(invd{0},1,ltd{1},0,E{1}){2}'.format(unique_propagators[prop_plus], prop_id[(li, pi)], powmod))
                        if p.power > 1:
                            # add derivative prescription
                            der.append('ltd{},{}'.format(prop_id[(li, pi)], p.power - 1))
                    elif (-1, (li, pi)) in co:
                        prop_min = '-2*E{0}'.format(prop_id[(li, pi)])
                        if prop_min not in unique_propagators:
                            unique_propagators[prop_min] = len(unique_propagators)
                        propagators.append('invd{},{}'.format(unique_propagators[prop_min], prop_min))
                        r.append('prop(invd{0},-1,ltd{1},0,-E{1}){2}'.format(unique_propagators[prop_min], prop_id[(li, pi)], powmod))
                        if p.power > 1:
                            der.append('ltd{},{}'.format(prop_id[(li, pi)], p.power - 1))
                    else:
                        sig = tuple(list(prop_mom_in_lmb[prop_id[(li, pi)]][1]) + list(prop_mom_in_lmb[prop_id[(li, pi)]][2]))
                        minsig = tuple(list(-prop_mom_in_lmb[prop_id[(li, pi)]][1]) + list(-prop_mom_in_lmb[prop_id[(li, pi)]][2]))
                        
                        if sig in on_shell_condition:
                            momp = '+{}'.format(on_shell_condition[sig])
                        elif minsig in on_shell_condition:
                            momp = '-{}'.format(on_shell_condition[minsig])
                        else:
                            momp = self.momenta_decomposition_to_string((prop_mom_in_lmb[prop_id[(li, pi)]][1], prop_mom_in_lmb[prop_id[(li, pi)]][2]), True)
                            momp = '' if momp == '' else '+energies({})'.format(momp)

                        prop_plus = '{}+E{}{}'.format(energy_full, prop_id[(li, pi)], momp)
                        prop_min = '{}-E{}{}'.format(energy_full, prop_id[(li, pi)], momp)

                        if prop_plus not in unique_propagators:
                            unique_propagators[prop_plus] = len(unique_propagators)
                        if prop_min not in unique_propagators:
                            unique_propagators[prop_min] = len(unique_propagators)

                        propagators.append('invd{},{}'.format(unique_propagators[prop_plus], prop_plus))
                        propagators.append('invd{},{}'.format(unique_propagators[prop_min], prop_min))

                        r.append('prop(invd{0},{1},E{2}{3}){4}*prop(invd{5},{1},-E{2}{3}){4}'.format(unique_propagators[prop_plus], ','.join(energy), prop_id[(li, pi)],
                            momp, powmod, unique_propagators[prop_min]))

            diagres.append('{}*{}{}{}({})'.format(
                ltd_closure_factor,
                'ltdcbtolmb({})*'.format(','.join(m)) if len(m) > 0 else '',
                'ltdenergy({})*'.format(','.join(ltdenergy)) if len(ltdenergy) > 0 else '',
                'der({})*'.format(','.join(der))  if len(der) > 0 else '',
                '*'.join(r) if len(r) > 0 else '1'))

        res.append('\n\t\t\t+'.join(diagres))

        res = '\n\t\t*'.join(['({})'.format(l) for l in res])

        propagators = list(sorted(set(propagators)))

        ltd_instr = '(-1)^{}*(2*pi*i_)^{}*constants({})*\n\tellipsoids({})*\n\tallenergies({})*(\n\t\t{}\n)'.format(g.n_loops,
            g.n_loops, ','.join(constants + ['1']), ','.join(propagators), ','.join(energies), res)

        if ltd_instr not in unique_ltd:
            unique_ltd[ltd_instr] = graph_id
            integrand_body = 'Fill ltdtopo({},{}) = {};\n'.format(*graph_id, ltd_instr)
        else:
            integrand_body = ''

        return (integrand_body, 'Fill ltdmap({},{}) = diag({},{});\n'.format(*graph_id, *unique_ltd[ltd_instr]))

    def generate_integrand(self, topo, integrand_type, workspace, numerator_call, progress_bar=None):
        """Construct a table of integrand descriptors"""
        unique_pf = {}
        integrand_body = ''
        topo_map = {}

        unique_ltd = {}
        ltd_topo_map = ''
        ltd_integrand_body = ''
        den_library = []

        for cut in topo.cuts:
            unique_energies = {}
            unique_propagators = {}
            for diag_set in cut['diagram_sets']:
                # collect all graphs
                graphs = []
                for diag_info in diag_set['diagram_info']:
                    for uv_structure in diag_info['uv']:
                        signature_offset = 0 # offset in the forest basis
                        for uv_subgraph in uv_structure['uv_subgraphs']:
                            if uv_subgraph['first_occurrence_id'] != uv_subgraph['id']:
                                signature_offset += uv_subgraph['derived_graphs'][0]['loop_topo'].n_loops
                                continue

                            for dg in uv_subgraph['derived_graphs']:
                                if dg['skip_pf']:
                                    continue

                                graphs.append((signature_offset, dg['id'], dg['loop_topo'], None))
                                if 'soft_ct_id' in dg:
                                   massct = (False, uv_subgraph['onshell']) if uv_subgraph['onshell'] is not None else None
                                   graphs.append((signature_offset, dg['soft_ct_id'], dg['loop_topo_orig_mass'], massct))
                                if 'onshell_ct_id' in dg:
                                   massct = (True, uv_subgraph['onshell']) if uv_subgraph['onshell'] is not None else None
                                   graphs.append((signature_offset, dg['onshell_ct_id'], dg['loop_topo_orig_mass'], massct))
                            graphs.append((signature_offset, uv_subgraph['integrated_ct_id'], uv_subgraph['integrated_ct_bubble_graph'], None))
                            signature_offset += uv_subgraph['derived_graphs'][0]['loop_topo'].n_loops
                        if len(uv_structure['bubble']) > 0:
                            for bubble in uv_structure['bubble']:
                                graphs.append((signature_offset, bubble['id'], bubble['remaining_graph_loop_topo'], None))
                        else:
                            graphs.append((signature_offset, uv_structure['remaining_graph_id'], uv_structure['remaining_graph_loop_topo'], None))

                for signature_offset, graph_id, g, onshell in graphs:
                    signatures, n_props, energies, constants, shift_map, unique_energy = [], [], [], [], [], []
                    pf_prefactor = ['1']
                    prop_mom_in_lmb = {}
                    prop_id = {}

                    on_shell_condition = {}
                    if onshell is not None:
                        flip_mass_sign, self_energy_ext_edges = onshell

                        # construct the on-shell condition, TODO: check sign
                        ext_edge = next(ee for ee in self.edges.values() if ee['name'] == self_energy_ext_edges[0])
                        on_shell_condition[tuple(ext_edge['signature'][0] +
                            ext_edge['signature'][1] + [0]*len( ext_edge['signature'][1]))] = '{}masses({})'.format('-' if flip_mass_sign else '', ext_edge['PDG'])
                    for li, l in enumerate(g.loop_lines):
                        is_constant = all(s == 0 for s in l.signature)
                        if not is_constant:
                            signatures.append(list(l.signature))
                            n_props.append(sum(p.power for p in (l.propagators)))
                        for pi, p in enumerate(l.propagators):
                            # contruct the momentum in the LMB
                            # it should be the same as self.edge['momentum'] apart from the UV and the sign
                            lmp = np.array([0]*topo.topo.n_loops)
                            for s, v in zip(l.signature, g.loop_momentum_map):
                                lmp += s * np.array(v[0])

                            if onshell is None:
                                # transport shift to LMB
                                shift = np.array([0]*topo.topo.n_loops)
                                extshift = np.array(p.parametric_shift[1])
                                for s, c in zip(p.parametric_shift[0], cut['cuts']):
                                    shift += s * np.array(c['signature'][0])
                                    extshift += s * np.array(c['signature'][1])
                                extshift = np.array(list(extshift[:len(extshift)//2]) + [0]*(len(extshift)//2)) +\
                                    np.array(list(extshift[len(extshift)//2:]) + [0]*(len(extshift)//2))
                                totalmom = self.momenta_decomposition_to_string((lmp + shift, extshift), True)
                            else:
                                # for on-shell graphs, the shifts are still in the graph basis
                                ext_edge_sigs = [next(ee for ee in self.edges.values() if ee['name'] == e)['signature'] for e in self_energy_ext_edges]
                                shift = np.array([0]*topo.topo.n_loops)
                                extshift = np.array([0]*len(ext_edge_sigs[0][1]))
                                for s, c in zip(p.parametric_shift[1], ext_edge_sigs):
                                    shift += s * np.array(c[0])
                                    extshift += s * np.array(c[1])
                                extshift = np.array(list(extshift) + [0]*len(extshift))
                                # drop the spatial shift
                                totalmom = self.momenta_decomposition_to_string((lmp, np.array([0]*len(ext_edge_sigs[0][1]))), True)

                            # recycle energy computations when there are duplicate edges
                            edge_mass = 'masses({})'.format(next(ee for ee in self.edges.values() if ee['name'] == p.name)['PDG'])
                            energy_instr = '{},{}'.format(totalmom, edge_mass if not p.uv else 'mUV')
                            if not is_constant:
                                if energy_instr not in unique_energies:
                                    unique_energies[energy_instr] = len(unique_energies)

                                pf_prefactor.append('2*E{}'.format(unique_energies[energy_instr]) if p.power == 1 else '(2*E{})^{}'.format(unique_energies[energy_instr], p.power))
                                e_str = 'E{},{}'.format(unique_energies[energy_instr], energy_instr)
                                if e_str not in energies:
                                    energies.append(e_str)
                                prop_id[(li, pi)] = unique_energies[energy_instr]
                                prop_mom_in_lmb[unique_energies[energy_instr]] = (lmp, shift, extshift)

                                for _ in range(p.power):
                                    shift_map.append(list(shift) + list(extshift))
                                    unique_energy.append(unique_energies[energy_instr])
                            else:
                                for _ in range(p.power):
                                    constants.append(energy_instr)

                    if integrand_type == "both" or integrand_type == "PF":
                        pf_prefactor = '*'.join(pf_prefactor)

                        if len(signatures) == 0:
                            # no loop dependence for this cut
                            res = '\t1\n'
                            resden = ''
                        else:
                            pf = LTD.partial_fractioning.PartialFractioning(n_props, signatures,
                                                    name=str(diag_set['id']), shift_map=np.array(shift_map).T,
                                                    n_sg_loops=topo.topo.n_loops, ltd_index=signature_offset,
                                                    progress_bar = progress_bar)
                            pf.shifts_to_externals()
                            res, used_props = pf.to_FORM(energy_index_map=unique_energy, den_library=den_library, on_shell_conditions=on_shell_condition)
                            res = '\n'.join(['\t' + l for l in res.split('\n')])
                            resden = ','.join('invd{},{}'.format(i, d) for i, d in enumerate(den_library) if i in used_props)

                        pf_instr = '(2*pi*i_)^{}*constants({})*\nallenergies({})*\nellipsoids({})*(\n{})'.format(g.n_loops, ','.join(constants + [pf_prefactor]),
                            ','.join(energies), resden, res)

                        if pf_instr not in unique_pf:
                            integrand_body += 'Fill pftopo({}) = {};\n'.format(len(unique_pf), pf_instr)
                            unique_pf[pf_instr] = len(unique_pf)

                        topo_map[graph_id] = unique_pf[pf_instr]

                    if integrand_type == "both" or integrand_type == "LTD":
                        (body_addition, topo_addition) = self.generate_ltd_integrand(g, graph_id, signature_offset, unique_ltd, constants, unique_propagators, prop_id, prop_mom_in_lmb,
                                energies, on_shell_condition)
                        ltd_integrand_body += body_addition
                        ltd_topo_map += topo_addition

            # the on-shell energy needs to be accessible to construct the input file
            diag_set['energies'] = unique_energies
        topo.topo_map = topo_map

        if integrand_body != 0:
            with open(pjoin(workspace, 'pftable_{}.h'.format(numerator_call)), 'w') as f:
                f.write("""
Auto S invd, E, shift, ltd;
S r, s;
CF a, num, ncmd, conf1, ellipsoids, allenergies, replace, constants;
NF energies;
CTable pftopo(0:{});

{}
""".format(len(unique_pf), integrand_body))

        if ltd_integrand_body != 0:
            with open(pjoin(workspace, 'ltdtable_{}.h'.format(numerator_call)), 'w') as f:
                f.write("""
Auto S invd, E, shift, ltd;
S r, s;
CF a, num, ncmd, conf1, replace, energies, ellipsoids, ltdcbtolmb, ltdenergy, constants;
NF allenergies;
CTable ltdtopo(0:{},0:{});
CTable ltdmap(0:{},0:{});

{}

{}
""".format(0, 0, 0, 0, ltd_topo_map, ltd_integrand_body))


    def get_edge_scaling(self, pdg, edge_name=None):

        if edge_name is not None and edge_name in forced_edge_scalings:
            return forced_edge_scalings[edge_name]

        # all scalings that deviate from -2
        scalings = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 11: -1, 12: -1, 13: -1}
        return scalings[abs(pdg)] if abs(pdg) in scalings else -2

    def get_node_scaling(self, pdgs_input, names_of_edges_connected=None):

        if names_of_edges_connected is not None and names_of_edges_connected in forced_node_scalings:
            return forced_node_scalings[names_of_edges_connected]
        
        # Remove dummy particles
        pdgs = tuple([ (-(abs(pdg)%1000) if pdg<0 else pdg%1000) for pdg in pdgs_input if pdg not in [1122,]])

        # only the triple gluon vertex and the ghost gluon vertex have a non-zero scaling
        if pdgs == (25, 21, 21):
            return 2
        elif pdgs == (25, 21, 21, 21) or pdgs == (21, 21, 21) or pdgs == (-82, 21, 82):
            return 1
        else:
            return 0

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, numerator_call, final_state_particle_ids=(),jet_ids=None, write_yaml=True, bar=None,
        integrand_type=None, workspace=None, include_integration_channel_info=True, selected_cuts=None):
        if workspace is None:
            workspace = pjoin(root_output_path, os.pardir, 'workspace')

        if bar:
            bar.update(i_lmb='1')
            max_lmb = 1
            if isinstance(self.additional_lmbs, list):
                max_lmb += len(self.additional_lmbs)
            bar.update(max_lmb='%d'%max_lmb)

        if self.is_zero:
            return False

        # Relabel edges according to alphaLoop conventions:
        for edge_key, edge_data in self.edges.items():
            edge_data['name'] = 'p' + edge_data['name'] if edge_data['type'] == 'virtual' else 'q' + edge_data['name'][1:]

        # TODO: sort such that the first 4 entries are external (it seems to happen by chance now every time)
        edge_map_lin = [(e['name'], e['vertices'][0], e['vertices'][1]) for e in self.edges.values()]
        assert(e[0] != 'q' or int(e[1:]) < 5 for e in edge_map_lin)

        particle_ids = { e['name']: e['PDG'] for e in self.edges.values() }
        particle_masses = { e['name']: model['parameter_dict'][model.get_particle(e['PDG']).get('mass')].real for e in self.edges.values() }

        num_incoming = sum(1 for e in edge_map_lin if e[0][0] == 'q') // 2

        if num_incoming == 1:
            external_momenta = {'q1': [500., 0., 0., 0.], 'q2': [500., 0., 0., 0.]}
            #external_momenta = {'q1': [1., 0., 0., 0.], 'q2': [1., 0., 0., 0.]}
            p = np.array(external_momenta['q1'])
        else:
            external_momenta = { 'q1': [500., 0., 0., 500.], 'q2': [500., 0., 0., -500.], 'q3': [500., 0., 0., 500.], 'q4': [500., 0., 0., -500.] }
            #external_momenta = {'q1': [1., 0., 0., 1.], 'q2': [1., 0., 0., -1.], 'q3': [1., 0., 0., 1.], 'q4': [1., 0., 0., -1.]}
            p = np.array(external_momenta['q1']) + np.array(external_momenta['q2'])

        loop_momenta = []
        n_loops = len(self.edges) - len(self.nodes) + 1
        for loop_var in range(n_loops):
            lm = next((ee['name'], ee['signature'][0][loop_var]) for ee in self.edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] == 1)
            loop_momenta.append(lm)

        if isinstance(self.additional_lmbs, list):
            call_signature_ID = numerator_call
        else:
            call_signature_ID = self.additional_lmbs*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+numerator_call

        # If effective_vertex_id is specified, we will only consider the CC cut that cuts all internal edges connecting to that vertex id
        cut_filter = []
        if self.effective_vertex_id is not None:
            cut_filter.append( tuple( edge_data['name' ] for edge_key, edge_data in self.edges.items() if self.effective_vertex_id in edge_data['vertices'] and edge_data['type']=='virtual' ) )
        
        if isinstance(selected_cuts,dict) and self.name in selected_cuts:
            cut_filter = [ tuple(sorted(sc)) for sc in selected_cuts[self.name] ]

        topo = LTD.squared_topologies.SquaredTopologyGenerator(edge_map_lin,
            self.name, ['q1', 'q2'][:num_incoming], n_jets, external_momenta,
            loop_momenta_names=tuple([l for l,s in loop_momenta]),
            loop_momenta_signs=tuple([s for l,s in loop_momenta]),
            particle_ids=particle_ids,
            masses=particle_masses,
            final_state_particle_ids=final_state_particle_ids,
            jet_ids=jet_ids,
            overall_numerator=1.0,
            numerator_structure={},
            FORM_numerator={'call_signature': {'id': call_signature_ID}},
            FORM_integrand={'call_signature': {'id': call_signature_ID}},
            edge_weights={e['name']: self.get_edge_scaling(e['PDG'],e['name']) for e in self.edges.values()},
            vertex_weights={nv: self.get_node_scaling( n['PDGs'], tuple(sorted([ self.edges[eID]['name'] for eID in n['edge_ids'] ])) ) for nv, n in self.nodes.items()},
            generation_options=FORM_processing_options,
            analytic_result=(self.benchmark_result if hasattr(self,"benchmark_result") else None),
            default_kinematics=self.default_kinematics,
            cut_filter=cut_filter,
            uv_test = FORM_processing_options['uv_test']
        )
        # check if cut is possible
        if len(topo.cuts) == 0:
            logger.info("No cuts for graph {}".format(self.name))
            return False

        self.squared_topology = topo

        # Also generate the squared topology yaml files of all additional LMB topologies used for cross-check
        if isinstance(self.additional_lmbs, list):
            for i_lmb, (_,_,_,other_lmb_supergraph) in enumerate(self.additional_lmbs):
                if bar:
                    bar.update(i_lmb='%d'%(i_lmb+2))
                other_lmb_supergraph.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call, 
                        final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, write_yaml=write_yaml,workspace=workspace,
                        bar=bar, integrand_type=integrand_type, include_integration_channel_info=False, selected_cuts=selected_cuts)

        self.generate_integrand(topo, integrand_type, workspace, call_signature_ID, progress_bar = bar)
        self.generate_form_input(topo)

        if write_yaml:
            if isinstance(self.additional_lmbs, int):
                topo.export(pjoin(root_output_path, "%s_LMB%d.yaml"%(self.name,self.additional_lmbs)),
                    model=model,
                    include_integration_channel_info=include_integration_channel_info,
                    optimize_channels=FORM_processing_options['optimise_integration_channels']
                )
            else:
                topo.export(pjoin(root_output_path, "%s.yaml"%self.name), 
                    model=model,
                    include_integration_channel_info=include_integration_channel_info, 
                    optimize_channels=FORM_processing_options['optimise_integration_channels']
                )

        return True


    def generate_form_input(self, topo):
        configurations = []
        uv_diagrams = []
        uv_forest = []
        topo_map = '#procedure uvmap()\n'

        pure_forest_counter = 0
        for cut_index, cut in enumerate(topo.cuts):
            for diag_set in cut['diagram_sets']:
                diag_set_uv_conf = []
                diag_momenta = []

                has_uv = FORM_processing_options['uv_test'] is None

                n_loops = len(self.edges) - len(self.nodes) + 1
                cmb_offset = len(cut['cuts']) - 1
                for diag_info in diag_set['diagram_info']:
                    # write the entire UV structure as a sum
                    conf = []
                    bubble_uv_derivative = ''
                    for uv_index, uv_structure in enumerate(diag_info['uv']):
                        forest_element = []

                        if FORM_processing_options['uv_test'] is not None:
                            if uv_structure['remaining_graph'].n_loops == 0 and len(uv_structure['uv_subgraphs']) > 0:
                                pure_forest_counter += 1

                                if FORM_processing_options['uv_test'] < 0 or FORM_processing_options['uv_test'] == pure_forest_counter - 1:
                                    conf = []
                                    has_uv = True
                                else:
                                    if uv_index > 0:
                                        continue
                            else:
                                if uv_index > 0:
                                    continue

                        # construct the map from the cmb/lmb to the forest basis
                        forest_to_cb = []
                        for lmb_index, (r, aff, extshift, _, _, _) in enumerate(zip(*uv_structure['forest_to_cb_matrix'])):
                            if all(x == 0 for x in r):
                                assert(all(x == 0 for x in aff))
                                continue
                            mom = ''.join('+{}*fmb{}'.format(a, forest_index + 1) for forest_index, a in enumerate(r) if a != 0)
                            # the shift should be subtracted
                            shift = ''
                            for cmb_index, a in enumerate(aff):
                                if a == 0:
                                    continue

                                # also subtract the external momenta
                                d = self.momenta_decomposition_to_string(([0] * n_loops, cut['cuts'][cmb_index]['signature'][1]), False)
                                if d != '':
                                    shift += '-{}*(c{}-({}))'.format(a, cmb_index + 1, d)
                                else:
                                    shift += '-{}*c{}'.format(a, cmb_index + 1)

                            eshift = self.momenta_decomposition_to_string(([0] * n_loops, extshift), True)
                            if eshift != '':
                                shift += '-({})'.format(eshift)

                            m = 'c{},{}{}'.format(lmb_index + cmb_offset + 1, mom, shift)
                            forest_to_cb.append(m)
                        if len(forest_to_cb) > 0:
                            forest_element.append('cbtofmb({})'.format(','.join(forest_to_cb)))
                        
                        # only needed for spatial part
                        cb_to_forest = []
                        for fmb_index, (_, _, _, r, aff, extshift) in enumerate(zip(*uv_structure['forest_to_cb_matrix'])):
                            if all(x == 0 for x in r):
                                assert(all(x == 0 for x in aff))
                                continue
                            mom = ''.join('+{}*cs{}'.format(a, cmb_index + cmb_offset + 1) for cmb_index, a in enumerate(r) if a != 0)
                            # the shift should be added
                            shift = ''
                            for cmb_index, a in enumerate(aff):
                                if a == 0:
                                    continue

                                d = self.momenta_decomposition_to_string(([0] * n_loops, cut['cuts'][cmb_index]['signature'][1]), False)
                                d = d.replace('p', 'ps')
                                if d != '':
                                    shift += '+{}*(cs{}-({}))'.format(a, cmb_index + 1, d)
                                else:
                                    shift += '+{}*cs{}'.format(a, cmb_index + 1)

                            eshift = self.momenta_decomposition_to_string(([0] * n_loops, extshift), True)
                            eshift = eshift.replace('p', 'ps')
                            if eshift != '':
                                shift += '+({})'.format(eshift)

                            m = 'fmbs{},{}{}'.format(fmb_index + 1, mom, shift)
                            cb_to_forest.append(m)
                        if len(cb_to_forest) > 0:                           
                            forest_element.append('fmbtocb({})'.format(','.join(cb_to_forest)))

                        diag_moms = ','.join(self.momenta_decomposition_to_string(lmm, False) for lmm in uv_structure['remaining_graph_loop_topo'].loop_momentum_map)

                        if len(uv_structure['bubble']) > 0:
                            # we have a bubble: replace the remaining graphs by the sum of bubble derivative graphs
                            bubbles = []
                            for bubble in uv_structure['bubble']:
                                trans = []

                                # if the cutkosky cut has a negative sign, we derive in -p^0 instead of p^0.
                                # here we compensate for this sign
                                trans.append(str(next(c for c in cut['cuts'] if c['edge'] == bubble['derivative'][0])['sign']))

                                ext_mom = next(ee for ee in self.edges.values() if ee['name'] == bubble['derivative'][0])['momentum']
                                der_mom = next(ee for ee in self.edges.values() if ee['name'] == bubble['derivative'][1])['momentum']
                                ext_sig = next(ee for ee in self.edges.values() if ee['name'] == bubble['derivative'][0])['signature']
                                der_sig = next(ee for ee in self.edges.values() if ee['name'] == bubble['derivative'][1])['signature']

                                if bubble['derivative'][0] == bubble['derivative'][1]:
                                    # numerator derivative
                                    bubble_uv_derivative = '{}*der({})*'.format('*'.join(trans), ext_mom)
                                    trans.append('der({})'.format(ext_mom))
                                else:
                                    # check if we pick up a sign change due to the external momentum flowing in the opposite direction
                                    signs = [se * sc for se, sc in zip(der_sig[0] + der_sig[1], ext_sig[0] + ext_sig[1]) if se * sc != 0]
                                    assert(len(set(signs)) == 1)
                                    trans.append('-2*{}der({},0)'.format('1*' if signs[0] == 1 else '-1*', der_mom))

                                if diag_moms != '':
                                    trans.append('diag({},{},{})'.format(diag_set['id'], topo.topo_map[bubble['id']], diag_moms))
                                else:
                                    trans.append('diag({},{})'.format(diag_set['id'], topo.topo_map[bubble['id']]))
                                bubbles.append('*'.join(trans))

                            forest_element.append('({})'.format('+'.join(bubbles)))
                        else:
                            if diag_moms != '':
                                forest_element.append('{}diag({},{},{})'.format(bubble_uv_derivative, diag_set['id'], topo.topo_map[uv_structure['remaining_graph_id']], diag_moms))
                            else:
                                forest_element.append('{}diag({},{})'.format(bubble_uv_derivative, diag_set['id'], topo.topo_map[uv_structure['remaining_graph_id']]))

                        if uv_index == 0:
                            if diag_moms != '':
                                diag_set_uv_conf.append('forestid({},{})'.format(len(uv_forest),diag_moms))
                            else:
                                diag_set_uv_conf.append('forestid({})'.format(len(uv_forest)))

                        for uv_subgraph in uv_structure['uv_subgraphs']:
                            if uv_subgraph['first_occurrence_id'] == uv_subgraph['id']:
                                for dg in uv_subgraph['derived_graphs']:
                                    if dg['skip_pf']:
                                        continue
                                    rp = '*'.join('t{}'.format(i) if raised_power == 1 else 't{}^{}'.format(i, raised_power)
                                            for i, (_,_,raised_power) in enumerate(dg['loop_topo'].uv_loop_lines[0]) if raised_power != 0)
                                    topo_map += '\tid uvtopo({},{},k1?,...,k{}?) = diag({},{},k1,...,k{});\n'.format(uv_subgraph['id'],
                                        '1' if rp == '' else rp, dg['graph'].n_loops, diag_set['id'], topo.topo_map[dg['id']], dg['graph'].n_loops)

                                    if 'soft_ct_id' in dg:
                                        if 'onshell_ct_id' in dg:
                                            topo_map += '\tid irtopo({},{},k1?,...,k{}?) = (1+gamma0)/2*diag({},{},k1,...,k{});\n'.format(uv_subgraph['id'],
                                                '1' if rp == '' else rp, dg['graph'].n_loops, diag_set['id'], topo.topo_map[dg['soft_ct_id']], dg['graph'].n_loops)
                                            topo_map += '\tid irtopo({},{},xneg,k1?,...,k{}?) = (1-gamma0)/2*diag({},{},k1,...,k{});\n'.format(uv_subgraph['id'],
                                                '1' if rp == '' else rp, dg['graph'].n_loops, diag_set['id'], topo.topo_map[dg['onshell_ct_id']], dg['graph'].n_loops)
                                        else:
                                            topo_map += '\tid irtopo({},{},k1?,...,k{}?) = diag({},{},k1,...,k{});\n'.format(uv_subgraph['id'],
                                            '1' if rp == '' else rp, dg['graph'].n_loops, diag_set['id'], topo.topo_map[dg['soft_ct_id']], dg['graph'].n_loops)

                            # construct the vertex structure of the UV subgraph
                            # TODO: are the LTD vertices reliable?
                            uv_loop_graph = uv_subgraph['derived_graphs'][0]['loop_topo'] # all derived graphs have the same topo
                            uv_diag_moms = ','.join(self.momenta_decomposition_to_string(lmm, False) for lmm in uv_loop_graph.loop_momentum_map)
                            vertex_structure = []
                            subgraph_vertices = set(v for ll in uv_loop_graph.loop_lines for v in (ll.start_node, ll.end_node))
                            for v in subgraph_vertices:
                                vertex = []
                                for ll in uv_loop_graph.loop_lines:
                                    for (dv, outgoing) in ((ll.start_node, 1), (ll.end_node, -1)):
                                        if v != dv:
                                            continue
                                        loop_mom_sig = ''
                                        for s, lmm in zip(ll.signature, uv_loop_graph.loop_momentum_map):
                                            if s != 0:
                                                loop_mom_sig += '{}({})'.format('+' if s * outgoing == 1 else '-', self.momenta_decomposition_to_string(lmm, False))
                                        vertex.append(loop_mom_sig)
                                vertex_structure.append('vxs({})'.format(','.join(vertex)))

                            uv_props = []
                            for i, (ll_sig, propagators, _raised_power) in enumerate(uv_loop_graph.uv_loop_lines[0]):
                                loop_mom_sig = ''
                                loop_mom_shift = ''
                                for s, lmm, lm_shift in zip(ll_sig, uv_loop_graph.loop_momentum_map, uv_loop_graph.uv_loop_lines[1]):
                                    if s != 0:
                                        loop_mom_sig += '{}({})'.format('+' if s == 1 else '-', self.momenta_decomposition_to_string(lmm, False))
                                        ds = self.momenta_decomposition_to_string(lm_shift, False)
                                        if ds != '':
                                            loop_mom_shift += '{}({})'.format('+' if s == 1 else '-', ds)

                                # the parametric shift is given in terms of external momenta of the subgraph
                                # translate the signature and param_shift to momenta of the supergraph
                                for (edge_name, param_shift, power) in propagators:
                                    ext_mom_sig = ''
                                    edge_mass = 'masses({})'.format(next(ee for ee in self.edges.values() if ee['name'] == edge_name)['PDG'])

                                    if all(s == 0 for s in param_shift[1]):
                                        for _ in range(power):
                                            if loop_mom_shift == '':
                                                uv_props.append('uvprop({},t{},0,{})'.format(loop_mom_sig, i, edge_mass))
                                            else:
                                                uv_props.append('uvprop({},t{},{},{})'.format(loop_mom_sig, i, loop_mom_shift, edge_mass))
                                        continue

                                    for (ext_index, s) in enumerate(param_shift[1]):
                                        if s != 0:
                                            ext_mom = uv_subgraph['graph'].edge_map_lin[uv_subgraph['graph'].ext[ext_index]][0]
                                            ext_edge = next(ee for ee in self.edges.values() if ee['name'] == ext_mom)
                                            ext_mom_sig += '{}({})'.format('+' if s == 1 else '-',
                                                self.momenta_decomposition_to_string(ext_edge['signature'], False))

                                    for _ in range(power):
                                        uv_props.append('uvprop({},t{},{},{})'.format(loop_mom_sig, i, loop_mom_shift + ext_mom_sig, edge_mass))
                            # it could be that there are no propagators with external momentum dependence when pinching duplicate edges
                            if uv_props == []:
                                uv_props = ['1']

                            if uv_diag_moms == '':
                                # should never happen!
                                logger.warn("No diag moms in UV graph")
                                uv_diag = 'uvtopo({})'.format(uv_subgraph['first_occurrence_id'])
                            else:
                                uv_diag = 'uvtopo({},{})'.format(uv_subgraph['first_occurrence_id'], uv_diag_moms)

                            if uv_subgraph['onshell'] is not None and FORM_processing_options['on_shell_renormalisation']:
                                ext_edge_sig = next(ee for ee in self.edges.values() if ee['name'] == uv_subgraph['onshell'][0])
                                totalmom = self.momenta_decomposition_to_string(ext_edge_sig['signature'], False)
                                uv_diag += '*onshell({},masses({}))'.format(totalmom, ext_edge_sig['PDG'])

                            if FORM_processing_options['generate_integrated_UV_CTs']:
                                uv_diag += '*intuv(1 - {}*diag({},{},{}))'.format('*'.join(vertex_structure), diag_set['id'], topo.topo_map[uv_subgraph['first_occurrence_id'] + 1], uv_diag_moms)

                            uv_conf_diag = '-tmax^{}*{}*{}'.format(uv_subgraph['taylor_order'],'*'.join(uv_props),uv_diag)
                            if uv_conf_diag not in uv_diagrams:
                                uv_diagrams.append(uv_conf_diag)

                            uv_conf = 'uvdiag({})'.format(uv_diagrams.index(uv_conf_diag))

                            sg_call = 'subgraph({}{},{},{})'.format(uv_subgraph['graph_index'],
                                (',' if len(uv_subgraph['subgraph_indices']) > 0 else '') + ','.join(str(si) for si in uv_subgraph['subgraph_indices']),
                                uv_conf, uv_diag_moms)

                            forest_element.append(sg_call)

                        conf.append('*'.join(forest_element))

                    cmb_offset += len(diag_info['uv'][0]['forest_to_cb_matrix'][0][0]) # add the amplitude loop count to the cmb start
                    uv_forest.append('+\n\t'.join(conf))

                # construct the map from the lmb to the cmb
                cmb_map = []
                for i in range(n_loops):
                    s = ''
                    for cc, cs in enumerate(diag_set['cb_to_lmb'][i * n_loops:i * n_loops+n_loops]):
                        if cs != 0:
                            s += '{}c{}'.format('+' if cs == 1 else '-', cc + 1)
                            if cc < len(cut['cuts'][:-1]):
                                d = self.momenta_decomposition_to_string(([0] * n_loops, cut['cuts'][cc]['signature'][1]), False)
                                if d != '':
                                    # note the sign inversion
                                    s += '{}({})'.format('-' if cs == 1 else '+', d)
                    cmb_map.append(('k' + str(i + 1), s))

                cmb_map = 'cmb({})'.format(','.join(d for c in cmb_map for d in c))

                # store which momenta are LTD momenta
                conf = 'c0' if n_loops == len(cut['cuts'][:-1]) else ','.join('c' + str(i + 1) for i in range(len(cut['cuts'][:-1]), n_loops) )

                if diag_momenta == []:
                    diag_momenta = ['1']

                tder = '*tder({})'.format(','.join(str(c['power']) for c in cut['cuts'] if c['power'] > 1)) if any(c['power'] > 1 for c in cut['cuts']) else ''
                conf = 'conf({},{},{},{}){}*{}'.format(diag_set['id'], cut_index, cmb_map, conf, tder, '*'.join(diag_momenta))
                if diag_set_uv_conf != []:
                    conf += '*{}'.format('*'.join(diag_set_uv_conf))

                if not has_uv:
                    continue

                configurations.append(conf)

            if len(configurations) > 0:
                configurations[-1] += '\n'

        topo_map += '#endprocedure\n'

        self.configurations = (topo_map, uv_diagrams, uv_forest, ' +'.join(configurations))


class FORMSuperGraphIsomorphicList(list):
    """ Container class for a list of FROMSuperGraph with the same denominator structure"""
    
    def __init__(self, graph_list):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraph instances or LTD2SuperGraph instances."""

        for g in graph_list:
            if isinstance(g,FORMSuperGraph):
                self.append(g)
            else:
                self.append(FORMSuperGraph.from_LTD2SuperGraph(g))

    def generate_numerator_form_input(self, additional_overall_factor='',):
        return '+'.join(g.generate_numerator_form_input(additional_overall_factor) for g in self)

    def analyze_FORM_output(self, FORM_stdout, FORM_vars):
        """ Analyze the FORM otuput to gather statistics about the generation."""

        code_generation_statistics = {}

        original_op_count_re = re.compile("^\*\*\* STATS: original\s*(?P<npower>\d+)P\s*(?P<nmult>\d+)M\s*(?P<nadd>\d+)A.*$")
        optimized_op_count_re = re.compile("^\*\*\* STATS: optimized\s*(?P<npower>\d+)P\s*(?P<nmult>\d+)M\s*(?P<nadd>\d+)A.*$")

        accumulation_mode = None
        current_optimized = None
        for line in FORM_stdout.split('\n'):
            if 'START integrand' in line:
                accumulation_mode = 'integrand_%s'%FORM_vars['INTEGRAND']
                continue
            if 'END integrand' in line:
                # register last displayed optimized count
                if current_optimized is not None:
                    if '%s_optimized_op_count'%accumulation_mode not in code_generation_statistics:
                        code_generation_statistics['%s_optimized_op_count'%accumulation_mode] = current_optimized
                    else:
                        for k, v in current_optimized.items():
                            code_generation_statistics['%s_optimized_op_count'%accumulation_mode][k] += v
                # Reset curent_optimized
                current_optimized = None
                accumulation_mode = None
                continue
            if 'START numerator' in line:
                accumulation_mode = 'numerator'
                continue
            if 'END numerator' in line:
                # register last displayed optimized count
                if current_optimized is not None:
                    if '%s_optimized_op_count'%accumulation_mode not in code_generation_statistics:
                        code_generation_statistics['%s_optimized_op_count'%accumulation_mode] = current_optimized
                    else:
                        for k, v in current_optimized.items():
                            code_generation_statistics['%s_optimized_op_count'%accumulation_mode][k] += v
                # Reset curent_optimized
                current_optimized = None
                accumulation_mode = None
                continue
            if accumulation_mode is None:
                continue

            original_op_count_match = original_op_count_re.match(line)
            if original_op_count_match is not None:
                # register last displayed optimized count
                if current_optimized is not None:
                    if '%s_optimized_op_count'%accumulation_mode not in code_generation_statistics:
                        code_generation_statistics['%s_optimized_op_count'%accumulation_mode] = current_optimized
                    else:
                        for k, v in current_optimized.items():
                            code_generation_statistics['%s_optimized_op_count'%accumulation_mode][k] += v
                # Then register new original op count
                orig_op_count = {
                    'additions': int(original_op_count_match.group('nadd')),
                    'multiplications': int(original_op_count_match.group('nmult')),
                    'powers': int(original_op_count_match.group('npower'))
                }
                if '%s_original_op_count'%accumulation_mode not in code_generation_statistics:
                    code_generation_statistics['%s_original_op_count'%accumulation_mode] = orig_op_count
                else:
                    for k, v in orig_op_count.items():
                        code_generation_statistics['%s_original_op_count'%accumulation_mode][k] += v
                # Reset curent_optimized
                current_optimized = None
                continue

            optimized_op_count_match = optimized_op_count_re.match(line)            
            if optimized_op_count_match is not None:
                # Overwrite the new value from the current optimized nop
                current_optimized = {
                    'additions': int(optimized_op_count_match.group('nadd')),
                    'multiplications': int(optimized_op_count_match.group('nmult')),
                    'powers': int(optimized_op_count_match.group('npower'))
                }

        #print(FORM_stdout)
        #print(code_generation_statistics)

        # Add the compression level metric if the necessary info is present.
        for info in ['numerator','integand_LTD','integand_PF']:
            if ('%s_original_op_count'%info in code_generation_statistics) and ('%s_optimized_op_count'%info in code_generation_statistics):
                for op_type in ['additions','multiplications','powers']:
                    compression = 0.
                    if code_generation_statistics['%s_original_op_count'%info][op_type]!=0:
                        compression = 100.0*(1.-(float(code_generation_statistics['%s_optimized_op_count'%info][op_type])/float(code_generation_statistics['%s_original_op_count'%info][op_type])))
                    if '%s_compression_percentage'%info not in code_generation_statistics:
                        code_generation_statistics['%s_compression_percentage'%info] = {}
                    code_generation_statistics['%s_compression_percentage'%info][op_type]=float('%.2f'%compression)

        return code_generation_statistics

    def generate_numerator_functions(self, model, additional_overall_factor='', output_format='c', workspace=None, FORM_vars=None, active_graph=None,process_definition=None, recycle=False, forced_options=None):
        """ Use form to plugin Feynman Rules and process the numerator algebra so as
        to generate a low-level routine in file_path that encodes the numerator of this supergraph."""

        FORM_source_to_run = 'numerator'
        if output_format in ['raw','pySecDec']:
            FORM_source_to_run = 'raw_numerator'

        _MANDATORY_FORM_VARIABLES = ['SGID','NINITIALMOMENTA','NFINALMOMENTA','SELECTEDEPSILONORDER','UVRENORMFINITEPOWERTODISCARD','OPTIMISATIONSTRATEGY']

        if forced_options is not None:
            FORM_processing_options = forced_options

        if FORM_vars is None:
            raise FormProcessingError("FORM_vars must be supplied when calling generate_numerator_functions.")
        FORM_vars = dict(FORM_vars)

        if active_graph is None:
            characteristic_super_graph = self[0]
        else:
            characteristic_super_graph = active_graph

        if 'NINITIALMOMENTA' not in FORM_vars:
            n_incoming = sum([1 for edge in characteristic_super_graph.edges.values() if edge['type'] == 'in'])
            FORM_vars['NINITIALMOMENTA'] = n_incoming
        if 'NFINALMOMENTA' not in FORM_vars:
            n_loops = len(characteristic_super_graph.edges) - len(characteristic_super_graph.nodes) + 1
            FORM_vars['NFINALMOMENTA'] = n_loops

        # set the order
        if FORM_processing_options['generate_integrated_UV_CTs']:
            FORM_vars['MAXPOLE'] = len(characteristic_super_graph.edges) - len(characteristic_super_graph.nodes) + 1 - min(len(c['cuts']) - 1 for c in characteristic_super_graph.squared_topology.cuts)
        else:
            FORM_vars['MAXPOLE'] = 0

        if FORM_vars is None or not all(opt in FORM_vars for opt in _MANDATORY_FORM_VARIABLES):
            raise FormProcessingError("The following variables must be supplied to FORM: %s"%str(_MANDATORY_FORM_VARIABLES))

        FORM_vars.update(FORM_processing_options['extra-options'])
        
        i_graph = int(FORM_vars['SGID'])

        # write the form input to a file
        if active_graph is None:
            form_input = self.generate_numerator_form_input(additional_overall_factor)
        else:
            form_input = characteristic_super_graph.generate_numerator_form_input(additional_overall_factor)

        if workspace is None:
            selected_workspace = FORM_workspace
            FORM_source = pjoin(plugin_path,"%s.frm"%FORM_source_to_run)
        else:
            selected_workspace = workspace
            FORM_source = pjoin(selected_workspace,'%s.frm'%FORM_source_to_run)

        with open(pjoin(selected_workspace,'input_%d.h'%i_graph), 'w') as f:
            (uv_map, uv_conf, uv_forest, conf) = characteristic_super_graph.configurations
            f.write(uv_map + '\n')
            f.write('CTable uvdiag(0:{});\n'.format(len(uv_conf)))
            f.write('{}\n\n'.format('\n'.join('Fill uvdiag({}) = {};'.format(i, uv) for i,uv in enumerate(uv_conf))))
            f.write('CTable forest(0:{});\n'.format(len(uv_forest)))
            f.write('{}\n\n'.format('\n'.join('Fill forest({}) = {};'.format(i, uv) for i,uv in enumerate(uv_forest))))
            f.write('L CONF =\n +{};\n\n'.format(conf))
            f.write('L F = {}\n;'.format(form_input))
        
        if platform.system() != 'Darwin':
            form_settings = formset.generate_form_settings(percentage=FORM_processing_options["FORM_use_max_mem_fraction"] * 100., ncpus=
                FORM_processing_options["FORM_parallel_cores"] * FORM_processing_options["cores"])
        else:
            form_settings = {}
        form_settings.update(FORM_processing_options["FORM_setup"])

        with open(pjoin(selected_workspace,'form.set'), 'w') as f:
            content = [ '%s %s'%(k,str(v)) for k,v in form_settings.items() ]
            f.write('\n'.join(content))

        if FORM_processing_options["FORM_parallel_cores"] == 1:
            FORM_cmd = [FORM_processing_options["FORM_path"]]
        else:
            FORM_cmd = [FORM_processing_options["tFORM_path"], '-w{}'.format(FORM_processing_options["FORM_parallel_cores"]), '-W']

        FORM_cmd = ' '.join(FORM_cmd +
                [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ] +
                [ '-M', '-l', '-C', '%s_%s.log'%(FORM_source_to_run,i_graph)] +
                [ FORM_source, ]
        )
        if recycle and os.path.isfile(pjoin(selected_workspace,'FORM_run_cmd_%d.exe'%i_graph)) and \
            os.path.isfile(pjoin(selected_workspace,'numerator_%d.log'%i_graph)) and \
            ( (FORM_vars['INTEGRAND'] not in ['PF','both']) or os.path.isfile(pjoin(selected_workspace,'out_integrand_PF_%d.proto_c'%i_graph)) ) and \
            ( (FORM_vars['INTEGRAND'] not in ['LTD','both']) or os.path.isfile(pjoin(selected_workspace,'out_integrand_LTD_%d.proto_c'%i_graph)) ):
            logger.info("Recycling already existing FORM output for graph #%d. Remove file out_integrand_*_%d.proto_c to force a reprocessing."%(i_graph, i_graph))
            characteristic_super_graph.code_generation_statistics = self.analyze_FORM_output(open(pjoin(selected_workspace,'numerator_%d.log'%i_graph),'r').read(), FORM_vars)
            return True

        with open(pjoin(selected_workspace,'FORM_run_cmd_%d.exe'%i_graph), 'w') as f:
            f.write(FORM_cmd)

        r = subprocess.run(FORM_cmd,
            shell=True,
            cwd=selected_workspace,
            capture_output=True)

        # TODO understand why FORM sometimes returns return_code 1 event though it apparently ran through fine
        characteristic_super_graph.code_generation_statistics = self.analyze_FORM_output(r.stdout.decode("utf-8"), FORM_vars)


        if r.returncode != 0:# or not os.path.isfile(pjoin(selected_workspace,'out_integrand_%d.proto_c'%i_graph)):
            raise FormProcessingError("FORM processing failed with error:\n%s\nFORM command to reproduce:\ncd %s; %s"%(
                    r.stdout.decode('UTF-8'),
                    selected_workspace, FORM_cmd
            ))

        return True

        # TODO: numerators are deprecated
        # return the code for the numerators
        #if not os.path.isfile(pjoin(selected_workspace,'out_%d.proto_c'%i_graph)):
        #    raise FormProcessingError(
        #            ( "FORM failed to produce an output for super graph ID=%d. Output file not found at '%s'."%
        #                                                            (i_graph,pjoin(selected_workspace,'out_%d.proto_c'%i_graph)))+
        #        "\nFORM command to reproduce:\ncd %s; %s"%(selected_workspace,FORM_cmd)
        #    )

        #with open(pjoin(selected_workspace,'out_%d.proto_c'%i_graph), 'r') as f:
        #    num_code = f.read()
#
        ## TODO Remove when FORM will have fixed its C output bug
        #num_code = temporary_fix_FORM_output(num_code)
#
        #return num_code

    def to_dict(self, file_path=None):
        """ Store that into a dict."""
        to_dump = [g.to_dict() for g in self]
        if file_path:
            open(file_path,'w').write(pformat(to_dump))
        else:
            return to_dump
    
    def multiplicity_factor(self, iso_id, workspace, form_source):
        
        #VH HACK bernard
        #return (self, 1)

        output_match = re.compile(r'isoF=(.*?);')
        reference = self[0].generate_numerator_form_input('', only_algebra=True)
        FORM_vars = {}
        FORM_vars['SGID'] = iso_id
        FORM_vars['NUMD'] = (len(self) if len(self)!=1 else 2)
        FORM_vars['FOURDIM'] = 1
        FORM_vars['MAXPOLE'] = 0 # no poles will be generated
        FORM_vars['SELECTEDEPSILONORDER'] = 0

        with open(pjoin(workspace,'iso_check_{}.frm'.format(iso_id)), 'w') as f:
            for i_graph, g in enumerate(self):
                mapped = g.generate_numerator_form_input('', only_algebra=True)
                f.write("L F{} = {};\n".format(i_graph + 1, mapped))
                if len(self) == 1:
                    f.write("L F2 = {};\n".format(mapped))
            
        cmd = ' '.join([
            FORM_processing_options["FORM_path"],
            ]+
            [ '-D %s=%s'%(k,v) for k,v in FORM_vars.items() ]+
            [ form_source, ]
        )

        with open(pjoin(workspace,'FORM_run_cmd_iso_check_%d.exe'%iso_id), 'w') as f:
            f.write(cmd)

        r = subprocess.run(cmd, shell=True, cwd=workspace, capture_output=True)
        if r.returncode != 0:
            raise FormProcessingError("FORM processing failed with error:\n%s"%(r.stdout.decode('UTF-8')))

        output = r.stdout.decode('UTF-8').replace(' ','').replace('\n','')
        if output_match.findall(output)[0] == "0":
            return(self, 0)
        if len(self) == 1:
            return(self, 1)

        factor = re.sub(r'rat\(([-0-9]+),([-0-9]+)\)', r'(\1)/(\2)', output_match.findall(output)[0])
        if "rat" in factor:
            raise FormProcessingError("Multiplicity not found: {} / {} is not rational (iso_check_%(SGID)d)".format(self[0].name,g.name )%FORM_vars)

        multiplicity = 1 + eval(factor)

        #logger.info("{} = {} * {} (id {})".format(self[0].name, factor, g.name, iso_id))
        return (self, multiplicity)

    @staticmethod
    def multiplicity_factor_helper(args):
        return args[0].multiplicity_factor(*args[1:])

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, numerator_call, final_state_particle_ids=(), jet_ids=None, workspace=None, bar=None,
            integrand_type=None,include_integration_channel_info=True, selected_cuts=None):
        for i, g in enumerate(self):
            # Now we generate the squared topology only for the first isomorphic graph
            # to obtain the replacement rules for the bubble.
            # The other elements of the isomorphic set are going to contribute only at the 
            # numerator 
            if i==0:
                r = g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call, 
                                final_state_particle_ids, jet_ids=jet_ids, write_yaml=i==0, workspace=workspace, bar=bar, integrand_type=integrand_type, 
                                include_integration_channel_info=include_integration_channel_info, selected_cuts=selected_cuts)
        #print(r)
        return r

    def generate_numerator_file(self, i_graph, model, root_output_path, additional_overall_factor, workspace, integrand_type,  process_definition, header_map, output_format, recycle, forced_options=None):
        timing = time.time()
        self.is_zero = True

        # add all numerators in one file and write the headers
        if forced_options is not None:
            FORM_processing_options=forced_options

        FORM_vars={
            'SELECTEDEPSILONORDER':'%d'%FORM_processing_options['selected_epsilon_UV_order'],
            'OPTIMISATIONSTRATEGY':FORM_processing_options['optimisation_strategy'],
        }
        if not FORM_processing_options['perform_msbar_subtraction']:
            FORM_vars['NOMSBARSUBTRACTION'] = 1
        if FORM_processing_options['uv_test']:
            FORM_vars['UVTEST'] = 1

        if FORM_processing_options['renormalisation_finite_terms']=='together':
            # Keep all terms, so set the discared power to 0
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = 0
        elif FORM_processing_options['renormalisation_finite_terms']=='only':
            # Discard all terms not prefixed with UVRenormFinite, so set discarded power to -1
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = -1
        elif FORM_processing_options['renormalisation_finite_terms']=='removed':
            # Discard all terms prefixed with UVRenormFinite, so set discarded power to 1
            FORM_vars['UVRENORMFINITEPOWERTODISCARD'] = 1
        else:
            raise FormProcessingError("The FORM processing option 'renormalisation_finite_terms' "+
                                      "can only take the following value: 'together', 'only' or 'removed', but not '%s'."%FORM_processing_options['renormalisation_finite_terms'])


        if integrand_type is not None:
            FORM_vars['INTEGRAND'] = integrand_type

        graphs_to_process = []
        if isinstance(self[0].additional_lmbs,list) and self[0].additional_lmbs != []:
            graphs_to_process.append( (0,i_graph, self[0]) )
            graphs_to_process.extend([(g.additional_lmbs,i_graph,g) for _,_,_,g in self[0].additional_lmbs])
        else:
            # By setting the active graph to None we will then sum overall members of the iso set.
            graphs_to_process.append( (0,i_graph, None) )

        max_buffer_size = 0
        all_ids = []
        for i_lmb, i_g, active_graph in graphs_to_process:
            i = i_lmb*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+i_g
            all_ids.append(i)
            FORM_vars['SGID']='%d'%i

            num = self.generate_numerator_functions(model, additional_overall_factor,
                workspace=workspace, FORM_vars=FORM_vars, active_graph=active_graph,
                process_definition=process_definition, output_format=output_format, recycle=recycle, forced_options=forced_options)

        return (i_graph, all_ids, max_buffer_size, not num, time.time() - timing, self[0].code_generation_statistics)

    def generate_numerator_file_helper(args):
        return args[0].generate_numerator_file(*args[1:-1],forced_options=args[-1])

class FORMSuperGraphList(list):
    """ Container class for a list of FORMSuperGraphIsomorphicList."""

    extension_names = {'c': 'c','raw': 'raw', 'pySecDec' : 'pySecDec'}

    def __init__(self, graph_list, name='all_QG_supergraphs'):
        """ Instantiates a list of FORMSuperGraphs from a list of either
        FORMSuperGraphIsomorphicList instances, FORMSuperGraph, or LTD2SuperGraph instances."""

        self.name = name

        for g in graph_list:
            if isinstance(g, FORMSuperGraph):
                self.append(FORMSuperGraphIsomorphicList([g]))
            elif isinstance(g, FORMSuperGraphIsomorphicList):
                self.append(g)
            else:
                self.append(FORMSuperGraphIsomorphicList([FORMSuperGraph.from_LTD2SuperGraph(g)]))

        self.code_generation_statistics = {}

    def aggregate_code_generation_statistics(self):
        """ Combine the FORM code generation statistics from all SG in this list."""
        
        for sg in self:
            for iso_sg in sg:
                for k, v in iso_sg.code_generation_statistics.items():
                    if 'compression_percentage' in k:
                        continue
                    if k in self.code_generation_statistics:
                        if not isinstance(v, dict):
                            self.code_generation_statistics[k] += v
                        else:
                            for kk, vv in v.items():
                                self.code_generation_statistics[k][kk] += vv
                    else:
                        self.code_generation_statistics[k] = v

        # Add the compression level metric if the necessary info is present.
        for info in ['numerator','integrand_LTD','integrand_PF']:
            if ('%s_original_op_count'%info in self.code_generation_statistics) and ('%s_optimized_op_count'%info in self.code_generation_statistics):
                for op_type in ['additions','multiplications','powers']:
                    compression = 0.
                    if self.code_generation_statistics['%s_original_op_count'%info][op_type]!=0:
                        compression = 100.0*(1.-(float(self.code_generation_statistics['%s_optimized_op_count'%info][op_type])/float(self.code_generation_statistics['%s_original_op_count'%info][op_type])))
                    if '%s_compression_percentage'%info not in self.code_generation_statistics:
                        self.code_generation_statistics['%s_compression_percentage'%info] = {}
                    self.code_generation_statistics['%s_compression_percentage'%info][op_type]=float('%.2f'%compression)

    @classmethod
    def from_squared_topology(cls, edge_map_lin, name, incoming_momentum_names, model, loop_momenta_names=None, particle_ids={},benchmark_result=0.0, 
                                default_kinematics=None, effective_vertex_id=None, overall_factor='1'):
        vertices = [v for e in edge_map_lin for v in e[1:]]

        topo_generator = LTD.ltd_utils.TopologyGenerator(edge_map_lin, {})
        topo_generator.generate_momentum_flow(loop_momenta_names)
        sig = topo_generator.get_signature_map()

        edges = {i: {
            'PDG': 3370 if e[0] not in particle_ids else particle_ids[e[0]],
            'indices': (1 + i * 2,) if vertices.count(e[1]) == 1 or vertices.count(e[2]) == 1 else (1 + i * 2, 1 + i * 2 + 1),
            'signature': [sig[e[0]][0],
                    [ i+o for i,o in zip(sig[e[0]] [1][:len(sig[e[0]] [1])//2], sig[e[0]] [1][len(sig[e[0]] [1])//2:])]],
            'name': e[0],
            'type': 'in' if e[0] in incoming_momentum_names else ('out' if vertices.count(e[1]) == 1 or vertices.count(e[2]) == 1 else 'virtual'),
            'vertices': tuple(e[1:]),
            } for i, e in enumerate(edge_map_lin)
        }

        nodes = {v: {
            'PDGs': tuple(e['PDG'] if (v == e['vertices'][1] or e['PDG'] in dummy_scalar_PDGs) else model.get_particle(e['PDG']).get_anti_pdg_code() for e in edges.values() if v in e['vertices']),
            'edge_ids': tuple(ei for ei, e in edges.items() if v in e['vertices']),
            'indices': tuple(e['indices'][0] if v == e['indices'][0] or len(e['indices']) == 1 else e['indices'][1] for e in edges.values() if v in e['vertices']),
            'momenta': tuple('DUMMY' for e in edges.values() if v in e['vertices']),
            'vertex_id': -1 if vertices.count(v) == 1 else 100,
            } for v in set(vertices)
        }

        for eid, e in edges.items():
            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)
            neg_mom =  FORMSuperGraph.momenta_decomposition_to_string([[-s for s in sp] for sp in e['signature']], set_outgoing_equal_to_incoming=False)

            # TODO: why was this here?
            #if len(e['vertices']) == 1:
            #    continue

            for i, vi in enumerate(e['vertices']):
                e_index = nodes[vi]['edge_ids'].index(eid)
                mom = list(nodes[vi]['momenta'])
                mom[e_index] = neg_mom if i == 0 else e['momentum']
                nodes[vi]['momenta'] = mom

        # set the correct particle ordering for all the edges
        for n in nodes.values():
            if len(n['PDGs']) > 1 and all(nn not in dummy_scalar_PDGs for nn in n['PDGs']):
                edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                    n[g] = tuple(n[g][eo['index']] for eo in edge_order)

        form_graph = FORMSuperGraph(
            name=name, edges=edges, nodes=nodes,
            overall_factor=overall_factor,
            multiplicity = 1,
            benchmark_result=benchmark_result,
            default_kinematics=default_kinematics,
            effective_vertex_id=effective_vertex_id
        )

        return FORMSuperGraphList([FORMSuperGraphIsomorphicList([form_graph])], name=name + '_set')


    @classmethod
    def from_dict(cls, dict_file_path, first=None, merge_isomorphic_graphs=False, verbose=False, model = None, workspace=None, cuts=None):
        """ Creates a FORMSuperGraph list from a dict file path."""
        from pathlib import Path
        p = Path(dict_file_path)
        sys.path.insert(0, str(p.parent))

        # compile the file first before importing it
        # with the same optimization flag as MG5
        # this avoids that memory isn't freed after compiling
        # when using the __import__ directly
        logger.info("Compiling imported supergraphs.")
        subprocess.run([sys.executable, '-O', '-m', p.stem], cwd=p.parent)
        m = __import__(p.stem)
        logger.info("Imported {} supergraphs.".format(len(m.graphs)))

        # Filter specific graphs by name 
        # filter_graphs = ['SG_QG3','SG_QG4']
        # m.graphs = [ g for (g,name) in zip(m.graphs, m.graph_names) if name in filter_graphs]
        # m.graph_names = [name for name in m.graph_names if name in filter_graphs ]

        # Now convert the vertex names to be integers according to QGRAF format:
        for i, g in enumerate(m.graphs):
            if isinstance(list(g['nodes'].keys())[0], str):
                new_nodes ={}
                node_names={}
                for i_node, (n_key, n) in enumerate(g['nodes'].items()):
                    node_names[n_key]=i_node+1
                    new_nodes[i_node+1] = n
                for n in g['nodes'].values():
                    n['edge_ids'] = tuple([ (node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2]) for edge_key in n['edge_ids'] ])
                new_edges = {}
                for edge_key, edge in g['edges'].items():
                    edge['vertices']=(node_names[edge_key[0]],node_names[edge_key[1]])
                    new_edges[(node_names[edge_key[0]],node_names[edge_key[1]],edge_key[2])]=edge
                g['nodes'] = new_nodes
                g['edges'] = new_edges
            else:
                # if the diagram comes from QGRAF then we want to rectify the flow of the loop momenta such 
                # that there is consistency within the same loop line
                topo_generator = LTD.ltd_utils.TopologyGenerator([(e['name'], e['vertices'][0], e['vertices'][1]) for e in g['edges'].values()])
                topo_generator.generate_momentum_flow()
                smap = topo_generator.get_signature_map()
                for e_key, e in g['edges'].items():
                    #print("\t",e['momentum'], " -> ", FORMSuperGraph.momenta_decomposition_to_string(smap[e['name']], set_outgoing_equal_to_incoming=True))
                    e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(smap[e['name']], set_outgoing_equal_to_incoming=True)
                    #print(m.graphs[i]['edges'][e_key])
                for n_key, n in g['nodes'].items():
                    if n['vertex_id'] < 0:
                        continue
                    new_moms = []
                    for idx, eid in zip(n['indices'], n['edge_ids']):
                        if g['edges'][eid]['type'] == 'in':
                            sgn = 1
                        elif g['edges'][eid]['type'] == 'out':
                            sgn = -1
                        else:
                            sgn = 2*g['edges'][eid]['indices'].index(idx)-1
                    
                        if not sgn in [1, -1]:
                            raise ValueError
                        
                        signature = [sgn * x for x in smap[g['edges'][eid]['name']][0]]
                        shifts = [sgn * x for x in smap[g['edges'][eid]['name']][1]]
                        new_moms += [FORMSuperGraph.momenta_decomposition_to_string([signature, shifts], set_outgoing_equal_to_incoming=True)]
                    #print("\t\t",n['momenta'])
                    #print("\t\t",tuple(new_moms))
                    n['momenta'] = tuple(new_moms)
                    #print("\t\t",m.graphs[i]['nodes'][n_key]['momenta'])

        full_graph_list = []
        for i, g in enumerate(m.graphs):
            if hasattr(m,'graph_names'):
                graph_name=m.graph_names[i]
            else:
                graph_name=p.stem + '_' + str(i)
            # convert to FORM supergraph
            form_graph = FORMSuperGraph(name=graph_name, edges = g['edges'], nodes=g['nodes'], 
                        overall_factor=g['overall_factor'], multiplicity = g.get('multiplicity',1) )
            form_graph.derive_signatures()
            full_graph_list.append(form_graph)

        # Adjust edge orientation of all repeated edges so that they have identical signatures
        for g in full_graph_list:

            # First detect edges that need flipping
            all_signatures = {}
            for e_key, e in g.edges.items():
                pos_signature = tuple(tuple(s for s in sp) for sp in e['signature'])
                neg_signature = tuple(tuple(-s for s in sp) for sp in e['signature'])
                if pos_signature in all_signatures:
                    all_signatures[pos_signature].append((e_key,+1))
                elif neg_signature in all_signatures:
                    all_signatures[neg_signature].append((e_key,-1))
                else:
                    if pos_signature[0].count(0) == len(pos_signature[0]) - 1 and\
                        pos_signature[1].count(0) == len(pos_signature[1]):
                        if next(s for s in pos_signature[0] if s != 0) < 0:
                            all_signatures[neg_signature] = [(e_key, -1)]
                            break
                        else:
                            all_signatures[pos_signature] = [(e_key,+1)]
                    else:
                        all_signatures[pos_signature] = [(e_key,+1)]

            for sig, edges in all_signatures.items():
                for edge_key, sign in edges:
                    if sign==+1:
                        continue
                    # For now we only allow to flip particles that are not fermionic otherwise 
                    # special care must be taken when flipping an edge and normally fermionic lines
                    # should already have been aligned at this stage.
                    flipped_part = model.get_particle(g.edges[edge_key]['PDG'])
                    if flipped_part.is_fermion():
                        raise FormProcessingError("A *fermionic* repeated edge with opposite signatures was found. This should happen for bosons only.")

                    g.edges[edge_key]['PDG'] = flipped_part.get_anti_pdg_code()
                    g.edges[edge_key]['signature'] = [[-s for s in sp] for sp in g.edges[edge_key]['signature']]
                    g.edges[edge_key]['indices'] = tuple([
                        g.edges[edge_key]['indices'][1],
                        g.edges[edge_key]['indices'][0],
                    ])
                    g.edges[edge_key]['vertices'] = tuple([
                        g.edges[edge_key]['vertices'][1],
                        g.edges[edge_key]['vertices'][0],
                    ])

            # Now adust the string momenta of edges and nodes accordingly.            
            g.impose_signatures()

        if first is not None:
            logger.info("Taking first {} supergraphs.".format(first))
            full_graph_list = full_graph_list[:first]

        import sympy as sp
        # group all isomorphic graphs
        iso_groups = []
        n_externals = max(len([1 for e in graph.edges.values() if e['type']=='in' or e['type']=='out']) for graph in full_graph_list)
        if model is None:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([1,2,3,4,5,6,11,12,13,21,22,25,82])}
        else:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([p['pdg_code'] for p in model['particles']])}
        
        with progressbar.ProgressBar(prefix='Group ISO graphs: {variables.graph_nr} -> {variables.iso_size} : ', max_value=len(full_graph_list),variables={'graph_nr':'0', 'iso_size':'0'}) as bar:
            for g_id, graph in enumerate(full_graph_list):
                g = igraph.Graph()
                g.add_vertices(len(graph.nodes))
                undirected_edges = set()

                for e in graph.edges.values():
                    undirected_edges.add(tuple(sorted(e['vertices'])))

                edge_colors = []
                ext_id = 0
                for ue in undirected_edges:
                    e_color = 1
                    for e in graph.edges.values():
                        if tuple(sorted(e['vertices'])) == ue:
                            #TODO: Reserve the first #externals primes for external edges
                            # with the current implementation it still allows swap 
                            # of final/initial state edges
                            if e['type'] == 'in':
                                ext_id += 1
                                e_color *= sp.prime(ext_id)
                            else:
                                e_color *= pdg_primes[abs(e['PDG'])]
                    edge_colors.append(e_color)
                    g.add_edges([tuple(sorted([x - 1 for x in ue]))])

                for (ref_graph, ref_colors), graph_list in iso_groups:
                    v_maps = ref_graph.get_isomorphisms_vf2(g, edge_color1=ref_colors, edge_color2=edge_colors)
                    if v_maps != [] and merge_isomorphic_graphs:
                        v_map = v_maps[0]
                        if verbose: 
                            logger.info("\033[1m{} <- {}\033[0m".format(graph_list[0].name,graph.name))

                        # map the signature from the new graph to the reference graph using the vertex map
                        edge_map = [next(re.index for re in ref_graph.es if (re.source, re.target) == (v_map[e.source], v_map[e.target]) or \
                            (re.source, re.target) == (v_map[e.target], v_map[e.source])) for e in g.es]

                        # go through all loop momenta and make the matrix
                        # and the shifts
                        mat = []
                        shifts = []
                        selected_edges = []
                        n_loops = len(graph_list[0].edges) - len(graph_list[0].nodes) + 1
                        for loop_var in range(n_loops):
                            # find the edge that has (+-)k`loop_var`
                            # TODO: Fix for higher degeneracy
                            orig_edge = next(ee for ee in graph_list[0].edges.values() if all(s == 0 for s in ee['signature'][1]) and \
                                sum(abs(s) for s in ee['signature'][0]) == 1 and ee['signature'][0][loop_var] != 0)
                            orig_edge_in_ref = next(ee for ee in ref_graph.es if tuple(sorted(orig_edge['vertices'])) == (ee.source + 1, ee.target + 1))

                            # now get the corresponding edge(s) in the new graph
                            en = g.es[edge_map.index(orig_edge_in_ref.index)]
                            new_edges = [ee for ee in graph.edges.values() \
                                if (ee['vertices'] == (en.source + 1, en.target + 1) or ee['vertices'] == (en.target + 1, en.source + 1)) and abs(ee['PDG'])== abs(orig_edge['PDG']) and ee not in selected_edges]

                            # in the case of a multi-edge, we simply take the one with the shortest momentum string
                            new_edge = sorted(new_edges, key=lambda x: len(x['momentum']))[0]
                            selected_edges.append(new_edge)
                            new_edge['name'] += 'LMB'

                            orientation_factor = 1
                            if model.get_particle(abs(new_edge['PDG']))['self_antipart']:
                                if v_map[new_edge['vertices'][0]-1]+1 != orig_edge['vertices'][0]:
                                    orientation_factor = -1
                            else:
                                match_direction = v_map[new_edge['vertices'][0]-1]+1 == orig_edge['vertices'][0]
                                match_particle = new_edge['PDG'] == orig_edge['PDG']
                                if not match_direction:# != match_particle:
                                    orientation_factor = -1

                            shifts.append([-orientation_factor * s * orig_edge['signature'][0][loop_var] for s in new_edge['signature'][1]])
                            mat.append([orientation_factor * s * orig_edge['signature'][0][loop_var] for s in new_edge['signature'][0]])

                        # Map all signatures and momentum labels of the new graph
                        shifts = np.array(shifts).T
                        mat = np.array(np.matrix(mat).T.I)
                        for e in graph.edges.values():
                            # skip external edges
                            if e["type"] != "virtual":
                                continue

                            new_sig = mat.dot(np.array(e['signature'][0]))
                            new_shift = e['signature'][1]+shifts.dot(new_sig)

                            e['signature'] = (list(new_sig), list(new_shift))
                            e['momentum'] = graph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)

                        if verbose:
                            parser = re.compile(r'(^|\+|-)(k|p)(\d*)')
                            n_incoming = sum([1 for edge in graph.edges.values() if edge['type'] == 'in'])
                            n_loops = len(graph.edges) - len(graph.nodes) + 1
                            old_momenta =['k%d'%(i+1) for i in range(n_loops)] 
                            new_momenta =[] 
                            for momentum in old_momenta:
                                parsed = [e.groups() for e in parser.finditer(momentum)]
                                signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
                                for pa in parsed:
                                    if pa[1] == 'p':
                                        signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                    else:
                                        signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                new_sig = mat.dot(np.array(signature[0]))
                                new_shift = signature[1]+shifts.dot(new_sig)
                                new_momenta.append(graph.momenta_decomposition_to_string((list(new_sig), list(new_shift)), set_outgoing_equal_to_incoming=False))
                            logger.info("Map: {} -> {}".format(tuple(old_momenta), tuple(new_momenta)))


                        # Map vertices
                        for v in graph.nodes.values():
                            if v['vertex_id']<0:
                                continue
                            # skip external edge        
                            n_incoming = sum([1 for edge in graph.edges.values() if edge['type'] == 'in'])
                            n_loops = len(graph.edges) - len(graph.nodes) + 1

                            # parse momentum
                            parser = re.compile(r'(^|\+|-)(k|p)(\d*)')
                            new_momenta = []
                            for momentum in v['momenta']:
                                parsed = [e.groups() for e in parser.finditer(momentum)]
                                signature = ([0 for _ in range(n_loops)], [0 for _ in range(n_incoming)])
                                for pa in parsed:
                                    if pa[1] == 'p':
                                        signature[1][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                    else:
                                        signature[0][int(pa[2]) - 1] = 1 if pa[0] == '' or pa[0] == '+' else -1
                                new_sig = mat.dot(np.array(signature[0]))
                                new_shift = signature[1]+shifts.dot(new_sig)
                                new_momenta.append(graph.momenta_decomposition_to_string((list(new_sig), list(new_shift)), set_outgoing_equal_to_incoming=False))
                            v['momenta'] = tuple(new_momenta)
                        graph_list.append(graph)
                        break
                else:
                    iso_groups.append(((g, edge_colors), [graph]))
                
                bar.update(g_id+1)
                bar.update(graph_nr=g_id+1)
                bar.update(iso_size=len(iso_groups))
        


        logger.info("\033[1m{} unique supergraphs\033[m".format(len(iso_groups)))
        if not cuts is None:
            if FORM_processing_options["cores"] == 1:
                cut_it = map(FORMSuperGraph.filter_valid_cuts_helper, 
                    list((iso_graph, cuts) for iso_graph in iso_groups))
            else:
                pool = multiprocessing.Pool(processes=FORM_processing_options["cores"])
                cut_it = pool.imap(FORMSuperGraph.filter_valid_cuts_helper, 
                    list((iso_graph, cuts) for iso_graph in iso_groups))
            
            graph_filtered = {'DUMP':[], 'KEEP':[]}
            with progressbar.ProgressBar(prefix='Filter SG with valid cuts: {variables.keep}\u2713  {variables.drop}\u2717 : ', max_value=len(iso_groups),variables={'keep': '0', 'drop':'0'}) as bar:
                for sgid, (iso_graph, valid_cutQ) in enumerate(cut_it):
                    if valid_cutQ: 
                        graph_filtered["KEEP"] += [iso_graph]
                    else: 
                        graph_filtered["DUMP"] += [iso_graph]
                    bar.update(sgid+1)
                    bar.update(keep=len(graph_filtered['KEEP']))
                    bar.update(drop=len(graph_filtered['DUMP']))
            #for k,v in graph_filtered.items():
            #    print(k,":")
            #    print("\t",np.array([g.name for g in v]))
            #full_graph_list[:] = graph_filtered['KEEP']
            iso_groups[:] = graph_filtered['KEEP']
            logger.info("\033[1mRemoved {} graphs with no valid cuts\033[0m".format(len(graph_filtered['DUMP'])))

        FORM_iso_sg_list = FORMSuperGraphList([FORMSuperGraphIsomorphicList(iso_group) for _, iso_group in iso_groups], name=p.stem)

        if workspace is not None:
            selected_workspace = workspace
            shutil.copy(pjoin(plugin_path,"multiplicity.frm"),pjoin(selected_workspace,'multiplicity.frm'))
            shutil.copy(pjoin(plugin_path,"numerator.frm"),pjoin(selected_workspace,'numerator.frm'))
            shutil.copy(pjoin(plugin_path,"raw_numerator.frm"),pjoin(selected_workspace,'raw_numerator.frm'))
            shutil.copy(pjoin(plugin_path,"tensorreduce.frm"),pjoin(selected_workspace,'tensorreduce.frm'))
            shutil.copy(pjoin(plugin_path,"Gstring.prc"),pjoin(selected_workspace,'Gstring.prc'))
            shutil.copy(pjoin(plugin_path,"integrateduv.frm"),pjoin(selected_workspace,'integrateduv.frm'))
            shutil.copy(pjoin(plugin_path,"diacolor.h"),pjoin(selected_workspace,'diacolor.h'))
            FORM_source = pjoin(selected_workspace,'multiplicity.frm')

            if FORM_processing_options["cores"] == 1:
                graph_it = map(FORMSuperGraphIsomorphicList.multiplicity_factor_helper, 
                    (list((iso_graphs, iso_id, selected_workspace, FORM_source))
                    for iso_id, iso_graphs in enumerate(FORM_iso_sg_list)))
            else:
                pool = multiprocessing.Pool(processes=FORM_processing_options["cores"])
                graph_it = pool.imap(FORMSuperGraphIsomorphicList.multiplicity_factor_helper, 
                    (list((iso_graphs, iso_id, selected_workspace, FORM_source))
                    for iso_id, iso_graphs in enumerate(FORM_iso_sg_list)))

            with progressbar.ProgressBar(
                prefix = 'Processing Isomorphic sets (ISO #{variables.processed}/%d, Zeros: {variables.zero_multiplicity_count})'%len(iso_groups),
                max_value=len(iso_groups),
                variables = {'processed' : '0', 'zero_multiplicity_count' : '0'}
                ) as bar:
            
                zero_multiplicity_count = 0
                processed_graphs = 0

                for iso_graphs,(isogs, multiplicity) in zip(FORM_iso_sg_list, graph_it):
                    if iso_graphs[0].name != isogs[0].name:
                        raise FormProcessingError("Failed to assign the multiplicity factor to the correct graph!")
                    processed_graphs += 1
                    bar.update(processed='%d'%processed_graphs, zero_multiplicity_count='%d'%zero_multiplicity_count)
                    if multiplicity == 0:
                        zero_multiplicity_count += 1
                        iso_graphs[:] = []
                    else:
                        iso_graphs[:] = iso_graphs[:1]
                        iso_graphs[0].multiplicity = multiplicity
                    bar.update(processed_graphs)
                bar.update(zero_multiplicity_count='%d'%zero_multiplicity_count)

                for iso_id in reversed(range(len(FORM_iso_sg_list))):
                    if FORM_iso_sg_list[iso_id] == []:
                        del FORM_iso_sg_list[iso_id]

            if FORM_processing_options["cores"] > 1:
                pool.close()
                
            if zero_multiplicity_count > 0 :
                logger.info("%d isomorphic sets have 0 multiplicity -> they are dropped!"%zero_multiplicity_count)
            
            #for iso_id, iso_graphs in enumerate(FORM_iso_sg_list):
            #    g = iso_graphs[0]
            #    print("#{}\t{:10}: {}".format(iso_id, g.name, g.multiplicity))

        # Generate additional copies for different LMB of the representative graph of each isomorphic set.
        # Note that depending on the option FORM_processing['representative_lmb'], the function below can also modify the LMB of the representative sg.
        for FORM_iso_sg in FORM_iso_sg_list:
            FORM_iso_sg[0].adjust_LMBs(model)

        # Export the drawings corresponding to each ISO supergraphs
        for i_graph, super_graphs in enumerate(FORM_iso_sg_list):
            super_graphs[0].draw(model, pjoin(workspace, '../../Drawings'), FORM_id=i_graph)
            # Draw supergraphs for additional LMBs
            if isinstance(super_graphs[0].additional_lmbs,list):
                for i_lmb,_,_,sg in super_graphs[0].additional_lmbs:
                    sg.draw(model, pjoin(workspace, '../../Drawings'), FORM_id=i_graph, lmb_id=i_lmb)
        return FORM_iso_sg_list

    def to_dict(self, file_path):
        """ Outputs the FORMSuperGraph list to a Python dict file path."""

        # TODO: Dump to Python dict
        pass
    
    @classmethod
    def minimal_dict_representation_of_a_topology_generator_graph(cls, topology_generator_graph, particle_PDGs):
        
        g = topology_generator_graph
        sig_map=g.get_signature_map()
        #assert(all(all(sm_i==0 for sm_i in sm[1][len(sm[1])//2:]) for sm in sig_map.values()))
        edges_description = []
        for e in g.edge_map_lin:
            edges_description.append({
                'name' : e[0],
                'PDG' : particle_PDGs[e[0]],
                'start_node' : e[1],
                'end_node' : e[2],
                'power' : (g.powers[e[0]] if e[0] in g.powers else 0),
                'loop_momentum_signature' : [list(sig_map[e[0]][0]),list(sig_map[e[0]][1])]
            })
        return {
            'massless_denominators_string' : '' if g.denominators_string is None else g.denominators_string,
            'loop_momenta' : [g.edge_map_lin[e][0] for e in g.loop_momenta],
            'edges' : edges_description
        }

    @classmethod
    def process_raw_numerator(cls,input_numerator_path):
        with open(input_numerator_path,'r') as f:
            return f.read().replace('rat','')

    def write_pySecDec_run(self, output_dir, numerator_input_path, topology_generator_graph, graph_name, graph_id, couplings_prefactor_str, additional_overall_factor, params, model, particle_PDGs):
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

        # Add gs which is typically not present as it is meant to be a dynamical parameter
        processed_params['gs'] = math.sqrt(model['parameter_dict']['aS'].real*(4.*math.pi))

        output_path = pjoin(output_dir,'run_graph_{}.py'.format(graph_id))
        output_numerator_path = pjoin(output_dir,'numerator_{}.txt'.format(graph_id))

        template_run = None
        with open(pjoin(plugin_path,'Templates','run_pySecDec_template.py'),'r') as f:
            template_run = f.read()

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
                    mass_param = model.get_particle(particle_PDGs[e[0]]).get('mass')
                    if mass_param.upper()!='ZERO':
                        default_externals.append( (model['parameter_dict'][mass_param].real, 0., 0., 0.) )
                    else:
                        # If there is a single incoming, do not put it onshell
                        if len(g.ext)==2:
                            default_externals.append( (1., 0., 0., 0.) )
                        else:
                            default_externals.append( (1., 0., 0., ((-1)**len(default_externals))*1.) )
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
                mass_param = model.get_particle(particle_PDGs[e[0]]).get('mass')
                if mass_param.upper()=='ZERO':
                    mass_str = ''
                else:
                    if mass_param not in real_parameters:
                        real_parameters.append(mass_param)
                        real_parameters_input.append(model['parameter_dict'][mass_param].real)
                    mass_str = '-%s**2'%mass_param

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

        # Drawing replacement
        repl_dict['drawing_input_internal_lines'] = str(internal_edges)
        repl_dict['drawing_input_power_list'] = str(power_list)
        repl_dict['drawing_input_external_lines'] = str(external_edges)

        # Loop package replacement
        repl_dict['propagators'] = str(propagators)
        repl_dict['loop_momenta'] = str(loop_momenta_str)
        repl_dict['external_momenta'] = str(external_momenta_str)
        repl_dict['lorentz_indices'] = str(lorentz_indices)
        repl_dict['power_list'] = str(power_list)
        repl_dict['numerator_path'] = './%s'%(os.path.basename(output_numerator_path))
        repl_dict['replacement_rules'] = str(replacement_rules)
        repl_dict['real_parameters'] = str(real_parameters)
        repl_dict['loop_additional_prefactor'] = '( I*(4*pi)**(-2+eps) )**(%d)'%n_loops
        repl_dict['contour_deformation'] = 'True'
        repl_dict['max_epsilon_order'] = 0

        # pySecDec integration replacement
        repl_dict['n_loops'] = str(g.n_loops)
        repl_dict['additional_overall_factor'] = '1.'+str(additional_overall_factor)
        repl_dict['real_parameters_input'] = str(real_parameters_input)
        repl_dict['complex_parameters_input'] = str([])
        repl_dict['couplings_prefactor'] = str(couplings_prefactor_str)
        repl_dict['couplings_values'] = str(processed_params)
        with open(output_path,'w') as f:
            f.write(template_run.format(**repl_dict))

        # Make the file executable
        os.chmod(output_path, os.stat(output_path).st_mode | stat.S_IEXEC)

    def generate_integrand_functions(self, root_output_path, model, additional_overall_factor='',
                                    params={}, output_format='c', workspace=None, header="", integrand_type=None, process_definition=None):
        header_map = {'header': header}
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if integrand_type is None:
            return

        if len(self)==0:
            raise FormProcessingError("Cannot generat numerators for an empty list of supergraphs.")

        if output_format not in self.extension_names:
            raise FormProcessingError("This FORMSuperGraphList instance requires at least one entry for generating numerators.")

        # add all numerators in one file and write the headers
        numerator_header = """
#include <complex>
#include <signal.h>
#include "%(header)snumerator.h"
#include "dual.h"
#include "dualt2.h"
#include "dualkt2.h"
#include "dualt3.h"
#include "dualkt3.h"
#include "dualklt3.h"
#include <mp++/complex128.hpp>

using namespace std;
using namespace std::complex_literals;
using namespace duals;
using namespace duals::literals;
using namespace dualst2;
using namespace dualskt2;
using namespace dualst3;
using namespace dualskt3;
using namespace dualsklt3;
using namespace mppp;
using namespace mppp::literals;

const complex<double> I{ 0.0, 1.0 };

"""

        var_pattern = re.compile(r'Z\d*_')
        lm_pattern = re.compile(r'lm(\d*)')
        energy_pattern = re.compile(r'E(\d+)')
        denom_pattern = re.compile(r'invd(\d+)')
        conf_exp = re.compile(r'conf\(([^)]*)\)(\*tder\([^)]*\))?\n')
        return_exp = re.compile(r'return ([^;]*);\n')
        float_pattern = re.compile(r'((\d+\.\d*)|(\.\d+))')
        diag_pattern = re.compile(r'diag\((\d*)\)')
        forest_pattern = re.compile(r'forestid\((\d*)\)')
        square_pattern = re.compile(r'([^^*+\-\s]*)\^2')
        pow_pattern = re.compile(r'([^^*+\-\s]*)\^(\d+)')
        int_pattern = re.compile(r'(?:^|(?<=[^\.\w]))(\d+)(?=[^\.\w]|$)')
        empty_line_pattern = re.compile(r'\n\n\n')

        integrand_type_list = [integrand_type] if integrand_type != 'both' else ['PF', 'LTD']

        all_numerator_ids = {'PF': [], 'LTD': []}

        if output_format in ['pySecDec','raw']:
            if not os.path.isdir(pjoin(root_output_path, 'raw_SG_expressions')):
                os.mkdir(pjoin(root_output_path, 'raw_SG_expressions'))

        if output_format == 'pySecDec':
            if not os.path.isdir(pjoin(root_output_path, 'pySecDec_runs')):
                os.mkdir(pjoin(root_output_path, 'pySecDec_runs'))

        for itype in integrand_type_list:
            with progressbar.ProgressBar(
                prefix = 'Processing %s integrand (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}) with FORM ({variables.timing} ms / supergraph) : '%(itype, len(self)),
                max_value=len(self),
                variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb': '0'}
                ) as bar:
                total_time = 0.

                for i_graph, graph in enumerate(self):
                    graph.is_zero = True

                    graphs_to_process = []
                    if isinstance(graph[0].additional_lmbs,list) and graph[0].additional_lmbs != []:
                        graphs_to_process.append( (0,i_graph,graph[0]) )
                        graphs_to_process.extend([(g.additional_lmbs,i_graph,g) for _,_,_,g in graph[0].additional_lmbs])
                    else:
                        # By setting the active graph to None we will then sum overall members of the iso set.
                        graphs_to_process.append( (0,i_graph, None) )
                    bar.update(max_lmb='%d'%len(graphs_to_process))
                    for i_lmb, i_g, active_graph in graphs_to_process:
                        bar.update(i_graph='%d'%(i_graph+1), i_lmb='%d'%(i_lmb+1))
                        i = i_lmb*FORM_processing_options['FORM_call_sig_id_offset_for_additional_lmb']+i_g
                        all_numerator_ids[itype].append(i)
                        time_before = time.time()

                        if output_format in ['raw','pySecDec']:

                            if active_graph is None:
                                characteristic_graph = graph[0]
                            else:
                                characteristic_graph = active_graph
                            
                            particle_ids = { e['name']: e['PDG'] for e in characteristic_graph.edges.values() }
                            with open(pjoin(root_output_path, 'raw_SG_expressions', 'numerator_{}.txt'.format(i)), 'w') as f:
                                f.write(self.process_raw_numerator(pjoin(root_output_path, 'workspace', 'out_raw_numerator_{}.txt'.format(i))))

                            couplings_prefactor = '1'
                            with open(pjoin(root_output_path, 'workspace', 'out_raw_prefactor_{}.txt'.format(i)),'r') as f:
                                couplings_prefactor = f.read().strip()


                            characteristic_graph_dict = self.minimal_dict_representation_of_a_topology_generator_graph(characteristic_graph.squared_topology.topo,particle_ids)
                            characteristic_graph_dict['numerator']='numerator_{}.txt'.format(i)
                            characteristic_graph_dict['name']=characteristic_graph.squared_topology.name
                            characteristic_graph_dict['couplings_prefactor']=couplings_prefactor
                            with open(pjoin(root_output_path, 'raw_SG_expressions', 'graph_{}.yaml'.format(i)), 'w') as f:
                                f.write(yaml.dump(characteristic_graph_dict, Dumper=Dumper))

                            if output_format == 'pySecDec':
                                
                                self.write_pySecDec_run(
                                    pjoin(root_output_path, 'pySecDec_runs'),
                                    pjoin(root_output_path, 'workspace', 'out_raw_numerator_{}.txt'.format(i)),
                                    characteristic_graph.squared_topology.topo,
                                    characteristic_graph.squared_topology.name,
                                    i,
                                    couplings_prefactor,
                                    additional_overall_factor,
                                    params,
                                    model,
                                    particle_ids
                                )

                            continue

                        with open(pjoin(root_output_path, 'workspace', 'out_integrand_{}_{}.proto_c'.format(itype, i))) as f:
                            num = f.read()
                        

                        # TODO Remove when FORM will have fixed its C output bug
                        num = temporary_fix_FORM_output(num)

                        total_time += time.time()-time_before
                        num = num.replace('i_', 'I')
                        num = num.replace('\nZ', '\n\tZ') # nicer indentation
                        num = num.replace('\\\n', '') # remove breaks in long numbers
                        num = re.sub(r'(\d)\n *\.', r'\1.', num) # work around form output bug where the . is disconnected from the floating number
                        num = re.sub(r'\(([^)]*)\n *', r'(\1', num) # connect functions over a newline

                        energies_per_cut = {}
                        propagators_per_cut = {}
                        confs = []
                        integrand_main_code = ''
                        integrand_f128_main_code = ''
                        integrand_mpfr_main_code = ''
                        conf_secs = num.split('#CONF')
                        for conf_sec in conf_secs[1:]:
                            conf_sec = conf_sec.replace("#CONF\n", '')
                            conf_part = list(conf_exp.finditer(conf_sec))[0].groups()

                            # parse the configuration id
                            conf = conf_part[0].split(',')
                            conf_sec = conf_exp.sub('', conf_sec)
                            conf_sec = conf_sec.replace('\n    ', "\n\t")
                            conf_sec = conf_sec.replace('\n   ', "\n\t")

                            denominator_mode = 'NUM' if int(conf[0]) >= 0 else ('FOREST' if len(conf) == 2 else 'DIAG')
                            dual_num_mode = conf_part[1] is not None
                            if conf_part[1] == "*tder(2)":
                                dual_num_type = "dual"
                            elif conf_part[1] == "*tder(3)":
                                dual_num_type = "dualt2"
                            elif conf_part[1] == "*tder(2,2)":
                                dual_num_type = "dualkt2"
                            elif conf_part[1] == "*tder(4)":
                                dual_num_type = "dualt3"
                            elif conf_part[1] == "*tder(2,3)":
                                dual_num_type = "dualkt3"
                            elif conf_part[1] == "*tder(2,2,2)":
                                dual_num_type = "dualklt3"
                            elif not dual_num_mode:
                                dual_num_type = ""
                                pass
                            else:
                                raise AssertionError("Unknown dual configuration: %s" % conf_part[1])
                            
                            if denominator_mode == 'DIAG':
                                cut_id = int(conf[1])

                                # parse the constants
                                const_secs = conf_sec.split('#CONSTANTS')[1]
                                const_secs = const_secs.replace('\n', '').replace(' ', '').replace('\t', '')
                                const_code = '*'.join('({})'.format(e) for e in const_secs.split(',') if e != '')

                                # parse the energies
                                energy_secs = conf_sec.split('#ENERGIES')[1]
                                energy_secs = energy_secs.replace('\n', '').replace(' ', '').replace('\t', '')

                                e = [e for e in energy_secs.split(',') if e != '']
                                assert(len(e) % 2 == 0)
                                for j in range(0, len(e), 2):
                                    energy_index = int(e[j][1:]) # skip E prefix
                                    energy_instr = 'sqrt({})'.format(lm_pattern.sub(r'lm[\1]', square_pattern.sub(r'\1*\1', e[j+1]))).replace('2*', '2.*').replace('4*', '4.*').replace('8*', '8.*')
                                    assert((cut_id, energy_index) not in energies_per_cut or energies_per_cut[(cut_id, energy_index)] == energy_instr)
                                    energies_per_cut[(cut_id, energy_index)] = energy_instr

                                const_code = pow_pattern.sub(r'pow(\1,\2)', const_code).replace('(1)', '(1.)')
                                const_code = square_pattern.sub(r'\1*\1', const_code)

                                # convert ints to floats to prevent mpfr and dual number conversion issues
                                const_code = re.sub(r'([\+\-(])(\d+)\*', r'\1\2.*', const_code)
                                if const_code == '':
                                    const_code = '1.'

                                # parse the denominators
                                denom_secs = conf_sec.split('#ELLIPSOIDS')[1]
                                denom_secs = denom_secs.replace('\n', '').replace('\t', '')
                                d = [d for d in denom_secs.split(',') if d != '']

                                for j in range(0, len(d), 2):
                                    denom_index = int(d[j][4:]) # strip invd prefix
                                    denom_instr = '1./({})'.format(energy_pattern.sub(r'E[\1]', lm_pattern.sub(r'lm[\1]', int_pattern.sub(r'\1.', d[j+1]))))
                                    assert((cut_id, denom_index) not in propagators_per_cut or propagators_per_cut[(cut_id, denom_index)] == denom_instr)
                                    propagators_per_cut[(cut_id, denom_index)] = denom_instr

                                conf_sec = conf_sec.split('#ELLIPSOIDS')[-1]

                            const_code = denom_pattern.sub(r'invd[\1]', energy_pattern.sub(r'E[\1]', lm_pattern.sub(r'lm[\1]', const_code)))
                            conf_sec = denom_pattern.sub(r'invd[\1]', energy_pattern.sub(r'E[\1]', lm_pattern.sub(r'lm[\1]', conf_sec)))
                            conf_sec = re.sub(r'(^|[+\-*=])(\s*\d+)($|[^.\d\w])', r'\1\2.\3', conf_sec)
                            returnval = list(return_exp.finditer(conf_sec))[0].groups()[0]
                            returnval = re.sub(r'(^|[+\-*=]|\*\()(\s*\d+)($|[^.\d\w])', r'\1\2.\3', returnval)

                            if dual_num_mode and all(x not in returnval for x in ['E', 'lm', 'Z']):
                                returnval = '{}({})'.format(dual_num_type, returnval)

                            if denominator_mode == 'FOREST':
                                conf_sec = return_exp.sub('return {};\n'.format(returnval), conf_sec)
                            elif denominator_mode == 'DIAG':
                                conf_sec = return_exp.sub('return 1./({})*({});\n'.format(const_code, returnval), conf_sec)
                            else:
                                conf_sec = return_exp.sub('*out = {};\n'.format(returnval), conf_sec)

                            temp_vars = list(sorted(set(var_pattern.findall(conf_sec))))

                            # write out all integer powers as multiplications to prevent slow pow evaluation with floating exponent
                            conf_sec = re.sub(r'pow\(([^,]+),(\d+)\)', lambda x: '*'.join([x.group(1)]*int(x.group(2))) , conf_sec)

                            base_type = 'complex<double>'
                            dual_base_type = '{}<{}>'.format(dual_num_type, base_type) if dual_num_mode else base_type
                            base_type_f128 = 'complex128'
                            dual_base_type_f128 = '{}<{}>'.format(dual_num_type,base_type_f128) if dual_num_mode else base_type_f128
                            base_type_mpfr = 'mpcomplex'
                            dual_base_type_mpfr = '{}<{}>'.format(dual_num_type,base_type_mpfr) if dual_num_mode else base_type_mpfr
                            

                            if denominator_mode == 'FOREST':
                                main_code = conf_sec.replace('logmUVmu', 'log(mUV*mUV/(mu*mu))').replace('logmUV', 'log(mUV*mUV)').replace('logmu' , 'log(mu*mu)').replace('logmt' , 'log(masst*masst)')
                                main_code_with_diag_call = diag_pattern.sub(r'diag_\1(lm, params, E, invd)', main_code)
                                integrand_main_code += '\nstatic {0} forest_{2}({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(
                                    dual_base_type, base_type, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type, ','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_with_diag_call
                                )

                                main_code_f128 = main_code.replace('pi', 'mppp::real128_pi()').replace('complex<double>', 'complex128')
                                main_code_f128 = float_pattern.sub(r'real128(\1q)', main_code_f128)
                                main_code_f128 = diag_pattern.sub(r'diag_\1_f128(lm, params, E, invd)', main_code_f128)
                                integrand_f128_main_code += '\n' + '\nstatic {0} forest_{2}_f128({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(
                                    dual_base_type_f128, base_type_f128, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type_f128, ','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_f128
                                )

                                main_code_mpfr = main_code.replace('pi', 'mpreal(mpfr::const_pi())').replace('complex<double>', 'mpcomplex')
                                main_code_mpfr = float_pattern.sub(r'mpreal("\1")', main_code_mpfr)
                                main_code_mpfr = re.sub(r'sqrt\(([^)]+)\)', r'sqrt(mpcomplex(\1))', main_code_mpfr)
                                main_code_mpfr = re.sub(r'pow\(([^,]+)', r'pow(mpcomplex(\1)', main_code_mpfr)
                                main_code_mpfr = diag_pattern.sub(r'diag_\1_mpfr(lm, params, E, invd)', main_code_mpfr)
                                integrand_mpfr_main_code += '\n' + '\nstatic {0} forest_{2}_mpfr({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(
                                    dual_base_type_mpfr, base_type_mpfr, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type_mpfr, ','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_mpfr
                                )
                            elif denominator_mode == 'DIAG':
                                main_code = conf_sec
                                main_code = main_code.replace('logmUVmu', 'log(mUV*mUV/(mu*mu))').replace('logmUV', 'log(mUV*mUV)').replace('logmu' , 'log(mu*mu)').replace('logmt' , 'log(masst*masst)')
                                integrand_main_code += '\nstatic {0} diag_{2}({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(dual_base_type, base_type, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type, ','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code
                                )

                                main_code_f128 = main_code.replace('pi', 'mppp::real128_pi()').replace('complex<double>', 'complex128')
                                main_code_f128 = float_pattern.sub(r'real128(\1q)', main_code_f128)
                                integrand_f128_main_code += '\n' + '\nstatic {0} diag_{2}_f128({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(dual_base_type_f128, base_type_f128, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type_f128,','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_f128
                                )

                                main_code_mpfr = main_code.replace('pi', 'mpreal(mpfr::const_pi())').replace('complex<double>', 'mpcomplex')
                                main_code_mpfr = re.sub(r'sqrt\(([^)]+)\)', r'sqrt(mpcomplex(\1))', main_code_mpfr)
                                main_code_mpfr = re.sub(r'pow\(([^,]+)', r'pow(mpcomplex(\1)', main_code_mpfr)
                                main_code_mpfr = re.sub(r'log\(([^)]+)\)', r'log(mpcomplex(\1))', main_code_mpfr)
                                main_code_mpfr = float_pattern.sub(r'mpreal("\1")', main_code_mpfr)
                                integrand_mpfr_main_code += '\n' + '\nstatic {0} diag_{2}_mpfr({0} lm[], {1} params[], {0} E[], {0} invd[]) {{{3}\n{4}}}'.format(dual_base_type_mpfr, base_type_mpfr, abs(int(conf[0])),
                                    '\n\t{} {};'.format(dual_base_type_mpfr, ','.join(temp_vars)) if len(temp_vars) > 0 else '', main_code_mpfr
                                )
                            else:
                                cut_id = int(conf[0])
                                confs.append((cut_id, dual_num_type))

                                forests = []
                                if len(conf_sec) > 0:
                                    graph.is_zero = False

                                if len(temp_vars) == 0:
                                    # if generating in format C mode, construct the forests only once
                                    forests = list(sorted(set(forest_pattern.findall(conf_sec))))

                                # generate the code for the energies and invds
                                max_energy = max((eid + 1 for (cid, eid) in energies_per_cut if cid == cut_id), default=0)
                                energies = [('0' if not dual_num_mode else '{}()'.format(dual_base_type)) if (cut_id, jj) not in energies_per_cut else energies_per_cut[(cut_id, jj)] for jj in range(max_energy)]
                                max_denom = max((eid + 1 for (cid, eid) in propagators_per_cut if cid == cut_id), default=0)
                                denoms = [('0' if not dual_num_mode else '{}()'.format(dual_base_type)) if (cut_id, jj) not in propagators_per_cut else propagators_per_cut[(cut_id, jj)] for jj in range(max_denom)]

                                conf_sec = '\n\t{} E[] = {{{}}};\n'.format(dual_base_type, ('0' if not dual_num_mode else '{}()'.format(dual_base_type)) if len(energies) == 0 else ','.join(energies)) + \
                                           '\t{} invd[] = {{{}}};\n'.format(dual_base_type, ('0' if not dual_num_mode else '{}()'.format(dual_base_type)) if len(denoms) == 0 else ','.join(denoms)) + \
                                           ('\t{} {};'.format(dual_base_type, ','.join(temp_vars)) if len(temp_vars) > 0 else '') + \
                                           ('\n' + '\n'.join('\t{0} forest{1}=forest_{1}(lm, params, E, invd);'.format(dual_base_type, f) for f in forests) if len(forests) > 0 else '') + conf_sec

                                main_code = conf_sec.replace('logmUVmu', 'log(mUV*mUV/(mu*mu))').replace('logmUV', 'log(mUV*mUV)').replace('logmu' , 'log(mu*mu)').replace('logmt' , 'log(masst*masst)')
                                if len(temp_vars) == 0:
                                    main_code = forest_pattern.sub(r'forest\1', main_code)
                                else:
                                    main_code = forest_pattern.sub(r'forest_\1(lm, params, E, invd)', main_code)

                                integrand_main_code += '\nstatic inline void %(header)sevaluate_{2}_{3}_{4}({0} lm[], {1} params[], {0}* out) {{{5}}}'.format(dual_base_type, base_type, itype, i, int(conf[0]),
                                    diag_pattern.sub(r'diag_\1(lm, params, E, invd)', main_code)
                                )

                                main_code_f128 = main_code.replace('pi', 'mppp::real128_pi()').replace('complex<double>', 'complex128')
                                main_code_f128 = float_pattern.sub(r'real128(\1q)', main_code_f128)
                                main_code_f128 = main_code_f128.replace('(lm,', '_f128(lm,')
                                integrand_f128_main_code += '\n' + '\nstatic inline void %(header)sevaluate_{2}_{3}_{4}_f128({0} lm[], {1} params[], {0}* out) {{{5}}}'.format(dual_base_type_f128, base_type_f128, itype, i, int(conf[0]),
                                    diag_pattern.sub(r'diag_\1_f128(lm, params, E, invd)', main_code_f128)
                                )

                                main_code_mpfr = main_code.replace('pi', 'mpreal(mpfr::const_pi())').replace('complex<double>', 'mpcomplex')

                                for p in params:
                                    if p != 'pi':
                                        main_code_mpfr = main_code_mpfr.replace(p, 'mpcomplex({})'.format(p))
                                main_code_mpfr = re.sub(r'pow\(([^,]+)', r'pow(mpcomplex(\1)', main_code_mpfr)
                                main_code_mpfr = re.sub(r'sqrt\(([^)]+)\)', r'sqrt(mpcomplex(\1))', main_code_mpfr)
                                main_code_mpfr = re.sub(r'log\(([^)]+)\)', r'log(mpcomplex(\1))', main_code_mpfr)
                                main_code_mpfr = float_pattern.sub(r'mpreal("\1")', main_code_mpfr)
                                main_code_mpfr = forest_pattern.sub(r'forest_\1_mpfr(lm, params, E, invd)', main_code_mpfr)
                                main_code_mpfr = re.sub(r'\*out =([^;]*);', r'*out = (complex128)(\1);', main_code_mpfr)
                                integrand_mpfr_main_code += '\n' + '\nstatic inline void %(header)sevaluate_{2}_{3}_{4}_mpfr({0} lm[], {1} params[], {0}* out) {{{5}}}'.format(dual_base_type_mpfr, base_type_mpfr, itype, i, int(conf[0]),
                                    diag_pattern.sub(r'diag_\1_mpfr(lm, params, E, invd)', main_code_mpfr)
                                )

                        
                        for x in range(4):
                            integrand_main_code = empty_line_pattern.sub(r'\n\n', integrand_main_code)
                            integrand_f128_main_code = empty_line_pattern.sub(r'\n\n', integrand_f128_main_code)
                            integrand_mpfr_main_code = empty_line_pattern.sub(r'\n\n', integrand_mpfr_main_code)

                        integrand_main_code = numerator_header + integrand_main_code
                        integrand_f128_main_code = numerator_header + integrand_f128_main_code
                        integrand_mpfr_main_code = numerator_header + integrand_mpfr_main_code

                        integrand_f128_main_code = integrand_f128_main_code.replace('const complex<double> I{ 0.0, 1.0 };', 'constexpr complex128 I{ 0.0, 1.0 };')


                        conf_no_dual, conf_dual = [x for x in confs if x[1]] != "", [x for x in confs if x[1] == ""]

                        integrand_main_code += \
"""
extern "C" {{
void %(header)sevaluate_{0}_{1}(complex<double> lm[], complex<double> params[], int conf, complex<double>* out) {{
   switch(conf) {{
{2}
    }}
}}

void %(header)sevaluate_{0}_{1}_dual(double lm[], complex<double> params[], int conf, double* out) {{
   switch(conf) {{
{3}
    }}
}}
}}
""".format(itype, i,
        '\n'.join(
            ['\t\tcase {}: %(header)sevaluate_{}_{}_{}(lm, params, out); return;'.format(conf, itype, i, conf) for conf, is_dual in sorted(x for x in confs if not x[1])] +
            (['\t\tdefault: *out = 0.;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ),
        '\n'.join(
            ['\t\tcase {}: %(header)sevaluate_{}_{}_{}(({}<complex<double>>*)lm, params, ({}<complex<double>>*)out); return;'.format(conf, itype, i, conf, is_dual, is_dual) for conf, is_dual in sorted(x for x in confs if x[1])] +
            (['\t\tdefault: *out = 0.;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        )
        )


                        integrand_f128_main_code += \
"""
extern "C" {{
void %(header)sevaluate_{0}_{1}_f128(complex128 lm[], complex128 params[], int conf, complex128* out) {{
   switch(conf) {{
{2}
    }}
}}

void %(header)sevaluate_{0}_{1}_f128_dual(complex128 lm[], complex128 params[], int conf, complex128* out) {{
   switch(conf) {{
{3}
    }}
}}
}}
""".format(itype, i,
        '\n'.join(
                ['\t\tcase {}: %(header)sevaluate_{}_{}_{}_f128(lm, params, out); return;'.format(conf, itype, i, conf) for conf, is_dual in sorted(x for x in confs if not x[1])] +
                (['\t\tdefault: *out = 0.q;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        ),
        '\n'.join(
                ['\t\tcase {}: %(header)sevaluate_{}_{}_{}_f128(({}<complex128>*)lm, params, ({}<complex128>*)out); return;'.format(conf, itype, i, conf, is_dual, is_dual) for conf, is_dual in sorted(x for x in confs if x[1])] +
                (['\t\tdefault: *out = real128(0.q);']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        )
        )



                        integrand_mpfr_main_code = \
"""
#include "mpreal.h"
#include "mpcomplex.h"

using mpfr::mpreal;
using mpfr::mpcomplex;

""" + integrand_mpfr_main_code
                        integrand_mpfr_main_code += \
"""
extern "C" {{
void %(header)sevaluate_{0}_{1}_mpfr(complex128 lm[], complex128 params[], int conf, int prec, complex128* out) {{
   mpfr_set_default_prec((mpfr_prec_t)(ceil(prec * 3.3219280948873624)));
   switch(conf) {{
{2}
    }}
}}

void %(header)sevaluate_{0}_{1}_mpfr_dual(complex128 lm[], complex128 params[], int conf, int prec, complex128* out) {{
   mpfr_set_default_prec((mpfr_prec_t)(ceil(prec * 3.3219280948873624)));
   switch(conf) {{
{3}
    }}
}}
}}
""".format(itype, i,
        '\n'.join(
            ['\t\tcase {}: %(header)sevaluate_{}_{}_{}_mpfr(lm, params, out); return;'.format(conf, itype, i, conf) for conf, is_dual in sorted(x for x in confs if not x[1])] +
            (['\t\tdefault: *out = 0.;']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else
        ),
        '\n'.join(
            ['\t\tcase {}: %(header)sevaluate_{}_{}_{}_mpfr(({}<complex128>*)lm, params, ({}<complex128>*)out); return;'.format(conf, itype, i, conf, is_dual, is_dual) for conf, is_dual in sorted(x for x in confs if x[1])] +
            (['\t\tdefault: *out = real128(0.q);']) # ['\t\tdefault: raise(SIGABRT);'] if not graph.is_zero else 
        )
        )

                        bar.update(timing='%d'%int((total_time/float(i_graph+1))*1000.0))
                        bar.update(i_graph+1)

                        writers.CPPWriter(pjoin(root_output_path, '%(header)sintegrand_{}_{}_f64.c'%header_map).format(itype, i)).write((integrand_main_code)%header_map)
                        writers.CPPWriter(pjoin(root_output_path, '%(header)sintegrand_{}_{}_f128.c'%header_map).format(itype, i)).write((integrand_f128_main_code)%header_map)

                        if FORM_processing_options["generate_arb_prec_output"]:
                            writers.CPPWriter(pjoin(root_output_path, '%(header)sintegrand_{}_{}_mpfr.cpp'%header_map).format(itype, i)).write((integrand_mpfr_main_code)%header_map)

        integrand_C_source_size = 0.
        for fpath in glob_module.glob(pjoin(root_output_path,'*integrand*')):
            integrand_C_source_size += os.path.getsize(fpath)

        self.code_generation_statistics['integrand_C_source_size_in_kB'] = int(integrand_C_source_size/1000.0)

    def generate_numerator_functions(self, root_output_path, model=None, additional_overall_factor='', params={}, 
                                    output_format='c', workspace=None, header="", integrand_type=None, process_definition=None, recycle=False):

        start_time = time.time()
        header_map = {'header': header}
        """ Generates optimised source code for the graph numerator in several
        files rooted in the specified root_output_path."""

        if len(self)==0:
            raise FormProcessingError("Cannot generat numerators for an empty list of supergraphs.")

        if output_format not in self.extension_names:
            raise FormProcessingError("Output format specified not supported. %s not in %s."%(output_format, str(self.extension_names)))

        # copy all the source files before processing the numerators
        if workspace is not None:
            shutil.copy(pjoin(plugin_path,"numerator.frm"),pjoin(workspace,'numerator.frm'))
            shutil.copy(pjoin(plugin_path,"raw_numerator.frm"),pjoin(workspace,'raw_numerator.frm'))
            shutil.copy(pjoin(plugin_path,"tensorreduce.frm"),pjoin(workspace,'tensorreduce.frm'))
            shutil.copy(pjoin(plugin_path,"Gstring.prc"),pjoin(workspace,'Gstring.prc'))
            shutil.copy(pjoin(plugin_path,"integrateduv.frm"),pjoin(workspace,'integrateduv.frm'))
            shutil.copy(pjoin(plugin_path,"diacolor.h"),pjoin(workspace,'diacolor.h'))

        max_buffer_size = 0
        all_numerator_ids = []
        processed_graph_counter = 0
        total_time = 0.
        with progressbar.ProgressBar(
            prefix = 'Processing numerators with FORM (graph {variables.i_graph}/%d, {variables.timing} ms / supergraph) : '%len(self),
            max_value=len(self),
            variables = {'timing' : 0, 'i_graph': 0}
            ) as bar:

            if FORM_processing_options["cores"] == 1:
                graph_it = map(FORMSuperGraphIsomorphicList.generate_numerator_file_helper, 
                    list((graph, i, model, root_output_path, additional_overall_factor, workspace, integrand_type, process_definition, header_map, output_format, recycle, FORM_processing_options)
                    for i, graph in enumerate(self)))
            else:
                pool = multiprocessing.Pool(processes=FORM_processing_options["cores"])
                graph_it = pool.imap_unordered(FORMSuperGraphIsomorphicList.generate_numerator_file_helper, 
                    list((graph, i, model, root_output_path, additional_overall_factor, workspace, integrand_type, process_definition, header_map, output_format, recycle, FORM_processing_options)
                    for i, graph in enumerate(self)))

            for (graph_index, num_ids, max_buffer_graph, is_zero, timing, code_generation_statistics) in graph_it:
                max_buffer_size = max(max_buffer_size, max_buffer_graph)
                if is_zero:
                    self[graph_index].is_zero = True
                else:
                    all_numerator_ids += num_ids

                self[graph_index][0].code_generation_statistics = code_generation_statistics

                total_time += timing
                processed_graph_counter += 1
                bar.update(processed_graph_counter, i_graph=processed_graph_counter, timing='%d'%int((total_time/float(processed_graph_counter))*1000.0))

            if FORM_processing_options["cores"] > 1:
                pool.close()

        params = copy.deepcopy(params)
        params['mUV'] = 'params[0]'
        params['mu'] = 'params[1]'
        params['gs'] = 'params[2]'
        params['small_mass_sq'] = 'params[3]'
        params['uv_cutoff_scale_sq'] = 'params[4]'

        header_code = \
"""
#ifndef NUM_H
#define NUM_H

#include <tgmath.h>

{}

#endif
""".format('\n'.join('#define  {} {}'.format(k, v) for k, v in params.items()))

        writers.CPPWriter(pjoin(root_output_path, '%(header)snumerator.h'%header_map)).write(header_code%header_map)

        numerator_C_source_size = 0.
        for fpath in glob_module.glob(pjoin(root_output_path,'*numerator*')):
            numerator_C_source_size += os.path.getsize(fpath)
        self.code_generation_statistics['numerator_C_source_size_in_kB'] = int(numerator_C_source_size/1000.0)

        if integrand_type is not None:
            self.generate_integrand_functions(root_output_path, model, additional_overall_factor=additional_overall_factor, 
                            params=params, output_format=output_format, workspace=None, integrand_type=integrand_type,process_definition=process_definition)

        self.post_process_source_code(root_output_path)
        
        generation_time = time.time() - start_time
        self.code_generation_statistics['generation_time_in_s'] = float('%.1f'%generation_time)

    def post_process_source_code(self, root_output_path):
        """ Possibly split source files if necessary and generates makefile targets."""

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield tuple(lst[i:i + n])

        # Deduce the SG ids from the files "integrand_PF_" or "integrand_LTD_"
        SG_ids = set([])
        for fpath in glob_module.glob(pjoin(root_output_path,'integrand_PF_*_f64.c')):
            match_re = re.match(r"integrand_PF_(?P<SGID>\d+)_f64.c",os.path.basename(fpath))
            if match_re is not None:
                SG_ids.add(int(match_re.group('SGID')))
        for fpath in glob_module.glob(pjoin(root_output_path,'integrand_LTD_*_f64.c')):
            match_re = re.match(r"integrand_LTD_(?P<SGID>\d+)_f64.c",os.path.basename(fpath))
            if match_re is not None:
                SG_ids.add(int(match_re.group('SGID')))

        # Now for each SG ID process files and build a target
        repl_dict = {
            'all_sg_file_names' : [],
            'all_sg_targets' : []
        }
        split_file_name_re = re.compile(r"^integrand_.*_s(\d+)_.*\.(c|cpp)$")
        dependency_re = re.compile(r"^integrand_(?P<integrand_type>[A-Za-z]*)_(?P<SG_id>\d+)(?P<split_id>_s\d+)?_(?P<precision>.*)\.o$")
        for SG_id in sorted(SG_ids):

            dependencies = []
            for code_type in ['PF','LTD',]:
                all_matches = list(glob_module.glob(pjoin(root_output_path,'integrand_%s_%d_*.c'%(code_type,SG_id))))+list(glob_module.glob(pjoin(root_output_path,'integrand_%s_%d_*.cpp'%(code_type,SG_id))))
                for fpath in all_matches:
                    # Make sure ignore split files that may have been previously generated
                    if split_file_name_re.match(os.path.basename(fpath)):
                        continue
                    dependencies.extend(self.split_source_file(root_output_path,os.path.basename(fpath),FORM_processing_options['max_n_lines_in_C_source']))
            formatted_dependencies = [d.replace('.cpp','.o').replace('.c','.o') for d in dependencies]

            repl_dict['all_sg_file_names'].append('$(LIBBPATH)/libFORM_sg_%d.so'%SG_id)
            if FORM_processing_options['max_n_files_in_library'] is None or FORM_processing_options['max_n_files_in_library']<=0 or FORM_processing_options['max_n_files_in_library'] > len(formatted_dependencies):
                repl_dict['all_sg_targets'].append(
    """$(LIBBPATH)/libFORM_sg_%(SGID)d.so: %(dependencies)s
\t$(GPP) --shared -fPIC $(CFLAGS) -o $@ $^"""%{'SGID' : SG_id, 'dependencies' : ' '.join(formatted_dependencies) })
            else:
                type_ordering = {'PF' : 0, 'LTD': 1}
                precision_ordering = {'f64': 0, 'f128': 1, 'mpfr' :2}
                parsed_dependencies = {}
                for d in formatted_dependencies:
                    d_specs = dependency_re.match(d).groupdict()
                    d_specs['split_id'] = 0 if d_specs['split_id'] is None else int(d_specs['split_id'][2:])
                    if (d_specs['integrand_type'],d_specs['precision']) in parsed_dependencies:
                        parsed_dependencies[(d_specs['integrand_type'],d_specs['precision'])].append( (d_specs['split_id'],d) )
                    else:
                        parsed_dependencies[(d_specs['integrand_type'],d_specs['precision'])] = [ (d_specs['split_id'],d), ]
                for k in parsed_dependencies.keys():
                    parsed_dependencies[k] = [ d[1] for d in sorted(parsed_dependencies[k], key=lambda v: v[0]) ]
                sorted_dependency_keys = sorted( list(parsed_dependencies.keys()), key=lambda k: ( type_ordering[k[0]], precision_ordering[k[1]] ) )

                # Collect the overall files containing the primary called functions
                primary_source_files = []
                libs_used_overall = []
                for (int_type, prec) in sorted_dependency_keys:
                    primary_source_files.append(parsed_dependencies[(int_type, prec)].pop(-1))
                    libs_used_thus_far = []
                    for i_chunk, chunk in enumerate(chunks(parsed_dependencies[(int_type, prec)],FORM_processing_options['max_n_files_in_library'])):
                        libs_used_thus_far.append('FORM_sg_%d_%s_%s_s%d'%(SG_id,int_type,prec,i_chunk+1))
                        repl_dict['all_sg_targets'].append(
    """$(LIBBPATH)/lib%(lib_name)s.so: %(dependencies)s %(lib_deps)s
\t$(GPP) --shared -fPIC $(CFLAGS) -L$(LIBBPATH) %(lib_links)s-o $@ %(dependencies)s"""%{
                            'lib_name' : libs_used_thus_far[-1], 'dependencies' : ' '.join(chunk),
                            'lib_deps' : ' '.join('$(LIBBPATH)/lib%s.so'%l for l in libs_used_thus_far[:-1]),
                            'lib_links': ' '.join('-l%s'%l for l in libs_used_thus_far[:-1])+(' 'if len(libs_used_thus_far)>1 else '')
                        })
                    libs_used_overall.extend(libs_used_thus_far)
                # Add the main library now
                repl_dict['all_sg_targets'].append(
    """$(LIBBPATH)/libFORM_sg_%(SGID)d.so: %(dependencies)s %(lib_deps)s
\t$(GPP) --shared -fPIC $(CFLAGS) -L$(LIBBPATH) %(lib_links)s-o $@ %(dependencies)s"""%{
                            'SGID' : SG_id, 'dependencies' : ' '.join(primary_source_files),
                            'lib_deps' : ' '.join('$(LIBBPATH)/lib%s.so'%l for l in libs_used_overall),
                            'lib_links': ' '.join('-l%s'%l for l in libs_used_overall)+(' 'if len(libs_used_overall)>1 else '')
                        })

        repl_dict['all_sg_file_names'] = ' '.join(repl_dict['all_sg_file_names'])
        repl_dict['all_sg_targets'] = '\n'.join(repl_dict['all_sg_targets'])

        makefile_targets_template = open(pjoin(plugin_path,'Templates','FORM_output_makefile_targets.inc'),'r').read()
        with open(pjoin(root_output_path,'makefile_targets.inc'),'w') as out:
            out.write(makefile_targets_template.format(**repl_dict))

    def split_function(self, function):

        if (FORM_processing_options['max_n_lines_in_C_function'] is None or FORM_processing_options['max_n_lines_in_C_function']<=20):# and False:
            return [function,]

        # Trim function head and trail empty lines
        while len(function)>0 and function[0].strip()=='':
            function.pop(0)
        while len(function)>0 and function[-1].strip()=='':
            function.pop(-1)

        if len(function) <= FORM_processing_options['max_n_lines_in_C_function']:
            return [function,]

        # Check if function is candidate for being split
        # Be strict intentionally here so as to easily be able to tweak what is allowed to be split and not in the C-code generation
        combined_last_line = None
        for func_line in function[::-1]:
            if func_line.strip().endswith(';'):
                if combined_last_line is not None:
                    break
                else:
                    combined_last_line = [func_line.strip(),]
            else:
                if combined_last_line is not None and func_line.strip()!='':
                    combined_last_line.insert(0,func_line.strip())
        combined_last_line = (''.join(combined_last_line)).strip()
        if not (combined_last_line.startswith('*out =') or combined_last_line.startswith('return ')):
            return [function,]

        function_prototype_re = re.compile(r'^(?P<type>\S*)\s(?P<name>[^\(]*)\((?P<args>[^\)]*)\)(?P<suffix>.*)$')
        zs_re = re.compile(r'^(?P<type>[^ ]*)\s(?P<zs>(Z\d+\_,?)*)\;$')
        z_def_re = re.compile(r'^[^ ]*Z\d+\_\s*[+|-]?=\s*')
        var_def_re = re.compile(r'^(?P<type>[^ ]*)\s(?P<name>[^ ]*)\s*=\s*')

        a_Z = re.compile(r'Z(\d+)_')

        # First combine all multiline into a single line and remove trailing newlines
        new_function_lines = []
        current_multiline_content = []
        for line in function:
            if line.endswith('\n'):
                line = line[:-1]
            if line.strip()=='':
                new_function_lines.append('')
                continue
            if not any(line.endswith(c) for c in ';}{'):
                if len(current_multiline_content)==0:
                    current_multiline_content.append(line)
                else:
                    current_multiline_content.append(line.strip())
            else:
                if len(current_multiline_content)>0:
                    current_multiline_content.append(line.strip())
                    new_function_lines.append(''.join(current_multiline_content))
                    current_multiline_content = []
                else:
                    new_function_lines.append(line)

        function = new_function_lines

        # Obtain function prototype
        function_prototype_match = function_prototype_re.match(function[0])
        if function_prototype_match is None:
            raise FormProcessingError("Parsing error when attempting to split function with prototype:\n%s"%function[0])

        mother_function_header_lines = []
        mother_function_body_lines = []
        mother_function_trail_lines = []
        function_splits = []
        args_for_splits = []
        Z_type = None
        # Now detect the "header" of the function
        in_body = False
        did_function_end = False

        def process_line(line):
            return re.sub(a_Z,r'Zs_[\1]',line)

        for line in function:
            if did_function_end:
                mother_function_trail_lines.append(process_line(line))
                continue
            if not in_body:
                if line.strip()=='':
                    mother_function_header_lines.append(line)
                    continue
                if z_def_re.match(line):
                    in_body = True
                else:
                    zs_match = zs_re.match(line)
                    if zs_match:
                        Z_type = zs_match['type'].strip()
                        max_z_num = max( int(z_i) for z_i in re.findall('Z(\d+)_',zs_match['zs']) )
                        mother_function_header_lines.append('\t%s Zs_[%d];'%(Z_type, max_z_num+1))
                        continue
                    
                    var_def_match = var_def_re.match(line)
                    if var_def_match:
                        args_for_splits.append( ( var_def_match.group('type').strip(), var_def_match.group('name') ) )
                    mother_function_header_lines.append(process_line(line))
                    continue
            
            if Z_type is None:
                raise FormProcessingError("Body of main function to split reached even though Zs definition has not been found yet.")

            if line.strip()=='':
                continue

            if z_def_re.match(line):
                function_splits.append(process_line(line))
                continue
            else:
                if not (line.strip().startswith('*out =') or line.strip().startswith('return ')):
                    raise FormProcessingError("Function body contained a line that is not a Z definition:\n'%s'"%line)
                mother_function_trail_lines.append(process_line(line))
                did_function_end = True
                continue
        
        new_functions = [ ]

        # Now time to split the function!
        # First define the prototype of the split functions and the call signature
        all_args = [ tuple(arg.strip().split(' ')) for arg in function_prototype_match.group('args').split(',')]
        
        # Remove the "out" pointer that is of no use for the split functions
        all_args = [a for a in all_args if a[1]!='out']

        # Add the necessary arguments
        all_args = [ arg for arg in args_for_splits if arg not in all_args ] + [ (Z_type, 'Zs_[]') ] + all_args
        if len(set(a[1] for a in all_args))!=len(all_args):
            raise FormProcessingError("Split function constructed contains arguments with same names but different types: %s"%str(all_args))

        split_function_prototype = 'void {}_split_%d({}){}'.format(
            function_prototype_match.group('name'),
            ', '.join( ('%s %s'%a for a in  all_args) ),
            function_prototype_match.group('suffix')
        )
        split_function_call = '\t{}_split_%d({});'.format(
            function_prototype_match.group('name'),
            ', '.join( (('%s'%a[1]).replace('[]','') for a in  all_args) )
        )

        def chunks(lst, n):
            """Yield successive n-sized chunks from lst."""
            for i in range(0, len(lst), n):
                yield lst[i:i + n]
        
        for i_chunk, chunk in enumerate(chunks(function_splits,FORM_processing_options['max_n_lines_in_C_function'])):
            new_functions.append([split_function_prototype%i_chunk,])
            mother_function_body_lines.append(split_function_call%i_chunk)

            new_functions[-1].extend(chunk)
            new_functions[-1].append('}')

        new_functions.append(mother_function_header_lines+mother_function_body_lines+mother_function_trail_lines)

        # Reinstate the new line character at the end of each line
        new_functions = [ [fl+'\n' for fl in f] for f in new_functions]
        
        return new_functions

    def split_source_file(self, root_output_path, source_code_name, max_n_lines):
        
        raw_source = open(pjoin(root_output_path, source_code_name),'r').readlines()

        start_file_re = re.compile(r'^// Source post-processed and split into the following files: (?P<dependencies>.*)$')
        start_file_match = start_file_re.match(raw_source[0].strip())
        if start_file_match is not None:
            try:
                dependencies=eval(start_file_match.group('dependencies'))
            except:
                raise FormProcessingError("Could not read dependencies from first line of file '%s'."%pjoin(root_output_path, source_code_name))
            return dependencies

        num_lines = len(raw_source)
        
        #logger.critical("I GOT THIS MANY LINES for '%s': %d (max=%s)"%(pjoin(root_output_path, source_code_name),num_lines, str(max_n_lines)))
        # First make sure that any splitting is necessary in the first place.
        # For safety impose a minimum of a thousand lines per SG otherwise it may start spewing too many files.
        if (max_n_lines is None or max_n_lines<=200 or num_lines <= max_n_lines):# and False:
            source_code_name_split = source_code_name.split('_')
            # Clear out any additional file that may have existed from a previous run
            split_file_name_template = '%s%s%s'%('_'.join(source_code_name_split[:-1]), '_s*_', source_code_name_split[-1])
            for fpath in glob_module.glob(pjoin(root_output_path,split_file_name_template)):
                os.remove(fpath)
            open(pjoin(root_output_path, source_code_name),'w').write(''.join(
                ['// Source post-processed and split into the following files: [%s,]'%source_code_name,]+raw_source
            ))
            return [source_code_name,]

        extern_c_re = re.compile(r'^extern \"C\" {$')
        start_func_re = re.compile(r"^(static|void).*{$")
        end_func_re = re.compile(r"}\s*$")
        header = []
        functions = []
        in_function = False
        in_extern_C = False
        extern_C_code = []
        was_las_match_a_closing_bracket = False
        for line in open(pjoin(root_output_path, source_code_name),'r').readlines():
            if end_func_re.match(line):
                if was_las_match_a_closing_bracket:
                    # This denotes the end of an extern C
                    extern_C_code.append(line)
                    in_extern_C = False
                    was_las_match_a_closing_bracket=False
                else:
                    was_las_match_a_closing_bracket = True
                    if in_extern_C:
                        extern_C_code.append(line)
                    else:
                        if not in_function or len(functions)==0:
                            raise FormProcessingError("Found function end before start in %s: %s"%(pjoin(root_output_path, source_code_name),line))
                        else:
                            in_function = False
                        functions[-1].append(line)
                continue

            was_las_match_a_closing_bracket=False
            if extern_c_re.match(line):
                if in_function:
                    raise FormProcessingError("Extern C function found while in the body of a function %s: %s"%(pjoin(root_output_path, source_code_name),line))
                if len(extern_C_code) > 0:
                    raise FormProcessingError("Multiple extern C blocks in %s not supported in code splitter."%pjoin(root_output_path, source_code_name))
                extern_C_code.append(line,)
                in_extern_C = True
            elif in_extern_C:
                extern_C_code.append(line)
            elif start_func_re.match(line):
                in_function = True
                functions.append([line.replace('static','').replace('inline','').strip()+'\n',])
            else:
                if not in_function:
                    if len(functions)==0:
                        header.append(line)
                    else:
                        if line.strip()=='':
                            functions[-1].append(line)
                        else:
                            raise FormProcessingError("Error: found out-of-function code that is not a header in %s: '%s'"%(pjoin(root_output_path, source_code_name),line))
                else:
                    functions[-1].append(line)

        # Possibly split functions into smaller subfunctions
        functions = sum([self.split_function(function) for function in functions],[])
        forward_declarations = [ f[0].replace('static','').replace('inline','').replace(' {',';').strip()+'\n' for f in functions]

        #print("In source code %s, found %d header lines and %d functions totalling %d lines."%(pjoin(root_output_path, source_code_name),len(header),len(functions), sum(len(f) for f in functions)))
        i_file = 0
        i_func = 0
        files_written = []
        source_code_name_split = source_code_name.split('_')
        while True:
            split_file_name = '%s%s%s'%('_'.join(source_code_name_split[:-1]), ('_s%d_'%i_file if i_file>0 else '_'), source_code_name_split[-1])
            files_written.append(split_file_name)
            with open(pjoin(root_output_path,split_file_name),'w') as out:
                out.write(''.join(header))
                if i_func>0:
                    out.write(''.join(forward_declarations[:i_func])+'\n')
                n_function_lines_accumulated = 0
                while n_function_lines_accumulated<max_n_lines and len(functions)>0:
                    i_func += 1
                    n_function_lines_accumulated += len(functions[0])
                    out.write(''.join(functions.pop(0)))
                if len(functions) == 0:
                    # Now write the extern_c code at the end
                    if len(extern_C_code)>0:
                        out.write(''.join(extern_C_code))
                    break
                else:
                    i_file += 1

        split_file_name_re = re.compile(r"^integrand_.*_s(\d+)_.*\.c$")
        # Clear out any additional file that may have existed from a previous run
        split_file_name_template = '%s%s%s'%('_'.join(source_code_name_split[:-1]), '_s*_', source_code_name_split[-1])
        for fpath in glob_module.glob(pjoin(root_output_path,split_file_name_template)):
            if int(split_file_name_re.match(os.path.basename(fpath)).group(1))>i_file:
                os.remove(fpath)

        # Flag the source file as processed and specify the list of files split into as comment of the first line
        raw_source = open(pjoin(root_output_path, source_code_name),'r').readlines()
        open(pjoin(root_output_path, source_code_name),'w').write(''.join(
            ['// Source post-processed and split into the following files: [%s%s]'%(
                ','.join(files_written), ',' if len(files_written)==1 else ''
            ),]+raw_source
        ))
        return files_written

    def generate_squared_topology_files(self, root_output_path, model, process_definition, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True, workspace=None,
        integrand_type=None, include_integration_channel_info=True, selected_cuts=None):
        if workspace is None:
            workspace = pjoin(root_output_path, os.pardir, 'workspace')
        topo_collection = {
            'name': self.name,
            'topologies': []
        }

        # Set sensible jet_ids if none
        if jet_ids is None:
            jet_ids=tuple(list(range(1,6))+list(range(-1,-6,-1))+[21,82,-82]+list(dummy_scalar_PDGs.keys()))

        contributing_supergraphs = []

        non_zero_graph = 0
        with progressbar.ProgressBar(prefix='Generating squared topology files (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}, PF #{variables.PF_config}, {variables.timing} ms / supergraph) : '%len(self), 
                max_value=len(self), variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb' : '0', 'PF_config': 'N/A'} ) as bar:
            total_time=0.0
            for i, g in enumerate(self):
                time_before = time.time()
                bar.update(i_graph='%d'%(i+1))
                if g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call=non_zero_graph, 
                                                            final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids, 
                                                            workspace=workspace, bar=bar, integrand_type=integrand_type,
                                                            include_integration_channel_info=include_integration_channel_info, 
                                                            selected_cuts=selected_cuts):
                    topo_collection['topologies'].append({
                        'name': g[0].name,
                        # Let us not put it there but in the topology itself
                        #'benchmark_result': g[0].benchmark_result,
                        'multiplicity': g[0].multiplicity
                        ,'additional_LMBs': [
                            {
                                'name' : '%s_LMB%d'%(other_supergraph.name,other_supergraph.additional_lmbs),
                                'defining_lmb_to_this_lmb' : original_lmb_to_other_lmb_affine_transfo,
                                'this_lmb_to_defining_lmb' : other_lmb_to_original_lmb_affine_transfo
                            }
                            for i_lmb, other_lmb_to_original_lmb_affine_transfo, 
                                original_lmb_to_other_lmb_affine_transfo, other_supergraph in g[0].additional_lmbs
                        ]

                    })
                    non_zero_graph += 1
                    contributing_supergraphs.append(g)

                total_time += time.time()-time_before
                bar.update(timing='%d'%int((total_time/float(i+1))*1000.0))
                bar.update(i+1)
        
        if filter_non_contributing_graphs:
            self[:] = contributing_supergraphs

        if FORM_processing_options['generate_renormalisation_graphs']:
            renormalization_graphs = self.generate_renormalization_graphs(model,process_definition)
        else:
            renormalization_graphs = []
        self.extend([FORMSuperGraphIsomorphicList([g]) for g in renormalization_graphs])

        with progressbar.ProgressBar(prefix='Generating renormalization squared topology files (graph #{variables.i_graph}/%d, LMB #{variables.i_lmb}/{variables.max_lmb}, PF #{variables.PF_config}, {variables.timing} ms / supergraph) : '%len(renormalization_graphs), 
                max_value=len(renormalization_graphs), variables = {'timing' : '0', 'i_graph' : '0', 'i_lmb': '0', 'max_lmb' : '0', 'PF_config': 'N/A'} ) as bar:
            total_time=0.0
            for i, g in enumerate(renormalization_graphs):
                time_before = time.time()
                bar.update(i_graph='%d'%(i+1))
                if g.generate_squared_topology_files(root_output_path, model, process_definition, n_jets, numerator_call=non_zero_graph, 
                                                            final_state_particle_ids=final_state_particle_ids,jet_ids=jet_ids,
                                                            workspace=workspace, bar=bar, integrand_type=integrand_type,
                                                            include_integration_channel_info=include_integration_channel_info,
                                                            selected_cuts=selected_cuts):
                    topo_collection['topologies'].append({
                        'name': g.name,
                        'multiplicity': g.multiplicity,
                        'additional_LMBs': []
                    })
                    non_zero_graph += 1
                    contributing_supergraphs.append(g)

                total_time += time.time()-time_before
                bar.update(timing='%d'%int((total_time/float(i+1))*1000.0))
                bar.update(i+1)

        open(pjoin(root_output_path, self.name + '.yaml'), 'w').write(yaml.dump(topo_collection, Dumper=Dumper))

    def get_renormalization_vertex(self, in_pdgs, loop_count, model, process_definition, is_external_bubble=False):
        
        is_heft = (model.get('name')=='SM_HEFT')

        # This overall factor is an overall fudge which will need to be understood together with how exactly we must
        # inherit the overall phase of the reference loop SG that this vertex is supposed to renormalised. 
        overall_factor = '(-1)'

        if loop_count != 1:
            raise FormProcessingError("Renormalisation is currently implemented only at one-loop.")
            
        if any( model.get_particle(pdg).get('ghost') for pdg in in_pdgs):
            # For now we absorbed all renormalisation constants into the gluon contribution so that we set the ghost 
            # renormalisation conditions to zero. This implies that the local UV pole cancellation check must be done
            # by combining the ghost and gluon supergraph with appropriate routing.
            return '0'

        # WARNING: All the things below are a SHAME and must be generalised as soon as we're done with the publication....

        # All renormalisation 2-point vertices are for now coded up directly in numerator.frm
        if len(in_pdgs)==2: #and not is_external_bubble:
            return overall_factor

        # No renormalisation for shrunk loop-induced vertices involving only gluons and the Higgs
        if not is_heft and all(pdg in [21,25] for pdg in in_pdgs):
            return '0'

        hardcoded_mass_parameters = {
            6   : 'masst',
            -6  : 'masst',
        }
        
        hardcoded_log_quark_mass = {
            6   : 'logmt',
            -6  : 'logmt'
        }

        def get_particle_mass(pdg):
            if pdg in hardcoded_mass_parameters:
                return hardcoded_mass_parameters[pdg]
            return model.get_particle(pdg).get('mass')

        pdgs = sorted([abs(pdg) for pdg in in_pdgs],reverse=True)
        quark_pdgs = range(1,7)

        active_quarks = [q_pdg for q_pdg in quark_pdgs if q_pdg not in process_definition.get('forbidden_particles')]
        n_massless_quarks = len([1 for q_pdg in active_quarks if get_particle_mass(q_pdg).upper()=='ZERO' ])
        massive_quark_pdgs = [ q_pdg for q_pdg in active_quarks if get_particle_mass(q_pdg).upper()!='ZERO' ]

        symbol_replacement_dict = {
            'C_F' : '(4/3)',
            'T_F' : '(1/2)',
            'C_A' : '3',
            'pi'  : 'pi',
            'gs'  : 'gs',
            'ep'  : 'ep',
            'mu_r': 'mu_r',
            'n_f' : '%d'%n_massless_quarks,
            'UVRenormFINITE_marker': 'UVRenormFINITE'
        }
        symbol_replacement_dict['beta0'] = '(11-(2/3)*%(n_f)s)'%symbol_replacement_dict

        delta_Z_massless_quark = '%(C_F)s*%(gs)s^2/16/%(pi)s^2*(-1/%(ep)s)'
        delta_Z_massive_quark = '%(C_F)s*%(gs)s^2/16/%(pi)s^2*((-1/%(ep)s) + %(UVRenormFINITE_marker)s*(-4-3*(logmu - %(log_quark_mass)s)) )'

        # Note the many overall factors (1/2) come from the fact that the expressions are derived from the alpha couplings (i.e. g^2/4\pi)
        # and not for the coupling g directly, so that a factor 1/2 is warranted.

        delta_Z_gluon = []
        # Add the gluon and ghost contribution to the gluon wavefunction renormalisation
        delta_Z_gluon.append('( (1/2)*%(C_A)s*%(gs)s^2/48/%(pi)s^2 * ( 5/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massless quark contributions to the gluon wavefunction renormalisation
        delta_Z_gluon.append('( (1/2)*%(n_f)s*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( -4/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massive quark contribution to the gluon wavefunction renormalisation
        for q_pdg in massive_quark_pdgs:

            if q_pdg not in hardcoded_log_quark_mass:
                raise FormProcessingError( 
                    ("There is not hard-coded symbol for the logarithm of the mass of particle with PDG=%d.\n"%q_pdg)+
                    "Are you sure you did not to specify forbidden particles in your process definition, using the / x y z syntax?")

            delta_Z_gluon.append('( (1/2)*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( -4/%(ep)s + %(UVRenormFINITE_marker)s*(-4*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the gluon wavefunction renormalisation
        delta_Z_gluon = '(%s)'%('+'.join(delta_Z_gluon))

        delta_g_s = []
        # Add the gluon and ghost contribution to the strong coupling renormalisation
        delta_g_s.append('( (1/2)*%(C_A)s*%(gs)s^2/48/%(pi)s^2 * ( -11/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massless quark contributions to the strong coupling renormalisation
        delta_g_s.append('( (1/2)*%(n_f)s*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( 4/%(ep)s ) )'%symbol_replacement_dict)
        # Add the massive quark contribution to the strong coupling renormalisation
        for q_pdg in massive_quark_pdgs:
            delta_g_s.append('( (1/2)*%(T_F)s*%(gs)s^2/48/%(pi)s^2 * ( 4/%(ep)s + %(UVRenormFINITE_marker)s*(4*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the strong coupling renormalisation
        delta_g_s = '(%s)'%('+'.join(delta_g_s))

        delta_g_heft = '( (%(gs)s^2/4/%(pi)s^2)*( -(1/%(ep)s)*%(beta0)s + 11/4 ) )'%symbol_replacement_dict

        # Combine all terms for the yukawa renormalisation
        # WARNING: He consider here only contribution from the top and bottom as these are the only fermions we consider as
        # having a yukawa interaction in the model, but we do allow for them to be massless and adjust the renormalisation accordingly.
        potential_yukawa_quark_pdgs = [5,6]
        active_yukawa_quarks = [q_pdg for q_pdg in potential_yukawa_quark_pdgs if q_pdg not in process_definition.get('forbidden_particles')]
        n_massless_yukawa_quarks = len([1 for q_pdg in active_yukawa_quarks if get_particle_mass(q_pdg).upper()=='ZERO' ])
        massive_yukawa_quark_pdgs = [ q_pdg for q_pdg in active_yukawa_quarks if get_particle_mass(q_pdg).upper()!='ZERO' ]

        delta_yukawa = []
        # Add the massless quark contributions to the yukawa renormalisation
        delta_yukawa.append('( -(1/2)*%(n_f_yukawa)s*2*%(C_F)s*%(gs)s^2/16/%(pi)s^2 * ( 3/%(ep)s ) )'%(
            dict(symbol_replacement_dict,**{
                'n_f_yukawa':'%d'%n_massless_yukawa_quarks,
            })
        ))
        # Add the massive quark contribution to the yukawa renormalisation
        for q_pdg in massive_yukawa_quark_pdgs:
            delta_yukawa.append('( - (1/2)*2*%(C_F)s*%(gs)s^2/16/%(pi)s^2 * ( 3/%(ep)s + %(UVRenormFINITE_marker)s*(4 + 3*(logmu - %(log_quark_mass)s)) ) )'%(
                dict(symbol_replacement_dict,**{
                    'log_quark_mass':hardcoded_log_quark_mass[q_pdg],
                })
            ))
        # Combine all terms for the yukawa renormalisation
        delta_yukawa = '(%s)'%('+'.join(delta_yukawa))

        ########################################################
        # Now build the specified UV renormalisation vertices
        ########################################################

#       Self-energies are for now coded up direclty in numerator.frm

        # quark self-energy
#        if len(pdgs) == 2 and pdgs[0] in quark_pdgs and pdgs[1]==pdgs[0]:
#            quark_mass = get_particle_mass(pdgs[0])
#            if quark_mass=='ZERO':
#                res = ('(-(%s))'%delta_Z_massless_quark)%symbol_replacement_dict
#                res = res
#                return res
#            else:
#                res = ('(-(%s))'%delta_Z_massive_quark)%(
#                    dict(symbol_replacement_dict,**{
#                        'quark_mass':quark_mass,
#                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[0]],                            
#                    }))
#                return res

        # gluon self-energy
#        if len(pdgs) == 2 and pdgs[0]==21 and pdgs[1]==pdgs[0]:
#            return '(-(%s))'%delta_Z_gluon

        # a/z qq vertex
        if len(pdgs) == 3 and (pdgs[0] in [22,23]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 1
            quark_mass = get_particle_mass(pdgs[1])
            res = '%s*'%overall_factor
            if quark_mass=='ZERO':
                res += ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict
            else:
                res += ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    }))
            return '((1/%d)*%s)'%(loop_multiplicity_factor,res)

        # hqq vertex
        if len(pdgs) == 3 and (pdgs[0] in [25,]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 1
            res = [ '(+%s)'%delta_yukawa, ]
            quark_mass = get_particle_mass(pdgs[1])
            if quark_mass=='ZERO':
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict )
            else:
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    })) )
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor,overall_factor,('+'.join(res)))

        # gqq vertex
        if len(pdgs) == 3 and (pdgs[0] in [21,]) and (pdgs[1] in quark_pdgs) and (pdgs[2]==pdgs[1]):
            loop_multiplicity_factor = 2
            res = [ 
                '(+%s)'%delta_g_s, 
                '(+%s)'%delta_Z_gluon
            ]
            quark_mass = get_particle_mass(pdgs[1])
            if quark_mass=='ZERO':
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massless_quark)%symbol_replacement_dict )
            else:
                res.append( ('(+(2)*(1/2)*(%s))'%delta_Z_massive_quark)%(
                    dict(symbol_replacement_dict,**{
                        'quark_mass':quark_mass,
                        'log_quark_mass':hardcoded_log_quark_mass[pdgs[1]],
                    })) )
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor, overall_factor,('+'.join(res)))

        # ggg vertex
        if len(pdgs) == 3 and all(pdg==21 for pdg in pdgs):
            # 2 per n_f, 2 per ghost and 1 for the gluon
            loop_multiplicity_factor = 1 + 2 + 2 * n_massless_quarks
            res = [ 
                '(+%s)'%delta_g_s,
                '(+3*(%s))'%delta_Z_gluon,
            ]
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor, overall_factor,('+'.join(res)))

        # gggg vertex
        if len(pdgs) == 4 and all(pdg==21 for pdg in pdgs):
            # 6 per n_f, 6 per ghost and 3 for the gluon
            loop_multiplicity_factor = 3 + 6 + 6 * n_massless_quarks
            res = [ 
                '(+2*(%s))'%delta_g_s,
                '(+4*(%s))'%delta_Z_gluon,
            ]
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor,overall_factor,('+'.join(res)))

        # hgg, hggg and hgggg vertex
        if all(pdg in [21,25] for pdg in in_pdgs) and in_pdgs.count(25)==1 and len(in_pdgs)<=5:
            if len(in_pdgs)==3:
                loop_multiplicity_factor = 4
            elif len(in_pdgs)==4:
                loop_multiplicity_factor = 15
            elif len(in_pdgs)==5:
                loop_multiplicity_factor = 69
            else:
                raise FormProcessingError("Unreachable")
            res = [ 
                '(+%s)'%delta_g_heft, 
                '(+%d*(%s))'%(len(in_pdgs)-1,delta_Z_gluon)
            ]
            return '((1/%d)*%s*(%s))'%(loop_multiplicity_factor, overall_factor,('+'.join(res)))

        return None

    @classmethod
    def shrink_edges(cls, edges, nodes, edges_to_shrink, bubble_external_shrinking_step=False):
        subgraph_pdgs = set()
        for ek in list(edges.keys()):
            ee = edges[ek]
            if ee['name'] not in edges_to_shrink:
                continue
            subgraph_pdgs.add(ee['PDG'])

            # remove the edge from the vertices
            for na in ee['vertices']:
                node = nodes[na]
                edge_index = node['edge_ids'].index(ek)
                for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                    node[g] = tuple(ind for i, ind in enumerate(node[g]) if i != edge_index)

            # fuse vertices
            if ee['vertices'][0] != ee['vertices'][1]:
                node = nodes[ee['vertices'][0]]
                # START TEMPORARY HACK (Just for the node rendering purposes)
                if not bubble_external_shrinking_step:
                    # Mark this node as being shrunk. Eventually there will be a more general way of handling
                    # multiple effective shrunk renormalisation nodes. For now, I just flag it here as one-loop
                    node['renormalisation_vertex_n_loops'] = 1
                    node['renormalisation_vertex_n_shrunk_edges'] = len(edges_to_shrink)
                # END TEMPORARY HACK
                for na in ee['vertices'][1:]:
                    n = nodes[na]
                    for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                        node[g] = tuple(list(node[g]) + list(n[g]))
                    # START TEMPORARY HACK (Just for the node rendering purposes)
                    if bubble_external_shrinking_step:
                        for g in ['renormalisation_vertex_n_loops', 'renormalisation_vertex_n_shrunk_edges']:
                            if g not in node and g in n:
                                node[g] = n[g]
                    # END TEMPORARY HACK
                    del nodes[na]

                vert_to_replace = ee['vertices'][1]
                vert_to_replace_with = ee['vertices'][0]
                for e in edges.values():
                    e['vertices'] = tuple(vert_to_replace_with if ind == vert_to_replace else ind for ind in e['vertices'])

            del edges[ek]

        return subgraph_pdgs

    def generate_renormalization_graphs(self, model, process_definition):
        """Generate all required renormalization graphs for the graphs in the supergraph list"""

        import sympy as sp
        n_externals = max(len([1 for e in graph[0].edges.values() if e['type']=='in' or e['type']=='out']) for graph in self)
        if model is None:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([1,2,3,4,5,6,11,12,13,21,22,25,82])}
        else:
            pdg_primes = {pdg : sp.prime(i + n_externals + 1) for i, pdg in enumerate([p['pdg_code'] for p in model['particles']])}

        renormalization_graphs = []
        for i_sg, gs in enumerate(self):
            squared_topology = gs[0].squared_topology
            for i_cut, cut_info in enumerate(squared_topology.cuts):
                for i_diag_set, diag_set in enumerate(cut_info['diagram_sets']):
                    if len(diag_set['uv_propagators']) == 0:
                        continue

                    # create the renormalization graph per diagram set
                    nodes = copy.deepcopy(gs[0].nodes)
                    edges = copy.deepcopy(gs[0].edges)
                    vertex_factors = []
                    uv_effective_vertex_edges = []

                    # BEGIN PIECE OF MULTIPLICITY HACK
                    has_pure_external_gluon_effective_vertex = False
                    # END PIECE OF MULTIPLICITY HACK

                    # Compute the bubble edges so as to be able to determine if the bubble to computed the renormalisation vertex for is "internal" or "external".
                    bubble_edges = set()
                    for c in cut_info['cuts']:
                        sig = next(ee['signature'] for ee in edges.values() if ee['name'] == c['edge'])
                        inv_sig = tuple([-s for s in x] for x in sig)
                        bubble_edges |= set(ee['name'] for ee in edges.values() if (ee['signature'] == sig or ee['signature'] == inv_sig) and ee['name'] != c['edge'])

                    for diag_info in diag_set['diagram_info']:
                        if diag_info['uv_vertices'] is not None and len(diag_info['uv_vertices']) > 0:
                            graph = diag_info['graph']
                            # look up every uv vertex in the table
                            for v, l in diag_info['uv_vertices']:
                                edges_on_vertex = [e[0] for e in graph.edge_map_lin if v in e[1:]]
                                uv_effective_vertex_edges.append((edges_on_vertex, l))
                                edge_pdgs = tuple(next(ee for ee in edges.values() if ee['name'] == e)['PDG'] for e in edges_on_vertex)
                                # BEGIN PIECE OF MULTIPLICITY HACK
                                if not has_pure_external_gluon_effective_vertex:
                                    has_pure_external_gluon_effective_vertex = set(edge_pdgs)==set([21,])
                                # END PIECE OF MULTIPLICITY HACK
                                is_external_bubble = ( (len(edges_on_vertex)==2) and any( (e in bubble_edges) for e in edges_on_vertex) )
                                #misc.sprint(edge_pdgs)
                                #misc.sprint(l)
                                #misc.sprint(process_definition.nice_string())
                                #misc.sprint(is_external_bubble)
                                vertex_contrib = self.get_renormalization_vertex(edge_pdgs, l, model, process_definition, is_external_bubble=is_external_bubble)
                                #misc.sprint(vertex_contrib)

                                if vertex_contrib is None:
                                    logger.warning("WARNING: unknown renormalization vertex {} at {} loops".format(edge_pdgs, l))
                                    #raise AssertionError("Unknown renormalization vertex {} at {} loops".format(edge_pdgs, l))
                                    # It can often happen, for example for (but not only) Loop-Induced processes, that there must be
                                    # not renormalisation vertex associated to an integrated UV CT. For instance g g H H.
                                    # Ideally one should code those explicitly in the function `get_renormalization_vertex`, but for 
                                    # now it is also OK to simply let it return None with the above warning and use `0` as the numerator
                                    # for this renormalisation supergraph, so that it hopefully gets removed automatically at the time 
                                    # of explicitly building its numerator.
                                    vertex_factors.append('(0)')
                                else:
                                    vertex_factors.append('(' + vertex_contrib + ')')

                    # If this cut and diag set receives no contribution, then ignore it
                    if all(v=='(0)' for v in vertex_factors):
                        continue

                    # shrink the UV subgraph in the edge list
                    subgraph_pdgs = self.shrink_edges(edges, nodes, diag_set['uv_propagators'],bubble_external_shrinking_step=False)

                    # BEGIN PIECE OF MULTIPLICITY HACK
                    is_made_of_gluons_only = set(subgraph_pdgs)==set([21,])
                    # END PIECE OF MULTIPLICITY HACK

                    # Recomputed the bubble edges now the that UV subgraph has been shrunk, they will be removed at the end
                    bubble_edges = set()
                    for c in cut_info['cuts']:
                        sig = next(ee['signature'] for ee in edges.values() if ee['name'] == c['edge'])
                        inv_sig = tuple([-s for s in x] for x in sig)
                        bubble_edges |= set(ee['name'] for ee in edges.values() if (ee['signature'] == sig or ee['signature'] == inv_sig) and ee['name'] != c['edge'])

                    # compute the proper multiplicity by dividing out the symmetry factor of the UV divergent subgraphs
                    # TODO: make multi-loop and multi-uv subgraph compatible
                    # In the new paradigm this is done automatically as it is divided out in the corresponding renormalisation counterterms.
                    multiplicity = gs[0].multiplicity
#                    if len(set([21, 22]) & subgraph_pdgs) == 0:
#                        multiplicity /= 2

                    # set the correct particle ordering for all the edges
                    for n in nodes.values():
                        if len(n['PDGs']) > 1:
                            edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                            for g in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                                n[g] = tuple(n[g][eo['index']] for eo in edge_order)

                    # check if the remaining graph is seen before
                    # NOTE: this must be exactly the same procedure as the graph isomorphic filtering,
                    # to obtain the correct multiplicity
                    # TODO: The multiplicity factor must be computed accordingly to what was already done
                    # at LO for the original loop SG, that is including the color factors as well.
                    g = igraph.Graph()
                    vertex_map = list(nodes.keys())
                    g.add_vertices(len(vertex_map))
                    undirected_edges = set()

                    for e in edges.values():
                        undirected_edges.add(tuple(sorted(vertex_map.index(v) for v in e['vertices'])))

                    edge_colors = []
                    ext_id = 0
                    for ue in undirected_edges:
                        e_color = 1
                        for e in edges.values():
                            if tuple(sorted(vertex_map.index(v) for v in e['vertices'])) == ue:
                                # find the associated edge to get the PDG
                                edge = next(ee for ee in edges.values() if ee['name'] == e['name'])
                                if edge['type'] == 'in':
                                    ext_id += 1
                                    e_color *= sp.prime(ext_id)
                                else:
                                    e_color *= pdg_primes[abs(edge['PDG'])]
                        edge_colors.append(e_color)
                        g.add_edges([ue])

                    v_colors = [0]*len(vertex_map)
                    for v, l in uv_effective_vertex_edges:
                        # find the vertex that is shrunk and give it a color
                        n = next(ni for ni, n in nodes.items() if set(v) == set(edges[e]['name'] for e in n['edge_ids']))
                        v_colors[vertex_map.index(n)] = l

                    for (ref_graph, ref_v_colors, ref_e_colors), ref_ren_graph in renormalization_graphs:

                        # BEGIN PIECE OF MULTIPLICITY HACK
                        # If this loop super graph is a multi-gluon one but does not contain only gluons in its loop content
                        # then we skip it here immediately as we want to instead use the loop graphs with gluon only as a 
                        # representative so that we can inherit the right symmetry factor
#                        if has_pure_external_gluon_effective_vertex and not is_made_of_gluons_only:
#                            break
                        # END PIECE OF MULTIPLICITY HACK

                        if ref_graph.get_isomorphisms_vf2(g,
                            color1=ref_v_colors,
                            color2=v_colors,
                            edge_color1=ref_e_colors,
                            edge_color2=edge_colors) != []:

                            # Accumulate the multiplicity in the renormalisation SG since the right quantity
                            # will be divided out in the definition of the renormalisation vertex
                            # However we must not aggregate diagram sets from the same supergraph and same uv propgators
                            if not any(
                                (   ( self[o_sg][0].name == gs[0].name ) and 
                                        ( set(self[o_sg][0].squared_topology.cuts[o_cut]['diagram_sets'][o_diag_set]['uv_propagators']) ==
                                          set(diag_set['uv_propagators']) )
                                )  for o_sg,o_cut,o_diag_set in ref_ren_graph.matching_loop_subdiags):
                                    ref_ren_graph.multiplicity += multiplicity
                                    ref_ren_graph.matching_loop_subdiags.append((i_sg,i_cut,i_diag_set))
                                    ref_ren_graph.name += '_%sDS%d'%(gs[0].name,diag_set['id'])

                            # BEGIN PIECE OF MULTIPLICITY HACK
#                            if not has_pure_external_gluon_effective_vertex:
#                                if multiplicity != ref_ren_graph.multiplicity:
#                                    raise FormProcessingError(
#                                        "Two loop SG mapped to the same renormalisation graph have "+
#                                        "different multiplicities even though it is expected that they should not.")
                            # END PIECE OF MULTIPLICITY HACK
                            
                            break
                    else:
                        # remove all bubble edges, since they will cancel with the effective vertex
                        #subgraph_pdgs = self.shrink_edges(edges, nodes, bubble_edges, bubble_external_shrinking_step=True)
                        
                        # set the correct particle ordering for all the edges
                        for n in nodes.values():
                            if len(n['PDGs']) > 1:
                                edge_order = FORMSuperGraph.sort_edges(model, [{'PDG': pdg, 'index': i} for i, pdg in enumerate(n['PDGs'])])
                                for n_keys in ('PDGs', 'indices', 'momenta', 'edge_ids'):
                                    n[n_keys] = tuple(n[n_keys][eo['index']] for eo in edge_order)

                        # WARNING: Double-check that it is indeed always the right thing to do to inherit the overall factor
                        # in the renormalisation graph from the reference one.

                        # add a new supergraph for this renormalization component
                        form_graph = FORMSuperGraph(name='renorm_{}DS{}'.format(gs[0].name,diag_set['id']),
                                edges=edges, nodes=nodes,
                                overall_factor='(%s)*(%s)'%(
                                    gs[0].overall_factor, ('*'.join(vertex_factors))
                                ), multiplicity=multiplicity)

                        # Add a container with all SG names mapped to this renormalisation graph
                        form_graph.matching_loop_subdiags = [(i_sg,i_cut,i_diag_set)]

                        # set a basis
                        topo_generator = LTD.ltd_utils.TopologyGenerator([(e['name'], e['vertices'][0], e['vertices'][1]) for e in edges.values()])
                        topo_generator.generate_momentum_flow()
                        sig = topo_generator.get_signature_map()

                        for eid, e in edges.items():
                            e['signature'] = [sig[e['name']][0],
                                [ i+o for i,o in zip(sig[e['name']] [1][:len(sig[e['name']] [1])//2],sig[e['name']] [1][len(sig[e['name']] [1])//2:]) ]]
                            e['momentum'] = FORMSuperGraph.momenta_decomposition_to_string(e['signature'], set_outgoing_equal_to_incoming=False)

                            if len(e['vertices']) == 1:
                                continue

                            for i, vi in enumerate(e['vertices']):
                                e_index = nodes[vi]['edge_ids'].index(eid)
                                mom = list(nodes[vi]['momenta'])
                                mom[e_index] = '-({})'.format(e['momentum']) if i == 0 else e['momentum']
                                nodes[vi]['momenta'] = mom

                        renormalization_graphs.append(((g, v_colors, edge_colors), form_graph))

        return [g for _, g in renormalization_graphs]

    def produce_output(self):
        form_processor = FORMProcessor(self, computed_model, process_definition)
        TMP_OUTPUT = pjoin(root_path, 'TEST_' + self.name.split('_')[0])
        Path(TMP_OUTPUT).mkdir(parents=True, exist_ok=True)
        TMP_workspace = pjoin(TMP_OUTPUT, 'FORM', 'workspace')
        Path(TMP_workspace).mkdir(parents=True, exist_ok=True)
        TMP_FORM = pjoin(TMP_OUTPUT, 'FORM')
        #FORM_processing_options['compilation-options'] += ['-e', "OPTIMIZATION_LVL=2"]
        form_processor.generate_squared_topology_files(TMP_OUTPUT, 0,
                    workspace=TMP_workspace,
                    integrand_type='PF',
                    include_integration_channel_info=True,
                    selected_cuts=None)
        form_processor.generate_numerator_functions(TMP_FORM, output_format='c',
                    workspace=TMP_workspace,
                    integrand_type='PF')
        # copy the makefile, stripping the first line
        source_file = open('Templates/FORM_output_makefile', 'r')
        source_file.readline()
        with open(pjoin(TMP_FORM, 'Makefile'), 'w') as target_file:
            shutil.copyfileobj(source_file, target_file)
        if os.path.isfile(pjoin('Templates','makefile_user_opts.inc')):
            shutil.copy(pjoin('Templates','makefile_user_opts.inc'), pjoin(TMP_FORM, 'makefile_user_opts.inc'))
        else:
            shutil.copy(pjoin('Templates','makefile_user_opts_default.inc'), pjoin(TMP_FORM, 'makefile_user_opts.inc'))

        for n in ('mpcomplex.h', 'mpreal.h', 'dual.h', 'dualt2.h', 'dualkt2.h', 'dualt3.h', 'dualkt3.h', 'dualklt3.h'):
            shutil.copy(n, pjoin(TMP_FORM, n))
        Path(pjoin(TMP_OUTPUT, 'lib')).mkdir(parents=True, exist_ok=True)
        FORMProcessor.compile(TMP_FORM)

        drawings_output_path = pjoin(TMP_OUTPUT, 'Drawings')
        Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                    pjoin(drawings_output_path,'Makefile'))
        form_processor.draw(drawings_output_path)


class FORMProcessor(object):
    """ A class for taking care of the processing of a list of FORMSuperGraphList.
    Useful because many aspects common to all supergraphs and function do not belong to FORMSuperGraphList.
    """

    forced_options = None

    def __init__(self, super_graphs_list, model, process_definition, force_options=None):
        """ Specify aditional information such as the model that is useful for FORM processing."""
        self.super_graphs_list = super_graphs_list
        self.model = model
        self.process_definition = process_definition
        self.forced_options = force_options

        if isinstance(self.process_definition, base_objects.ProcessDefinition):
            all_processes = list(proc for proc in self.process_definition)
            self.repr_process = all_processes[0]
        else:
            self.repr_process = self.process_definition

        if not FORM_processing_options['generate_integrated_UV_CTs']:
            logger.warning('%s\n\nGeneration of integrated UV CTs is disabled per user request. Physical results will be incorrect.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if not FORM_processing_options['generate_renormalisation_graphs']:
            logger.warning('%s\n\nGeneration of renormalisation contributions is disabled per user request. Physical results will be incorrect.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['UV_min_dod_to_subtract']>0:
            logger.warning('%s\n\nAs per user request, not all UV divegences will be locally subtracted. Numerical integration may be divergent in the UV.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['UV_min_dod_to_subtract']<0:
            logger.warning('%s\n\nAs per user request, subleading UV divegences will also be locally subtracted. Very complicated integrands may result from this.\n\n%s'%(utils.bcolors.RED,utils.bcolors.ENDC))
        if FORM_processing_options['selected_epsilon_UV_order']!=0:
            logger.warning('%s\n\nAs per user request, the selected epsilon order to be exported will be %d. This must be for pole cancellation check only. \n\n%s'%(
                                                        utils.bcolors.RED, FORM_processing_options['selected_epsilon_UV_order'], utils.bcolors.ENDC))
        if FORM_processing_options['renormalisation_finite_terms']!='together':
            logger.warning(("%s\n\nAs per user request, the finite part of the renormalisation counterterms have been elected to be included "+
                            "as '%s' (and not the default 'together'). Results will be incorrect if not post-processed.\n\n%s")%(
                utils.bcolors.RED, FORM_processing_options['renormalisation_finite_terms'], utils.bcolors.ENDC))

    def draw(self, output_dir):
        """ For now simply one Mathematica script per supergraph."""
        if self.forced_options:
            FORM_processing_options=self.forced_options

        for i_graph, super_graphs in enumerate(self.super_graphs_list):
            super_graphs[0].draw(self.model, output_dir, FORM_id=i_graph)
            # Draw supergraphs for additional LMBs
            if isinstance(super_graphs[0].additional_lmbs,list):
                for i_lmb,_,_,sg in super_graphs[0].additional_lmbs:
                    sg.draw(self.model, output_dir, FORM_id=i_graph, lmb_id=i_lmb)

    def generate_numerator_functions(self, root_output_path, output_format='c',workspace=None, header="", integrand_type=None, force_overall_factor=None, recycle=False, additional_params=None):
        assert(header in ['MG', 'QG', ''])

        if self.forced_options:
            FORM_processing_options=self.forced_options

        params = {
            'masst': self.model['parameter_dict'][self.model.get_particle(6).get('mass')].real,
            'massb': self.model['parameter_dict'][self.model.get_particle(5).get('mass')].real,
            'massh': self.model['parameter_dict'][self.model.get_particle(25).get('mass')].real,
            'gs': self.model['parameter_dict']['G'].real,
            'ge': math.sqrt(4. * math.pi / self.model['parameter_dict']['aEWM1'].real),
            'yukawat': self.model['parameter_dict']['mdl_yt'].real / math.sqrt(2.),
            'yukawab': self.model['parameter_dict']['mdl_yb'].real / math.sqrt(2.),
            'ghhh': 6. * self.model['parameter_dict']['mdl_lam'].real,
            'vev': self.model['parameter_dict']['mdl_vev'].real,
            'pi': 'M_PI',
        }
        if additional_params is not None:
            params.update(additional_params)

        if force_overall_factor is None:
            helicity_averaging_factor = 1
            for leg in self.repr_process.get('legs'):
                # Skip final states
                if leg.get('state') is True:
                    continue

                helicity_averaging_factor *= len(self.model.get_particle(leg.get('id')).get_helicity_states())
            helicity_averaging_factor = "/" + str(helicity_averaging_factor)
            additional_overall_factor = helicity_averaging_factor
        else:
            additional_overall_factor = '*(%s)'%force_overall_factor

        res = self.super_graphs_list.generate_numerator_functions(
            root_output_path,
            model=self.model,
            output_format=output_format,
            additional_overall_factor=additional_overall_factor,
            params=params,workspace=workspace, header=header,
            integrand_type=integrand_type,
            process_definition=self.process_definition,
            recycle=recycle
        )

        self.super_graphs_list.aggregate_code_generation_statistics()
        self.report_generation_statistics(root_output_path=root_output_path)
    
        return res

    def report_generation_statistics(self,root_output_path=None):

        # TODO improve formatting
        logger.info("Generation statistics:\n%s%s%s"%(
            utils.bcolors.GREEN,
            pformat(self.super_graphs_list.code_generation_statistics),
            utils.bcolors.ENDC
        ))
        
        if root_output_path is not None:
            open(pjoin(root_output_path,'generation_statistics.txt'),'w').write(pformat(self.super_graphs_list.code_generation_statistics))

    @classmethod
    def compile(cls, root_output_path, arg=[]):

        compile_options = cls.forced_options if cls.forced_options is not None else FORM_processing_options

        if os.path.isfile(pjoin(root_output_path,'Makefile')):
            try:
                logger.info("Now compiling FORM-generated numerators with options: %s ..."%(' '.join(compile_options['compilation-options'])))

                t = time.time()
                misc.compile(arg=compile_options['compilation-options'] ,cwd=root_output_path,mode='cpp', nb_core=compile_options["cores"])
                logger.info("Compilation time: {:.2}s".format(time.time() - t))
            except MadGraph5Error as e:
                logger.info("%sCompilation of FORM-generated numerator failed:\n%s%s"%(
                    utils.bcolors.RED,str(e),utils.bcolors.ENDC))
        else:
            logger.warning(("\n%sYou are running FORM_processing directly from the __main__ of FORM_processing.py.\n"+
                           "You will thus need to compile numerators.c manually.%s")%(utils.bcolors.GREEN, utils.bcolors.ENDC))

    def generate_squared_topology_files(self, root_output_path, n_jets, final_state_particle_ids=(), jet_ids=None, filter_non_contributing_graphs=True, workspace=None,
        integrand_type=None, include_integration_channel_info=None, selected_cuts=None):

        if self.forced_options:
            FORM_processing_options=self.forced_options

        self.super_graphs_list.generate_squared_topology_files(
            root_output_path, self.model, self.process_definition, n_jets, final_state_particle_ids, jet_ids=jet_ids, filter_non_contributing_graphs=filter_non_contributing_graphs, workspace=workspace,
            integrand_type=integrand_type, include_integration_channel_info=(FORM_processing_options['include_integration_channel_info'] if include_integration_channel_info is None else include_integration_channel_info),
            selected_cuts = selected_cuts
        )

if __name__ == "__main__":
   
    parser = argparse.ArgumentParser(description='Generate numerators with FORM and yaml for Rust.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--diagrams_python_source', default='None', type=str,
                        help='path to the python diagram output files.')
    parser.add_argument('--model', default='sm-no_widths', type=str,
                        help='Path to UFO model to load.')
    parser.add_argument('--process', default='e+ e- > a > d d~', type=str,
                        help='Process definition to consider.')
    parser.add_argument('--cores', default=1, type=int,
                        help='Number of FORM cores')
    parser.add_argument('--restrict_card',
        default=pjoin(os.environ['MG5DIR'],'models','sm','restrict_no_widths.dat'), 
                        type=str, help='Model restriction card to consider.')
    args = parser.parse_args()

    if args.cores > 1:
        FORM_processing_options['cores'] = args.cores

    import alpha_loop.interface as interface
    cli = interface.alphaLoopInterface()

    cli.do_import('model %s'%args.model)
    computed_model = model_reader.ModelReader(cli._curr_model)
    computed_model.set_parameters_and_couplings(args.restrict_card)        
    process_definition=cli.extract_process(args.process, proc_number=0)

    FORM_processing_options['generate_renormalisation_graphs'] = False
    FORM_processing_options['generate_integrated_UV_CTs'] = False
    #FORM_processing_options['UV_min_dod_to_subtract'] = -2

    nestedbubbles = FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 2, 3),
            ('p4', 3, 4), ('p5', 2, 4), ('p6', 4, 1), ('q2', 4, 5)], "NestedBubbles", ['q1'], computed_model,
            loop_momenta_names=('p1', 'p2', 'p5'),overall_factor='k2.k2').produce_output() 

    #nestedbubbles = FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 2), ('p3', 2, 3),
    #        ('p4', 3, 2), ('p5', 3, 1), ('q2', 3, 4)], "NestedBubbles", ['q1'], computed_model,
    #        loop_momenta_names=('p1', 'p2', 'p3'),overall_factor='1').produce_output() 

    exit("A")

    # result is -2 Zeta[3] 3 Pi/(16 Pi^2)^3 = -5.75396*10^-6
    mercedes = FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 3), ('p3', 3, 6),
            ('p4', 6, 5), ('p5', 5, 1), ('p6', 2, 4), ('p7', 3, 4), ('p8', 4, 5), ('q2', 6, 7)], "Mercedes", ['q1'], computed_model,
            loop_momenta_names=('p1', 'p2', 'p3')).produce_output()

    # result is -5 Zeta[5] 4 Pi/(16 Pi^2)^4 = -1.04773*10^-7
    FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 2), ('p2', 2, 7), ('p3', 7, 3), ('p4', 3, 6),
        ('p5', 6, 5), ('p6', 5, 1), ('p7', 2, 4), ('p8', 3, 4), ('p9', 4, 5), ('p10', 7, 4), ('q2', 6, 8)],
        "DoubleMercedes", ['q1'], computed_model, loop_momenta_names=('p1', 'p2', 'p3', 'p4')).produce_output()

    # result is -5 /2 Zeta[5] 4 Pi/(16 Pi^2)^4 = -5.23865e-08
    FORMSuperGraphList.from_squared_topology([('q1', 0, 1), ('p1', 1, 3), ('p2', 1, 4), ('p3',2, 3), ('p4', 2, 5),
        ('p5', 3, 6), ('p6',4, 7), ('p7', 4, 8), ('p8', 5, 7), ('p9', 5, 8), ('p10', 6, 7), ('p11',6, 8), ('q2', 2, 9)],
        "STF4L17", ['q1'], computed_model).produce_output()
