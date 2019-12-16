#generate Latex table template

import json
import os
from uncertainties import *
import numpy
import argparse
import romanclass as roman
from copy import copy


FISHNET_CITATION = r'Basso:2017jwq'
FRANCESCO_CITATION = r'Frellesvig:2019byn'
SECDEC_CITATION = r'Borowka:2017idc'
WEINZIERL_CITATION = r'BeckerMultiLoop2012'
LADDER_CITATION = r'Usyukina:1992jd'


file_path = os.path.dirname(os.path.realpath( __file__ ))
#file_path = '/Users/Dario/Desktop/PhD/Code/pynloop/'
pjoin = os.path.join

ps_pts_per_topology = [2,3,2,2,1] #number of ps points per topology
columns = []
columns += [{'header': r'Topology',
			'shared': True,
			'json_link': lambda sample: sample.get('graph') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'Kin.',
			'shared': True,
			'json_link': lambda sample: sample.get('ps_point') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{C}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_cuts') if sample is not None else '?',
			'align': 'r'}]
columns += [{'header': r'$\mathtt{N_{\text{E}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_unique_existing_E_surface') if sample is not None else '?',
			'align': 'r'}]
columns += [{'header': r'$\mathtt{N_{\text{S}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_sources') if sample is not None else '?',
			'align': 'r'}]
columns += [{'header': r'$\mathtt{L_{\text{max}}}$',
			'shared': True,
			'json_link': lambda sample: format_overlap(sample.get('overlap_multiplicity')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{p}} \ [10^9]}$',
			'shared': True,
			'json_link': lambda sample: format_n_points(sample.get('num_samples')) if sample is not None else '?',
			'align': 'r'}]
columns += [{'header': r'$\mathtt{\sfrac{t}{p} \ [\mu \text{s}]}$',
			'shared': True,
			'json_link': lambda sample: format_timing(sample.get('t_per_ps_point_in_s')) if sample is not None else '?',
			'align': 'r'}]
columns += [{'header': r'Phase',
			'shared': False,
			'json_link': lambda sample: [r'Re',r'Im'],
			'align': 'c'}]
columns += [{'header': r'Exp.',
			'shared': True,
			'json_link': lambda sample: format_exponent(sample.get('result'), sample.get('error'), sample.get('analytical_result')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'Reference',
			'shared': False,
			'json_link': lambda sample: format_reference(sample.get('result'), sample.get('error'), sample.get('analytical_result')) if sample is not None else ['?','?'],
			'align': 'l'}]
columns += [{'header': r'Numerical LTD',
			'shared': False,
			'json_link': lambda sample: format_ltd(sample.get('result'), sample.get('error'), sample.get('analytical_result')) if sample is not None else ['?','?'],
			'align': 'c'}]
columns += [{'header': r'$\mathtt{\Delta \ [\sigma]}$',
			'shared': False,
			'json_link': lambda sample: format_delta(sample.get('accuracy')) if sample is not None else ['?','?'],
			'align': 'r'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%]}$',
			'shared': False,
			'json_link': lambda sample: format_delta(sample.get('percentage')) if sample is not None else ['?','?'],
			'align': 'r'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%] |\cdot|}$',
			'shared': True,
			'json_link': lambda sample: format_delta(sample.get('abs_error')) if sample is not None else '?',
			'align': 'r'}]

NHEADER = len(columns)

def format_overlap(overlap):
	if overlap is None:
		return '?'
	if overlap == []:
		return ''
	overlap_str = str(overlap)
	#overlap_str = overlap_str.replace(' ','')
	#half = 40.
	#center_index = 0
	#distance = 1000
	#for i,s in enumerate(overlap_str[:-1]):
	#	if s == ']' and overlap_str[i+1]==',':
	#		if abs(half-(i+1)) < distance:
	#			distance = abs(half-(i+1))
	#			center_index = i+1
	#if center_index == 0 or len(overlap_str) < half: 
	#	overlap_str = r'\mathtt{'+overlap_str+'}'
	#else:
	#	overlap_str = r'&\mathtt{'+overlap_str[:center_index+1]+r'}\\&\mathtt{'+overlap_str[center_index+1:]+r'}'
	#overlap_str = overlap_str.replace('[[','(')
	#overlap_str = overlap_str.replace('[','(')
	#overlap_str = overlap_str.replace(']]',')')
	#overlap_str = overlap_str.replace(']',')')
	return r'$\mathtt{'+overlap_str+'}$'

def format_n_points(n_points):
	if n_points is None:
		return '?'
	n_points_B = float(n_points)/10**9 
	if n_points % 10**9 == 0:
		return '{:d}'.format(int(n_points_B))+r'\phantom{.0}'
	else:
		return '{:.1f}'.format(n_points_B)
	#n_points_B = float(n_points)/10**9
	#if n_points_B>=1.:
	#	return '{:d}'.format(int(n_points_B))
	#else:
	#	return '{:.1f}'.format(n_points_B)

def format_timing(time):
	if time is None:
		return '?'
	return '{:.0f}'.format(time*10**6)

def determine_max_exp(results,errors,references):
	if not isinstance(references[0], str) and not isinstance(references[1], str):
		exponents_res = [int('{:e}'.format(results[0]).split('e')[1]),int('{:e}'.format(results[1]).split('e')[1])]
		exponents_err = [int('{:e}'.format(errors[0]).split('e')[1]),int('{:e}'.format(errors[1]).split('e')[1])]
		if results[0] == 0.:
			assert(errors[0]==0.)
			max_exp_res = exponents_res[1]
			max_exp_err = exponents_err[1]
		elif results[1] == 0.:
			assert(errors[1]==0.)
			max_exp_res = exponents_res[0]
			max_exp_err = exponents_err[0]
		else:
			max_exp_res = max(exponents_res)
			max_exp_err = max(exponents_err)
		max_exp_ltd = max(max_exp_res,max_exp_err)
		exponents_ref = [int('{:e}'.format(references[0]).split('e')[1]),int('{:e}'.format(references[1]).split('e')[1])]
		if references[0] == 0.:
			max_exp_ref = exponents_ref[1]
		elif references[1] == 0.:
			max_exp_ref = exponents_ref[0]
		else:
			max_exp_ref = max(exponents_ref)
		max_exp = max(max_exp_ref,max_exp_ltd)
	else:
		exponents_res = [int('{:e}'.format(results[0]).split('e')[1]),int('{:e}'.format(results[1]).split('e')[1])]
		max_exp_res = max(exponents_res)
		exponents_err = [int('{:e}'.format(errors[0]).split('e')[1]),int('{:e}'.format(errors[1]).split('e')[1])]
		max_exp_err = max(exponents_err)
		max_exp = max([max_exp_res,max_exp_err])
	return max_exp

def format_exponent(results,errors,references):
	max_exp = determine_max_exp(results,errors,references)
	return '{:e}'.format(10**max_exp).split('e')[1]

def format_reference(results,errors,references):
	references_str = []
	if not isinstance(references[0], str) and not isinstance(references[1], str):
		max_exp = determine_max_exp(results,errors,references)
	for phase in [0,1]: 
		reference = references[phase]
		if isinstance(reference, str):
			references_str += [r'\multicolumn{1}{c}{'+reference+'}']
			#references_str += ['~'+reference]
		else:
			if reference == 0.:
				references_str += [r'\phantom{+}\texttt{0}\phantom{.00000}']
			else:
				if isinstance(reference,UFloat):
					reference_str = r'{: .1u}'.format(reference/10**max_exp)
					if len(reference_str) < 17:
						diff = int((17 - len(reference_str))/2)
						reference_str = reference_str.replace('+/-',r'\phantom{'+'0'*diff+'}'+'+/-')
					reference_str = reference_str.replace('+/-','~+/-~')
					reference_str += r'\phantom{'+'0'*diff+'}'
				else:
					reference_str = r'{: .5f}'.format(reference/10**max_exp) #+ 'e'+'{:e}'.format(10**max_exp).split('e')[1]
					#reference_str = reference_str.replace(r'e',r'$\mathtt{\cdot 10^{')
					#reference_str += '}}$'
				references_str += [r'\texttt{'+reference_str+r'}']
	return references_str

def format_ltd(results,errors,references):
	results_str = []
	max_exp = determine_max_exp(results,errors,references)
	digits = 5
	for phase in [0,1]:
		if results[phase] == 0:
			assert(errors[phase]==0)
			results_str += ['']
		else:
			res = ufloat(results[phase], errors[phase])
			if res is None:
				return '?'
			#res_str = '{:ue}'.format(res)
			res_str = ('{:.'+str(digits)+'f}').format(res/10**max_exp) #+ 'e'+'{:e}'.format(10**max_exp).split('e')[1]
			
			#res_str = ('{:.1u}').format(res/10**max_exp) #+ 'e'+'{:e}'.format(10**max_exp).split('e')[1]
			#len_res_str = len(res_str)
			#if len_res_str < 2*digits + 9:
			#	assert(((2*digits + 9) - len_res_str) % 2 == 0)
			#	diff = ((2*digits + 9) - len_res_str)/2
			#	res_str = res_str.replace('+/-',r'\phantom{'+'0'*diff+'}')
			#if res_str[0] == r'(':
			#	if res_str[1] != r'-':
			#		res_str = r'(~' + res_str[1:]
			if not res_str[0] == r'-':
				res_str = r'~'+ res_str
			counter = 0
			dot = False
			for s in res_str:
				if s == '.':
					dot = True
				elif s == '+' and dot:
					break
				elif dot:
					counter += 1
			if '.' not in res_str:
				counter -= 1
			res_str = res_str.replace(r'+/-',r'~'*(digits+1-counter)+r'+/-~')
			#if ')e' in res_str:
			#	res_str = res_str.replace(r')e',r'~'*(digits-counter)+r')'+r'$\mathtt{\cdot 10^{')
			#	res_str += '}}$'
			#else:
			#	res_str += r'~'*(digits-counter)+r')'+r'$\mathtt{\phantom{{}\cdot{}10^{-10}}}$'
			results_str += [r'\texttt{'+res_str+r'}']
	return results_str

def format_delta(deltas):
	deltas_str = []
	phases = [0,1]
	if not isinstance(deltas,list):
		deltas = [deltas]
		phases = [0]
	for phase in phases:
		delta = deltas[phase]
		if isinstance(delta, str):
			deltas_str += [delta]
		else:
			if delta < 5e-4:
				deltas_str += ['{:.0e}'.format(delta)]
			elif delta >= 10. and delta < 100.:
				deltas_str += ['{:.2f}'.format(delta)]
			elif delta >= 100.:
				deltas_str += [r'{\color{red}'+'{:.0e}'.format(delta)+r'}']
			else:
				deltas_str += [r'{:.3f}'.format(delta)]
	if len(phases)==1:
		return deltas_str[0]
	return deltas_str

def set_header(columns):
	header_string = r'\hline' + '\n'
	for i,column in enumerate(columns):
		if i != len(columns)-1:
			header_string += column['header'] + r' & '
		else:
			header_string += column['header'] + r' \\' + '\n'
	return header_string

def set_topology(n_ps_points,data=None,add_real_ref=False):
	for i in range(n_ps_points):
		if i == 0:
			if add_real_ref:
				topology_string = r'\hline' + '\n'
				topology_string += r'\multirow{%i}{*}{'%(2*n_ps_points+1) + r'%' + '\n'
				topology_string += r'\begin{tabular}{@{}c@{}}'+ '\n'
				topology_string += r'\input{diagrams/'+data['diagram']+r'.tex} \\' + '\n'
				topology_name = str(columns[0]['json_link'](data)).replace('_',r'\_')
				topology_string += topology_name
				topology_string += r'\end{tabular}}' + '\n'
				if data['diagram_multi'] == 1:
					topology_string = topology_string.replace(r'}{*}',r'}{*}[-2pt]')
					topology_string = topology_string.replace(r'\end{tabular}}'+'\n',r'\end{tabular}}'+'\n'+r'\rule[-5pt]{0pt}{20pt}')
			else:
				if data['diagram'] is not None:
					assert(data['diagram_multi'] is not None)
					topology_string = r'\hline' + '\n'
					topology_string += r'\multirow{%i}{*}{'%(2*n_ps_points*data['diagram_multi']) + r'%' +'\n'
					topology_string += r'\begin{tabular}{@{}c@{}}'+ '\n'
					topology_string += r'\input{diagrams/'+data['diagram']+r'.tex} \\' + '\n'
					topology_name = str(columns[0]['json_link'](data)).replace('_',r'\_')
					topology_string += topology_name
					topology_string += r'\end{tabular}}'+'\n'
					#topology_string += r'\multirow{%i}{*}{'%(2*n_ps_points) +topology_name+'}' +'\n'
					if data['diagram_multi'] == 1:
						topology_string = topology_string.replace(r'}{*}',r'}{*}[7pt]')
						topology_string = topology_string.replace(r'\end{tabular}}'+'\n',r'\end{tabular}}'+'\n'+r'\rule[-5pt]{0pt}{27pt}')
				else:
					topology_string = r'\cline{2-%i}'%NHEADER +'\n'
		else:
			topology_string += r'\cline{2-%i}'%NHEADER +'\n'
		for column in columns[1:]:
			if add_real_ref:
				if column['shared']:
					column_str = str(column['json_link'](data))
					topology_string += r'& \multirow{3}{*}{'+column_str+'}'
				else:
					if column['header'] == 'Reference':
						column_str = str(column['json_link'](data)[0])
						topology_string +=  r'& ' + column_str
					else:
						column_str = str(column['json_link'](data)[0])
						topology_string +=  r'& \multirow{2}{*}{' + column_str	+'}'
			else:
				if column['shared']:
					column_str = str(column['json_link'](data))
					topology_string += r'& \multirow{2}{*}{'+column_str+'}'
				else:
					column_str = str(column['json_link'](data)[0])
					topology_string +=  r'& ' + column_str
			if column['header'] == 'Reference':
				if not add_real_ref:
					if data.get('citation') is not None:
						topology_string += r' \multirow{2}{*}{'+data.get('citation')+'}'
				else:
					if data.get('citation') is not None:
						topology_string += r' '+data.get('citation')
		topology_string += r'\\'+'\n'
		if add_real_ref:
			for column in columns[1:]:
				if column['header'] == 'Reference':
					column_str = format_reference(data['result'],data['error'],data['additional_reference'])[0]
					topology_string +=  r'& ' + column_str	+ ' ' + data.get('additional_citation')
				else:
					topology_string += r'& '
			topology_string += r'\\'+'\n'
		if data['diagram_multi'] == 1 and not add_real_ref:
			topology_string += r'\rule[-5pt]{0pt}{27pt}'
		if add_real_ref:
			topology_string += r'\rule[-5pt]{0pt}{20pt}'
		for column in columns[1:]:
			if column['shared']:
				topology_string += r'& '
			else:
				column_str = str(column['json_link'](data)[1])
				topology_string +=  r'& ' + column_str	
		topology_string += r'\\'+'\n'
	return topology_string

def set_table(columns,ps_pts_per_topology,data=None,exclude=None,nr=None):
	table_string = r'\begin{table}[tbp]'+'\n'
	table_string += r'\centering'+'\n'
	table_string += r'\resizebox{\columnwidth}{!}{%'+'\n'
	alignment = ''.join([column['align'] for column in columns])
	table_string += r'\texttt{%'+'\n'
	table_string += r'\begin{tabular}{@{}'+alignment+r'@{}}'+'\n'
	table_string += set_header(columns)
	if data is not None:
		for n_ps_points,sample in zip(ps_pts_per_topology,data):
			if exclude is None or sample['topology'] not in exclude:
				#print(topo_data['topology'])
				if sample.get('additional_reference') is not None:
					sample_string = set_topology(n_ps_points,data=sample,add_real_ref=True)
				else:
					sample_string = set_topology(n_ps_points,data=sample)
				if sample['topology'] in ['2L_4P_Ladder_PS3_massive','2L_4P_Ladder_PS1_massive','2L_4P_Ladder_PS2_massive_up']:
					sample_string = sample_string.replace(r'0}} \multi',r'}} \multi')
				table_string += sample_string
	else:
		for n_ps_points in ps_pts_per_topology:
			table_string += set_topology(n_ps_points,data=None)
	table_string += r'\end{tabular}'+'\n'
	table_string += r'}'+'\n'
	table_string += r'}'+'\n'
	if nr is None:
		table_string += r'\caption{\label{tab:i} Results}'+'\n'
	else:
		table_string += r'\caption{\label{tab:%i}}'%nr+'\n'
	table_string += r'\end{table}'+'\n'
	return table_string

def get_history(history_path):
    historical_data = []

    try:
        with open(history_path, 'r') as f:
            historical_data = json.load(f)
    except:
        pass

    return historical_data

def update_and_extract_data(data,ordered_topology_names):
	sorted_data = []
	found = False
	for topology_name in ordered_topology_names:
		for sample in data:
			if sample['topology'] == topology_name:
				#sample['ps_point'] = name_map['ps_point']
				#sample['graph'] = name_map['graph']
				sorted_data += [sample]
				data.remove(sample)
				found = True
		if not found:
			print(topology_name)
		found = False
	return sorted_data, data


def sort_and_split_data(data):
	current_data_length = len(data)

	topologies_in_order = [	'Box_4E_no_z_component_sources',
							'Pentagon_1s',
							'Pentagon_10E_1s',
							'Pentagon_2s',
							'Pentagon_3s',
							'Pentagon_6E_4s',
							'Pentagon_8E_5s',
							'Hexagon_1s',
							'Hexagon_6E_2s',
							'Hexagon_2s',
							'Hexagon_3s',
							'Hexagon_6E_4s',
							'Hexagon_9E_4s',
							'Hexagon_10E_4s',
							'Hexagon_10E_7s',
							'Hexagon_4s',
							'Hexagon_10E_5s']
	
	explore_1loop_3B, data = update_and_extract_data(data,topologies_in_order)
	assert(current_data_length-len(topologies_in_order)==len(data))
	current_data_length = len(data)

	topologies_in_order = [	'1L_4P_PS3',
							'1L_4P_PS1',
							'1L_4P_PS2',
							'1L_4P_PS3_massive',
							'1L_4P_PS1_massive',
							'1L_4P_PS2_massive',
							'1L_4P_PS2_massive_down',
							'1L_4P_PS2_massive_up',
							'1L_5P_PS3',
							'1L_5P_PS1',
							'1L_5P_PS2',
							'1L_5P_PS3_massive',
							'1L_5P_PS1_massive',
							'1L_5P_PS2_massive',
							'1L_5P_PS2_massive_down',
							'1L_5P_PS2_massive_up',
							'1L_6P_PS3',
							'1L_6P_PS1',
							'1L_6P_PS2',
							'1L_6P_PS3_massive',
							'1L_6P_PS1_massive',
							'1L_6P_PS2_massive',
							'1L_6P_PS2_massive_up',
							'1L_6P_PS2_massive_down',
							'1L_8P_PS3',
							'1L_8P_PS1',
							'1L_8P_PS2',
							'1L_8P_PS3_massive',
							'1L_8P_PS1_massive',
							'1L_8P_PS2_massive',
							'1L_8P_PS2_massive_rescaled_down',
							'1L_8P_PS2_massive_rescaled_up']

	systematic_1loop_3B, data = update_and_extract_data(data,topologies_in_order)
	assert(current_data_length-len(topologies_in_order)==len(data))
	current_data_length = len(data)

	systematic_1loop_3B_1 = systematic_1loop_3B[:16]
	systematic_1loop_3B_2 = systematic_1loop_3B[16:]

	topologies_in_order = [	'TM1_top',
							'TM1_bot',
							'T3_DoubleBox_Weinzierl',
							'T2_6P_2L_Weinzierl_A',
							'T2_6P_2L_Weinzierl_B',
							'T2_6P_2L_Weinzierl_C',
							'T2_6P_2L_Weinzierl_D',
							'T2_6P_2L_Weinzierl_E',
							'T2_6P_2L_Weinzierl_F',
							'T4_TripleBox_Weinzierl',
							'FISHNET_2x2',
							'T4_Quadruple_Box_Weinzierl',
							'FISHNET_1x5',
							'FISHNET_2x3',
							'FISHNET_1x6']

	explore_HigherLoop, data = update_and_extract_data(data,topologies_in_order)
	assert(current_data_length-len(topologies_in_order)==len(data))
	current_data_length = len(data)
	explore_HigherLoop_1 = explore_HigherLoop[:9]
	explore_HigherLoop_2 = explore_HigherLoop[9:]

	topologies_in_order = [	'2L_4P_Ladder_PS3',
							'2L_4P_Ladder_PS1',
							'2L_4P_Ladder_PS2',
							'2L_4P_Ladder_PS3_massive',
							'2L_4P_Ladder_PS1_massive',
							'2L_4P_Ladder_PS2_massive_down',
							'2L_4P_Ladder_PS2_massive_up',
							#'2L_4P_Ladder_PS2_massive',
							'2L_5P_Planar_PS3',
							'2L_5P_Planar_PS1',
							'2L_5P_Planar_PS2',
							'2L_5P_Planar_PS3_massive',
							'2L_5P_Planar_PS1_massive',
							'2L_5P_Planar_PS2_massive_down',
							'2L_5P_Planar_PS2_massive_up',
							#'2L_5P_Planar_PS2_massive',
							'2L_6P_A_PS3',
							'2L_6P_A_PS1',
							#'2L_6P_A_PS2',
							'2L_6P_A_PS3_massive',
							'2L_6P_A_PS1_massive',
							'2L_6P_A_PS2_massive_down',
							#'2L_6P_A_PS2_massive',
							'2L_6P_B_PS3',
							'2L_6P_B_PS1',
							'2L_6P_B_PS2',
							'2L_6P_B_PS3_massive',
							'2L_6P_B_PS1_massive',
							'2L_6P_B_PS2_massive_up',
							'2L_6P_B_PS2_massive_down',
							#'2L_6P_B_PS2_massive',
							'2L_6P_C_PS3',
							'2L_6P_C_PS1',
							'2L_6P_C_PS2',
							'2L_6P_C_PS3_massive',
							'2L_6P_C_PS1_massive',
							'2L_6P_C_PS2_massive_up',
							'2L_6P_C_PS2_massive_down',
							#'2L_6P_C_PS2_massive',
							'2L_6P_D_PS3',
							'2L_6P_D_PS1',
							'2L_6P_D_PS2',
							'2L_6P_D_PS3_massive',
							'2L_6P_D_PS1_massive',
							'2L_6P_D_PS2_massive_up',
							'2L_6P_D_PS2_massive_down',
							#'2L_6P_D_PS2_massive',
							'2L_6P_E_PS3',
							'2L_6P_E_PS1',
							'2L_6P_E_PS2',
							'2L_6P_E_PS3_massive',
							'2L_6P_E_PS1_massive',
							'2L_6P_E_PS2_massive_up',
							'2L_6P_E_PS2_massive_down',
							#'2L_6P_E_PS2_massive',
							'2L_6P_F_PS3',
							'2L_6P_F_PS1',
							'2L_6P_F_PS2',
							'2L_6P_F_PS3_massive',
							'2L_6P_F_PS2_massive_up',
							'2L_6P_F_PS2_massive_down',
							#'2L_6P_F_PS1_massive',
							#'2L_6P_F_PS2_massive',
							'2L_8P_PS3',
							#'2L_8P_PS1',
							#'2L_8P_PS2',
							'3L_4P_Ladder_PS3',
							'3L_4P_Ladder_PS1',
							'3L_4P_Ladder_PS2',
							'3L_4P_Ladder_PS3_massive',
							'3L_4P_Ladder_PS1_massive',
							'3L_4P_Ladder_PS2_massive_down',
							'3L_4P_Ladder_PS2_massive_up',
							#'3L_4P_Ladder_PS2_massive',
							'3L_5P_Planar_PS3',
							'3L_5P_Planar_PS1',
							'3L_5P_Planar_PS2',
							'3L_5P_Planar_PS3_massive',
							'3L_5P_Planar_PS1_massive',
							'3L_5P_Planar_PS2_massive_down',
							'3L_5P_Planar_PS2_massive_up',
							#'3L_5P_Planar_PS2_massive',
							'FISHNET_2x2_PS3',
							'FISHNET_2x2_PS3_massive',
							#'FISHNET_2x2_PS1_massive',
							#'FISHNET_2x2_PS2_massive',
							'4L_4P_Ladder_PS3',
							'4L_4P_Ladder_PS1',
							'4L_4P_Ladder_PS2',
							'4L_4P_Ladder_PS3_massive',
							'4L_4P_Ladder_PS1_massive',
							#'4L_4P_Ladder_PS2_massive',
							]

	systematic_HigherLoop, data = update_and_extract_data(data,topologies_in_order)

	systematic_HigherLoop_1 = systematic_HigherLoop[:14]
	systematic_HigherLoop_2 = systematic_HigherLoop[14:54]
	systematic_HigherLoop_3 = systematic_HigherLoop[54:]
	
	for sample in data:
		print(sample['topology'])
	assert(data==[])

	return [explore_1loop_3B, explore_HigherLoop_1,explore_HigherLoop_2, systematic_1loop_3B_1, systematic_1loop_3B_2, systematic_HigherLoop_1,systematic_HigherLoop_2,systematic_HigherLoop_3]

def get_all_data(paths):
	all_sample_data = []
	for path in all_paths:
		historical_data = get_history(path)
		for t, v in historical_data.items():
			samples = v['samples']
			for sample in samples:
				if sample.get('n_cuts') is not None:
					all_sample_data += [sample]
	return all_sample_data

def complete_data(data):
	for sample in data:
		assert(sample.get('maximal_overlap') is not None)
		if sample.get('overlap_multiplicity') is None:
			sample['overlap_multiplicity'] = [len(overlap) for overlap in sample['maximal_overlap']]
		if isinstance(sample['analytical_result'][0], UFloat):
			analytic_res_re = sample['analytical_result'][0].nominal_value
		else:
			analytic_res_re = sample['analytical_result'][0]
		if isinstance(sample['analytical_result'][1], UFloat):
			analytic_res_im = sample['analytical_result'][1].nominal_value
		else:
			analytic_res_im = sample['analytical_result'][1]
		analytic_res = [analytic_res_re,analytic_res_im]
		if ([analytic_res[0]-sample['result'][0],analytic_res[1]-sample['result'][1]] == [0.,0.]
			or isinstance(sample['analytical_result'][0], UFloat) or isinstance(sample['analytical_result'][1], UFloat)):
			if [analytic_res[0]-sample['result'][0],analytic_res[1]-sample['result'][1]] == [0.,0.]:
				sample['analytical_result'] = ['n/a','n/a']
			sample['accuracy'] = [' ',' ']
			sample['percentage'] = [100.0*(abs(sample['error'][i_phase]) / abs(sample['result'][i_phase])) for i_phase in [0,1]]
			sample['abs_error'] = 100.* (numpy.sqrt(
					(sample['error'][0])**2 + (sample['error'][1])**2
					)/numpy.sqrt(sample['result'][0]**2 + sample['result'][1]**2))
		else:
			sample['accuracy'] = [abs(analytic_res[i_phase] - sample['result'][i_phase]) / sample['error'][i_phase]
									if sample['error'][i_phase] != 0. else ' ' for i_phase in [0,1]]
			assert(sample['accuracy'] != [0.,0.])
			sample['percentage'] = [100.0*(abs(analytic_res[i_phase]-sample['result'][i_phase]) / abs(analytic_res[i_phase]))
					if abs(analytic_res[i_phase]) != 0. else ' ' for i_phase in [0,1]]
			assert(sample['percentage'] != [0.,0.])
			assert(numpy.sqrt(analytic_res[0]**2 + analytic_res[1]**2) != 0.)
			sample['abs_error'] = 100.* (numpy.sqrt(
					(analytic_res[0]-sample['result'][0])**2 + (analytic_res[1]-sample['result'][1])**2
					)/numpy.sqrt(analytic_res[0]**2 + analytic_res[1]**2))
			assert(sample['abs_error'] != 0.)

	for sample in data:
		topology_name = sample['topology']
		if 'Box_4E' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = 'Box4E'
		elif 'Pentagon' in topology_name:
			if 'Pentagon_1s' in topology_name:
				sample['ps_point'] = str(roman.Roman(1))
			elif '10E_1s' in topology_name:
				sample['ps_point'] = str(roman.Roman(2))
			elif '2s' in topology_name:
				sample['ps_point'] = str(roman.Roman(3))
			elif '3s' in topology_name:
				sample['ps_point'] = str(roman.Roman(4))
			elif '4s' in topology_name:
				sample['ps_point'] = str(roman.Roman(5))
			elif '5s' in topology_name:
				sample['ps_point'] = str(roman.Roman(6))
			sample['graph'] = '1L5P'
		elif 'Hexagon' in topology_name:
			if '1s' in topology_name:
				sample['ps_point'] = str(roman.Roman(1))
			elif '6E_2s' in topology_name:
				sample['ps_point'] = str(roman.Roman(2))
			elif '2s' in topology_name:
				sample['ps_point'] = str(roman.Roman(3))
			elif '3s' in topology_name:
				sample['ps_point'] = str(roman.Roman(4))
			elif '6E_4s' in topology_name:
				sample['ps_point'] = str(roman.Roman(5))
			elif '9E_4s' in topology_name:
				sample['ps_point'] = str(roman.Roman(6))
			elif '10E_4s' in topology_name:
				sample['ps_point'] = str(roman.Roman(7))
			elif '10E_7s' in topology_name:
				sample['ps_point'] = str(roman.Roman(8))
			elif '4s' in topology_name:
				sample['ps_point'] = str(roman.Roman(9))
			elif '5s' in topology_name:
				sample['ps_point'] = str(roman.Roman(10))
			sample['graph'] = '1L6P'
		elif '1L_4P' in topology_name:
			sample['graph'] = '1L4P'
		elif '1L_5P' in topology_name:
			sample['graph'] = '1L5P'
		elif '1L_6P' in topology_name:
			sample['graph'] = '1L6P'
		elif '1L_8P' in topology_name:
			sample['graph'] = '1L8P'
		elif 'DoubleBox' in topology_name or '2L_4P_Ladder' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L4P.b'
			if not 'massive' in topology_name:
				sample['citation'] = r'\cite{'+LADDER_CITATION+'}'
			else:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
		elif '6P_2L_Weinzierl_A' in topology_name or '2L_6P_A' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.a'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif '6P_2L_Weinzierl_B' in topology_name or '2L_6P_B' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.b'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif '6P_2L_Weinzierl_C' in topology_name or '2L_6P_C' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.c'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif '6P_2L_Weinzierl_D' in topology_name or '2L_6P_D' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.d'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif '6P_2L_Weinzierl_E' in topology_name or '2L_6P_E' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.e'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif '6P_2L_Weinzierl_F' in topology_name or '2L_6P_F' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '2L6P.f'
			if 'Weinzierl' in topology_name:
				sample['citation'] = r'\cite{'+SECDEC_CITATION+'}'
				sample['additional_citation'] = r'\cite{'+WEINZIERL_CITATION+'}'
		elif 'TripleBox' in topology_name or '3L_4P' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '3L4P'
			if not 'massive' in topology_name:
				sample['citation'] = r'\cite{'+LADDER_CITATION+'}'
		elif 'Quadruple_Box' in topology_name or '4L_4P_Ladder' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '4L4P.b'
			if not 'massive' in topology_name:
				sample['citation'] = r'\cite{'+LADDER_CITATION+'}'
		elif 'TM1' in topology_name:
			if 'TM1_top' in topology_name:
				sample['ps_point'] = str(roman.Roman(1))
			elif 'TM1_bot' in topology_name:
				sample['ps_point'] = str(roman.Roman(2))
			sample['graph'] = '2L4P.a'
			sample['citation'] = r'\cite{'+FRANCESCO_CITATION+'}'
		elif '2L_5P' in topology_name:
			sample['graph'] = '2L5P'
		elif '2L_8P' in topology_name:
			sample['graph'] = '2L8P'
		elif '3L_5P' in topology_name:
			sample['graph'] = '3L5P'
		elif 'FISHNET_2x2' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '4L4P.a'
			if not 'massive' in topology_name:
				sample['citation'] = r'\cite{'+FISHNET_CITATION+'}'
		elif 'FISHNET_1x5' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '5L4P'
			sample['citation'] = r'\cite{'+FISHNET_CITATION+'}'
		elif 'FISHNET_1x6' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '6L4P.b'
			sample['citation'] = r'\cite{'+FISHNET_CITATION+'}'
		elif 'FISHNET_2x3' in topology_name:
			sample['ps_point'] = str(roman.Roman(1))
			sample['graph'] = '6L4P.a'
			sample['citation'] = r'\cite{'+FISHNET_CITATION+'}'
		else:
			print('topology not found!')
			assert(1==2)
		if 'PS3' in topology_name:
			sample['ps_point'] = 'K1'
		elif 'PS1' in topology_name:
			sample['ps_point'] = 'K2'
		elif 'PS2' in topology_name:
			sample['ps_point'] = 'K3'
		if 'massive' in topology_name:
			sample['ps_point'] += r'$^*$'

	for sample in data:
		topology_name = sample['topology']
		if ('Box_4E' in topology_name or '1L_4P_PS3' in topology_name) and not 'massive' in topology_name:
			sample['diagram'] = 'Box'
			if 'Box_4E' in topology_name:
				sample['diagram_multi'] = 1
			elif '1L_4P_PS3' in topology_name:
				sample['diagram_multi'] = 6
		elif ('Pentagon_1s' in topology_name or '1L_5P_PS3' in topology_name) and not 'massive' in topology_name:
			sample['diagram'] = 'Pentagon'
			if 'Pentagon_1s' in topology_name:
				sample['diagram_multi'] = 6
			elif '1L_5P_PS3' in topology_name:
				sample['diagram_multi'] = 6
		elif ('Hexagon_1s' in topology_name or '1L_6P_PS3' in topology_name) and not 'massive' in topology_name:
			sample['diagram'] = 'Hexagon'
			if 'Hexagon_1s' in topology_name:
				sample['diagram_multi'] = 10
			elif '1L_6P_PS3' in topology_name:
				sample['diagram_multi'] = 6
		elif '1L_8P_PS3' in topology_name and not 'massive' in topology_name:
				sample['diagram'] = 'Octagon'
				sample['diagram_multi'] = 6
		elif 'TM1_top' in topology_name:
			sample['diagram'] = '2L_4P_AB'
			sample['diagram_multi'] = 2
		elif ('T3_DoubleBox_Weinzierl' in topology_name or '2L_4P_Ladder_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_4P'
				if 'T3_DoubleBox_Weinzierl' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_4P_Ladder_PS3' in topology_name:
					sample['diagram_multi'] = 6
		elif ('T2_6P_2L_Weinzierl_A' in topology_name or '2L_6P_A_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_A'
				if 'T2_6P_2L_Weinzierl_A' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_A_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T2_6P_2L_Weinzierl_B' in topology_name or '2L_6P_B_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_B'
				if 'T2_6P_2L_Weinzierl_B' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_B_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T2_6P_2L_Weinzierl_C' in topology_name or '2L_6P_C_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_C'
				if 'T2_6P_2L_Weinzierl_C' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_C_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T2_6P_2L_Weinzierl_D' in topology_name or '2L_6P_D_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_D'
				if 'T2_6P_2L_Weinzierl_D' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_D_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T2_6P_2L_Weinzierl_E' in topology_name or '2L_6P_E_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_E'
				if 'T2_6P_2L_Weinzierl_E' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_E_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T2_6P_2L_Weinzierl_F' in topology_name or '2L_6P_F_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '2L_6P_F'
				if 'T2_6P_2L_Weinzierl_F' in topology_name:
					sample['diagram_multi'] = 1
				elif '2L_6P_F_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif ('T4_TripleBox_Weinzierl' in topology_name or '3L_4P_Ladder_PS3' in topology_name) and not 'massive' in topology_name:
			sample['diagram'] = '3L_4P'
			if 'T4_TripleBox_Weinzierl' in topology_name:
				sample['diagram_multi'] = 1
			elif '3L_4P_Ladder_PS3' in topology_name:
				sample['diagram_multi'] = 6
		elif ('FISHNET_2x2' in topology_name or 'FISHNET_2x2_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '4L_4P_A'
				if 'FISHNET_2x2_PS3' in topology_name:
					sample['diagram_multi'] = 2
				else:
					sample['diagram_multi'] = 1
		elif ('T4_Quadruple_Box_Weinzierl' in topology_name or '4L_4P_Ladder_PS3' in topology_name) and not 'massive' in topology_name:
				sample['diagram'] = '4L_4P_B'
				if 'T4_Quadruple_Box_Weinzierl' in topology_name:
					 sample['diagram_multi'] = 1
				elif '4L_4P_Ladder_PS3' in topology_name:
					sample['diagram_multi'] = 2
		elif 'FISHNET_1x5' in topology_name:
			sample['diagram'] = '5L_4P'
			sample['diagram_multi'] = 1
		elif 'FISHNET_1x6' in topology_name:
			sample['diagram'] = '6L_4P_B'
			sample['diagram_multi'] = 1
		elif 'FISHNET_2x3' in topology_name:
			sample['diagram'] = '6L_4P_A'
			sample['diagram_multi'] = 1
		elif '2L_5P_Planar_PS3' in topology_name and not 'massive' in topology_name:
				sample['diagram'] = '2L_5P'
				sample['diagram_multi'] = 6
		elif '2L_8P_PS3' in topology_name:
			sample['diagram'] = '2L_8P'
			sample['diagram_multi'] = 1
		elif '3L_5P_Planar_PS3' in topology_name and not 'massive' in topology_name:
				sample['diagram'] = '3L_5P'
				sample['diagram_multi'] = 6
		else:
			sample['diagram'] = None
			sample['diagram_multi'] = None

	return data

if __name__ == "__main__":

	all_paths = []
	all_paths += [pjoin(file_path,"deformation_paper_results/explore_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS1PS2_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS3_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/explore_HigherLoop.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS1PS2_HigherLoop.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS3_HigherLoop.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/missing_loops.json")]

	all_sample_data = get_all_data(all_paths)

	manual_FISHNET_1x5_sample = {'num_samples': int(1.8e9), 
								'topology': 'FISHNET_1x5',
								'result': [0., 3.289e-16],
								'error': [0., 1.96409e-18],
								'analytical_result': [0., 3.31696678e-16],
								'n_unique_existing_E_surface': 0,
								'n_sources': 0,
								'maximal_overlap': [],
								'max_radius': 0.,
								'min_radius': 0.,
								'n_cuts': 780,
								't_per_ps_point_in_s': 127.542888*10**(-3)/500}

	manual_FISHNET_1x6_sample = {'num_samples': int(1e9), 
								'topology': 'FISHNET_1x6',
								'result': [1.09968e-18,0.],
								'error': [4.17287e-19,0.],
								'analytical_result': [9.060040029310959e-19,0.],
								'n_unique_existing_E_surface': 0,
								'n_sources': 0,
								'maximal_overlap': [],
								'max_radius': 0.,
								'min_radius': 0.,
								'n_cuts': 2911,
								't_per_ps_point_in_s': 600.14064*10**(-3)/500}

	manual_FISHNET_2x3_sample = {'num_samples': int(14.5e9), 
								'topology': 'FISHNET_2x3',
								'result': [8.36493e-19,0.],
								'error': [2.16724e-21,0.],
								'analytical_result': [8.4044862640909e-19,0.],
								'n_unique_existing_E_surface': 0,
								'n_sources': 0,
								'maximal_overlap': [],
								'max_radius': 0.,
								'min_radius': 0.,
								'n_cuts': 2415,
								't_per_ps_point_in_s': 598.039545*10**(-3)/500}

	manual_FISHNET_2x2_sample = {'num_samples': int(667000000), 
								'topology': 'FISHNET_2x2',
								'result': [8.38828e-05,-1.0278e-07],
								'error': [7.77242e-07,7.75435e-07],
								'analytical_result': [8.416099347763927e-5,0],
								'n_unique_existing_E_surface': 44,
								'n_sources': 280,
								'maximal_overlap': [],
								'overlap_multiplicity': [44],
								'max_radius': 0.31933513395572766,
								'min_radius': 0.00018338481862451912,
								'n_cuts': 192,
								't_per_ps_point_in_s': 0.}

	placeholder_FISHNET_2x2_PS3_sample = {'num_samples': int(0), 
								'topology': 'FISHNET_2x2_PS3',
								'result': [99.,99.],
								'error': [99.,99.],
								'analytical_result': [99.,99.],
								'n_unique_existing_E_surface': 0,
								'n_sources': 0,
								'maximal_overlap': [],
								'max_radius': 0.,
								'min_radius': 0.,
								'n_cuts': 0,
								't_per_ps_point_in_s': 0.}

	placeholder_4L_4P_Ladder_PS3_sample = {'num_samples': int(0), 
								'topology': '4L_4P_Ladder_PS3',
								'result': [99.,99.],
								'error': [99.,99.],
								'analytical_result': [99.,99.],
								'n_unique_existing_E_surface': 0,
								'n_sources': 0,
								'maximal_overlap': [],
								'max_radius': 0.,
								'min_radius': 0.,
								'n_cuts': 0,
								't_per_ps_point_in_s': 0.}

	all_sample_data += [	manual_FISHNET_1x5_sample,
							manual_FISHNET_1x6_sample,
							manual_FISHNET_2x3_sample,
							manual_FISHNET_2x2_sample,
							placeholder_FISHNET_2x2_PS3_sample,
							placeholder_4L_4P_Ladder_PS3_sample]

	
	for sample in all_sample_data:
		if sample['topology'] == 'T2_6P_2L_Weinzierl_A':
			sample['analytical_result'] = [ufloat(-86.0768134710165628, 0.0858487113390352118),
												0.]
			sample['additional_reference'] = [ufloat( -8.66, 0.08)*10,0.]
		elif sample['topology'] == 'T2_6P_2L_Weinzierl_B':
			sample['analytical_result'] = [ufloat(-118.862749862757883,0.0528974822691667579),
												0.]
			sample['additional_reference'] = [ufloat(-1.17,0.02)*10**2,0.]
		elif sample['topology'] == 'T2_6P_2L_Weinzierl_C':
			sample['analytical_result'] = [ufloat(-76.0653436805269237,0.0565739838497340792),
											0.]
			sample['additional_reference'] = [ufloat(-7.75, 0.13)*10,0.]
		elif sample['topology'] == 'T2_6P_2L_Weinzierl_D':
			sample['analytical_result'] = [ufloat(-18.3271249575150019, 0.0171281153680031097),
											0.]
			sample['additional_reference'] = [ufloat(-1.91, 0.02)*10,0.]
		elif sample['topology'] == 'T2_6P_2L_Weinzierl_E':
			sample['analytical_result'] = [ufloat(-45.9727585576039985, 0.0443921298216177022),
											0.]
			sample['additional_reference'] = [ufloat(-4.64, 0.08)*10,0.]
		elif sample['topology'] == 'T2_6P_2L_Weinzierl_F':
			sample['analytical_result'] = [ufloat(-102.713355708622871, 0.0288577968836020429),
											0.]
			sample['additional_reference'] = [ufloat(-1.03, 0.02)*10**2,0.]
		elif sample['topology'] == '2L_4P_Ladder_PS3_massive':
			sample['analytical_result'] = [ufloat(2.80171724275311214e-6, 7.79255366562339195e-9),
											ufloat(3.34481785517457328e-6, 7.72520711768175799e-9)]
		elif sample['topology'] == '2L_4P_Ladder_PS1_massive':
			sample['analytical_result'] = [ufloat(7.89148979543551809e-8, 7.08734798497054057e-9),
											ufloat(6.87085771146301916e-8, 7.07933376620586562e-9)]
		elif sample['topology'] == '2L_4P_Ladder_PS2_massive_up':
			sample['analytical_result'] = [ufloat(3.05451908291667266e-10, 1.15482681073316096e-11),
											ufloat(7.4245774920915505e-12, 1.15823697424355456e-11)]

	all_sample_data = complete_data(all_sample_data)

	for sample in all_sample_data:
		if sample['topology'] == '4L_4P_Ladder_PS3':
			if sample['result'] == [3.619624658469241e-12, -2.0527445803375738e-12]:
				all_sample_data.remove(sample)
			elif sample['result'] == [4.308707534909754e-12, -2.5066951276497768e-12]:
				all_sample_data.remove(sample)

			#print(sample)
	#stop
	#all_sample_data.remove({'revision': 'e58c4b7', 'diff': False, 'num_samples': 500000000, 'topology': '4L_4P_Ladder_PS3', 'result': [3.619624658469241e-12, -2.0527445803375738e-12], 'error': [9.950475725998295e-13, 9.897342748196721e-13], 'analytical_result': [4.03555168296597e-12, -1.0244887725916831e-12], 'n_unique_existing_E_surface': 28, 'n_sources': 270, 'maximal_overlap': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]], 'max_radius': 0.6782331861809063, 'min_radius': 0.10393659196076636, 'n_cuts': 209, 't_per_ps_point_in_s': 0.002752984665, 'overlap_multiplicity': [28], 'accuracy': [0.41799712491133195, 1.0389210860998392], 'percentage': [10.306571620736593, 100.36769901779185], 'abs_error': 26.64043635840804, 'ps_point': 'K1', 'graph': '4L4P.b', 'citation': 'Ladder', 'diagram': '4L_4P_B', 'diagram_multi': 2}
	#)
	#all_sample_data.remove({'revision': 'e58c4b7', 'diff': False, 'num_samples': 400000000, 'topology': '4L_4P_Ladder_PS3', 'result': [4.308707534909754e-12, -2.5066951276497768e-12], 'error': [1.0527262612060722e-12, 1.0918350285928435e-12], 'analytical_result': [4.03555168296597e-12, -1.0244887725916831e-12], 'n_unique_existing_E_surface': 28, 'n_sources': 270, 'maximal_overlap': [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]], 'max_radius': 0.6782331861809063, 'min_radius': 0.10393659196076636, 'n_cuts': 209, 't_per_ps_point_in_s': 0.003459343254, 'overlap_multiplicity': [28], 'accuracy': [0.25947472007665046, 1.357536913766506], 'percentage': [6.768736306779889, 144.67765725811782], 'abs_error': 36.198953754857264, 'ps_point': 'K1', 'graph': '4L4P.b', 'citation': 'Ladder', 'diagram': '4L_4P_B', 'diagram_multi': 2}
	#	)

	exclude_topologies = [	'2L_6P_A_PS1',			
							'2L_6P_A_PS1_massive',
							'2L_6P_B_PS1',		
							'2L_6P_B_PS2',			
							'2L_6P_B_PS1_massive',
							'2L_6P_C_PS1',			
							'2L_6P_C_PS2',			
							'2L_6P_C_PS1_massive',
							'2L_6P_D_PS1',			
							'2L_6P_D_PS2',			
							'2L_6P_D_PS1_massive',
							'2L_6P_E_PS1',			
							'2L_6P_E_PS2',			
							'2L_6P_E_PS1_massive',
							'2L_6P_F_PS1',			
							'2L_6P_F_PS2',
							'1L_4P_PS2_massive',
							'1L_5P_PS2_massive',
							'1L_6P_PS2_massive',
							'1L_8P_PS2_massive',
							'1L_4P_PS2_massive_down',
							'1L_5P_PS2_massive_down',
							'1L_6P_PS2_massive_down',
							'1L_8P_PS2_massive_rescaled_down',
							'2L_4P_Ladder_PS2_massive_down',
							'2L_5P_Planar_PS2_massive_down',
							'2L_6P_A_PS2_massive_down',
							'2L_6P_B_PS2_massive_down',
							'2L_6P_C_PS2_massive_down',
							'2L_6P_D_PS2_massive_down',
							'2L_6P_E_PS2_massive_down',
							'2L_6P_F_PS2_massive_down',
							'2L_6P_B_PS2_massive_up',
							'2L_6P_C_PS2_massive_up',
							'2L_6P_D_PS2_massive_up',
							'2L_6P_E_PS2_massive_up',
							'2L_6P_F_PS2_massive_up',
							'3L_4P_Ladder_PS2_massive_down',
							'3L_5P_Planar_PS2_massive_down',
							'4L_4P_Ladder_PS1',
							'4L_4P_Ladder_PS2',
							'4L_4P_Ladder_PS1_massive']
	
	all_topology_names = [sample['topology'] for sample in all_sample_data]
	for topology_name in exclude_topologies:
		assert(topology_name in all_topology_names)

	duplicates = [topology_name for i, topology_name in enumerate(all_topology_names) if topology_name in all_topology_names[:i]]
	#print(duplicates)
	assert(duplicates==[])

	data_length = len(all_sample_data)
	data_sets = sort_and_split_data(all_sample_data)
	assert(data_length-sum([len(data_set) for data_set in data_sets]) == 0)

	table_string = ''
	for table_nr,data_set in enumerate(data_sets):
		table_string += set_table(columns,[1 for i in range(len(data_set))],data=data_set,exclude=exclude_topologies,nr=table_nr)
		table_string += '\n'+'\n'

	table_captions = [r'\caption{\label{tab:explore_1loop} Explore 1loop}',
					r'\caption{\label{tab:explore_higherloopA} Explore higherloop A}',
					r'\caption{\label{tab:explore_higherloopB} Explore higherloop B}',
					r'\caption{\label{tab:Kseries_1loopA} Kseries 1loop A}',
					r'\caption{\label{tab:Kseries_1loopB} Kseries 1loop B}',
					r'\caption{\label{tab:Kseries_2loopA} Results Kseries 2loop A}',
					r'\caption{\label{tab:Kseries_2loopB} Results Kseries 2loop B}',
					r'\caption{\label{tab:Kseries_HigherLoop} Results Kseries higher loop}']

	for table_nr, table_caption in enumerate(table_captions):
		table_string = table_string.replace(r'\caption{\label{tab:%i}}'%table_nr,table_caption)

	table_string = table_string.replace(r'Reference',r'\multicolumn{1}{c}{Reference}')

	print(table_string)









