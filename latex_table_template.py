#generate Latex table template

import json
import os
from uncertainties import ufloat
import numpy
import argparse
import romanclass as roman

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
			'json_link': lambda sample: sample.get('PS_point') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{C}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_cuts') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{E}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_unique_existing_E_surface') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{S}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_sources') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{L_{\text{max}}}$',
			'shared': True,
			'json_link': lambda sample: format_overlap(sample.get('overlap_multiplicity')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{p}} \ [10^9]}$',
			'shared': True,
			'json_link': lambda sample: format_n_points(sample.get('num_samples')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\sfrac{t}{p} \ [\mu \text{s}]}$',
			'shared': True,
			'json_link': lambda sample: format_timing(sample.get('t_per_ps_point_in_s')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'Phase',
			'shared': False,
			'json_link': lambda sample: [r'$\Re$',r'$\Im$'],
			'align': 'l'}]
columns += [{'header': r'Exp.',
			'shared': True,
			'json_link': lambda sample: format_exponent(sample.get('result'), sample.get('error'), sample.get('analytical_result')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'Reference',
			'shared': False,
			'json_link': lambda sample: format_reference(sample.get('analytical_result')) if sample is not None else ['?','?'],
			'align': 'l'}]
columns += [{'header': r'Numerical LTD',
			'shared': False,
			'json_link': lambda sample: format_ltd(sample.get('result'), sample.get('error')) if sample is not None else ['?','?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\sigma]}$',
			'shared': False,
			'json_link': lambda sample: format_delta(sample.get('accuracy')) if sample is not None else ['?','?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%]}$',
			'shared': False,
			'json_link': lambda sample: format_delta(sample.get('percentage')) if sample is not None else ['?','?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%] |\cdot|}$',
			'shared': True,
			'json_link': lambda sample: format_delta(sample.get('abs_error')) if sample is not None else '?',
			'align': 'l'}]

NHEADER = len(columns)

def format_overlap(overlap):
	if overlap is None:
		return '?'
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
	if n_points_B>=1.:
		return '{:d}'.format(int(n_points_B))
	else:
		return '{:.1f}'.format(n_points_B)

def format_timing(time):
	if time is None:
		return '?'
	return '{:.0f}'.format(time*10**6)

def format_exponent(results,errors,references):
	if not isinstance(references[0], str) and not isinstance(references[1], str):
		exponents_res = [int('{:e}'.format(results[0]).split('e')[1]),int('{:e}'.format(results[1]).split('e')[1])]
		max_exp_res = max(exponents_res)
		exponents_err = [int('{:e}'.format(errors[0]).split('e')[1]),int('{:e}'.format(errors[1]).split('e')[1])]
		max_exp_err = max(exponents_err)
		max_exp_ltd = max(max_exp_res,max_exp_err)
		exponents_ref = [int('{:e}'.format(references[0]).split('e')[1]),int('{:e}'.format(references[1]).split('e')[1])]
		if references[0] == 0.:
			max_exp_ref = exponents_ref[1]
		elif references[1] == 0.:
			max_exp_ref = exponents_ref[0]
		else:
			max_exp_ref = max(exponents_ref)
		assert(max_exp_ref==max_exp_ltd)
	else:
		exponents_res = [int('{:e}'.format(results[0]).split('e')[1]),int('{:e}'.format(results[1]).split('e')[1])]
		max_exp_res = max(exponents_res)
		exponents_err = [int('{:e}'.format(errors[0]).split('e')[1]),int('{:e}'.format(errors[1]).split('e')[1])]
		max_exp_err = max(exponents_err)
		max_exp_ltd = max([max_exp_res,max_exp_err])
	return '{:e}'.format(10**max_exp_ltd).split('e')[1]

def format_reference(references):
	references_str = []
	if not isinstance(references[0], str) and not isinstance(references[1], str):
		exponents = [int('{:e}'.format(references[0]).split('e')[1]),int('{:e}'.format(references[1]).split('e')[1])]
		if references[0] == 0.:
			max_exp = exponents[1]
		elif references[1] == 0.:
			max_exp = exponents[0]
		else:
			max_exp = max(exponents)
	for phase in [0,1]: 
		reference = references[phase]
		if isinstance(reference, str):
			references_str += [reference]
		else:
			reference_str = r'{: .5f}'.format(reference/10**max_exp) #+ 'e'+'{:e}'.format(10**max_exp).split('e')[1]
			#reference_str = reference_str.replace(r'e',r'$\mathtt{\cdot 10^{')
			#reference_str += '}}$'
			references_str += [r'\texttt{'+reference_str+r'}']
	return references_str

def format_ltd(results,errors):
	results_str = []
	exponents_res = [int('{:e}'.format(results[0]).split('e')[1]),int('{:e}'.format(results[1]).split('e')[1])]
	max_exp_res = max(exponents_res)
	exponents_err = [int('{:e}'.format(errors[0]).split('e')[1]),int('{:e}'.format(errors[1]).split('e')[1])]
	max_exp_err = max(exponents_err)
	max_exp = max([max_exp_err,max_exp_res])
	digits = 5
	for phase in [0,1]:
		res = ufloat(results[phase], errors[phase])
		if res is None:
			return '?'
		#res_str = '{:ue}'.format(res)
		res_str = ('{:.'+str(digits)+'f}').format(res/10**max_exp) #+ 'e'+'{:e}'.format(10**max_exp).split('e')[1]
		if res_str[0] == r'(':
			if res_str[1] != r'-':
				res_str = r'(~' + res_str[1:]
		elif res_str[0] == r'-':
			res_str = r'~'+ res_str
		else:
			res_str = r'~' + r'~'+ res_str
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
			if delta < 1e-3:
				deltas_str += [r'{\color{ForestGreen}'+'{:.0e}'.format(delta)+r'}']
			if delta >= 10.:
				deltas_str += [r'{\color{red}'+'{:.0e}'.format(delta)+r'}']
			else:
				deltas_str += ['{:.3f}'.format(delta)]
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

def set_topology(n_ps_points,data=None):
	topology_string = r'\hline' + '\n'
	for i in range(n_ps_points):
		if i == 0:
			topology_name = str(columns[0]['json_link'](data)).replace('_',r'\_')
			topology_string += r'\multirow{%i}{*}{'%(2*n_ps_points) +topology_name+'}' +'\n'
		else:
			topology_string += r'\cline{2-%i}'%NHEADER +'\n'
		for column in columns[1:]:
			if column['shared']:
				column_str = str(column['json_link'](data))
				topology_string += r'& \multirow{2}{*}{'+column_str+'}'
			else:
				column_str = str(column['json_link'](data)[0])
				topology_string +=  r'& ' + column_str	
		topology_string += r'\\'+'\n'
		for column in columns[1:]:
			if column['shared']:
				topology_string += r'& '
			else:
				column_str = str(column['json_link'](data)[1])
				topology_string +=  r'& ' + column_str	
		topology_string += r'\\'+'\n'
	return topology_string

def set_table(columns,ps_pts_per_topology,data=None):
	table_string = r'\begin{table}[tbp]'+'\n'
	table_string += r'\centering'+'\n'
	table_string += r'\resizebox{\columnwidth}{!}{%'+'\n'
	alignment = ''.join([column['align'] for column in columns])
	table_string += r'\texttt{%'+'\n'
	table_string += r'\begin{tabular}{'+alignment+'}'+'\n'
	table_string += set_header(columns)
	if data is not None:
		for n_ps_points,topo_data in zip(ps_pts_per_topology,data):
			table_string += set_topology(n_ps_points,data=topo_data)
	else:
		for n_ps_points in ps_pts_per_topology:
			table_string += set_topology(n_ps_points,data=None)
	table_string += r'\end{tabular}'+'\n'
	table_string += r'}'+'\n'
	table_string += r'}'+'\n'
	table_string += r'\caption{\label{tab:i} Results}'+'\n'
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

def update_and_extract_data(data,name_maps):
	sorted_data = []
	found = False
	for name_map in name_maps:
		for sample in data:
			if sample['topology'] == name_map['topology_name']:
				sample['PS_point'] = name_map['PS_point']
				sample['graph'] = name_map['graph']
				sorted_data += [sample]
				data.remove(sample)
				found = True
		if not found:
			print(topology_name)
		found = False
	return sorted_data, data


def sort_and_split_data(data):
	topologies_in_order = [	'Pentagon_1s',
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
	counter_5p = 1
	counter_6p = 1
	name_maps = []
	ps_point = None
	graph = None
	for topology_name in topologies_in_order:
		if 'Pentagon' in topology_name:
			ps_point = str(roman.Roman(counter_5p)).lower()
			graph = '1L5P'
			counter_5p += 1
		elif 'Hexagon' in topology_name:
			ps_point = str(roman.Roman(counter_6p)).lower()
			graph = '1L6P'
			counter_6p += 1
		name_maps += [{'topology_name': topology_name, 'PS_point': ps_point, 'graph': graph}]

	explore_1loop_3B, data = update_and_extract_data(data,name_maps)

	topologies_in_order = [	'1L_4P_PS3',
							'1L_4P_PS1',
							'1L_4P_PS2',
							'1L_4P_PS3_massive',
							'1L_4P_PS1_massive',
							'1L_4P_PS2_massive',
							'1L_5P_PS3',
							'1L_5P_PS1',
							'1L_5P_PS2',
							'1L_5P_PS3_massive',
							'1L_5P_PS1_massive',
							'1L_5P_PS2_massive',
							'1L_6P_PS3',
							'1L_6P_PS1',
							'1L_6P_PS2',
							'1L_6P_PS3_massive',
							'1L_6P_PS1_massive',
							'1L_6P_PS2_massive',
							'1L_8P_PS3',
							'1L_8P_PS1',
							'1L_8P_PS2',
							'1L_8P_PS3_massive',
							'1L_8P_PS1_massive',
							'1L_8P_PS2_massive']

	name_maps = []
	ps_point = None
	graph = None
	for topology_name in topologies_in_order:
		if '1L_4P' in topology_name:
			graph = '1L4P'
		elif '1L_5P' in topology_name:
			graph = '1L5P'
		elif '1L_6P' in topology_name:
			graph = '1L6P'
		elif '1L_8P' in topology_name:
			graph = '1L8P'
		if 'PS3' in topology_name:
			ps_point = str(roman.Roman(1))
		elif 'PS1' in topology_name:
			ps_point = str(roman.Roman(2))
		elif 'PS2' in topology_name:
			ps_point = str(roman.Roman(3))
		if 'massive' in topology_name:
			ps_point += r'$^*$'
		name_maps += [{'topology_name': topology_name, 'PS_point': ps_point, 'graph': graph}]
		
	systematic_1loop_3B, data = update_and_extract_data(data,name_maps)


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
							'T4_Quadruple_Box_Weinzierl']

	name_maps = []
	ps_point = None
	graph = None
	for topology_name in topologies_in_order:
		if 'DoubleBox' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L4P.c'
		elif '6P_2L_Weinzierl_A' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.a'
		elif '6P_2L_Weinzierl_B' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.b'
		elif '6P_2L_Weinzierl_C' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.c'
		elif '6P_2L_Weinzierl_D' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.d'
		elif '6P_2L_Weinzierl_E' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.e'
		elif '6P_2L_Weinzierl_F' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L6P.f'
		elif 'TripleBox' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '3L4P'
		elif 'Quadruple_Box' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '4L4P'
		elif 'TM1_top' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L4P.a'
		elif 'TM1_bot' in topology_name:
			ps_point = str(roman.Roman(1)).lower()
			graph = '2L4P.b'
		name_maps += [{'topology_name': topology_name, 'PS_point': ps_point, 'graph': graph}]

	explore_HigherLoop, data = update_and_extract_data(data,name_maps)


	topologies_in_order = [	'2L_4P_Ladder_PS3',
							'2L_4P_Ladder_PS1',
							'2L_4P_Ladder_PS2',
							'2L_4P_Ladder_PS3_massive',
							'2L_4P_Ladder_PS1_massive',
							#'2L_4P_Ladder_PS2_massive',
							'2L_5P_Planar_PS3',
							'2L_5P_Planar_PS1',
							'2L_5P_Planar_PS2',
							'2L_5P_Planar_PS3_massive',
							'2L_5P_Planar_PS1_massive',
							#'2L_5P_Planar_PS2_massive',
							'2L_6P_A_PS3',
							'2L_6P_A_PS1',
							#'2L_6P_A_PS2',
							'2L_6P_A_PS3_massive',
							'2L_6P_A_PS1_massive',
							#'2L_6P_A_PS2_massive',
							'2L_6P_B_PS3',
							'2L_6P_B_PS1',
							'2L_6P_B_PS2',
							'2L_6P_B_PS3_massive',
							'2L_6P_B_PS1_massive',
							#'2L_6P_B_PS2_massive',
							'2L_6P_C_PS3',
							'2L_6P_C_PS1',
							'2L_6P_C_PS2',
							'2L_6P_C_PS3_massive',
							'2L_6P_C_PS1_massive',
							#'2L_6P_C_PS2_massive',
							'2L_6P_D_PS3',
							'2L_6P_D_PS1',
							'2L_6P_D_PS2',
							'2L_6P_D_PS3_massive',
							'2L_6P_D_PS1_massive',
							#'2L_6P_D_PS2_massive',
							'2L_6P_E_PS3',
							'2L_6P_E_PS1',
							'2L_6P_E_PS2',
							'2L_6P_E_PS3_massive',
							'2L_6P_E_PS1_massive',
							#'2L_6P_E_PS2_massive',
							'2L_6P_F_PS3',
							'2L_6P_F_PS1',
							'2L_6P_F_PS2',
							'2L_6P_F_PS3_massive',
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
							#'3L_4P_Ladder_PS2_massive',
							'3L_5P_Planar_PS3',
							'3L_5P_Planar_PS1',
							'3L_5P_Planar_PS2',
							'3L_5P_Planar_PS3_massive',
							'3L_5P_Planar_PS1_massive',
							#'3L_5P_Planar_PS2_massive',
							'FISHNET_2x2_PS3_massive',
							#'FISHNET_2x2_PS1_massive',
							#'FISHNET_2x2_PS2_massive',
							'4L_4P_Ladder_PS3',
							'4L_4P_Ladder_PS1',
							'4L_4P_Ladder_PS2',
							'4L_4P_Ladder_PS3_massive',
							'4L_4P_Ladder_PS1_massive',
							#'4L_4P_Ladder_PS2_massive'
							]
	name_maps = []
	ps_point = None
	graph = None
	for topology_name in topologies_in_order:
		if '2L_4P' in topology_name:
			graph = '2L4P.c'
		elif '2L_5P' in topology_name:
			graph = '2L5P'
		elif '2L_6P_A' in topology_name:
			graph = '2L6P.a'
		elif '2L_6P_B' in topology_name:
			graph = '2L6P.b'
		elif '2L_6P_C' in topology_name:
			graph = '2L6P.c'
		elif '2L_6P_D' in topology_name:
			graph = '2L6P.d'
		elif '2L_6P_E' in topology_name:
			graph = '2L6P.e'
		elif '2L_6P_F' in topology_name:
			graph = '2L6P.f'
		elif '2L_8P' in topology_name:
			graph = '2L8P'
		elif '3L_4P' in topology_name:
			graph = '3L4P'	
		elif 'FISHNET_2x2' in topology_name:
			graph = '4L4P.a'	
		elif '4L_4P' in topology_name:
			graph = '4L4P.b'	
		if 'PS3' in topology_name:
			ps_point = str(roman.Roman(1))
		elif 'PS1' in topology_name:
			ps_point = str(roman.Roman(2))
		elif 'PS2' in topology_name:
			ps_point = str(roman.Roman(3))
		if 'massive' in topology_name:
			ps_point += r'$^*$'
		name_maps += [{'topology_name': topology_name, 'PS_point': ps_point, 'graph': graph}]

	systematic_HigherLoop, data = update_and_extract_data(data,name_maps)

	systematic_HigherLoop_1 = systematic_HigherLoop[:19]
	systematic_HigherLoop_2 = systematic_HigherLoop[19:34]
	systematic_HigherLoop_3 = systematic_HigherLoop[34:]

	assert(data==[])

	return [explore_1loop_3B, systematic_1loop_3B, explore_HigherLoop, systematic_HigherLoop_1,systematic_HigherLoop_2,systematic_HigherLoop_3]




if __name__ == "__main__":

	all_paths = []
	all_paths += [pjoin(file_path,"deformation_paper_results/explore_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS1PS2_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS3_1loop_3B.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/explore_HigherLoop.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS1PS2_HigherLoop.json")]
	all_paths += [pjoin(file_path,"deformation_paper_results/PS3_HigherLoop.json")]

	#DEFAULT_PATH = pjoin(file_path,"deformation_paper_results/PS3_HigherLoop.json")

	#parser = argparse.ArgumentParser(description='Tootl to tabulate json data')
	#parser.add_argument('--path', default=DEFAULT_PATH ,help='Specify the path to the json file.')
	#args = parser.parse_args()

	#PATH = args.path

    #print(set_table(columns,ps_pts_per_topology))
	all_sample_data = []
	for path in all_paths:
		historical_data = get_history(path)
		sample_data = []
		for t, v in historical_data.items():
			samples = v['samples']
			sample_data += samples

		for sample in sample_data:
			sample['accuracy'] = [abs(sample['analytical_result'][i_phase] - sample['result'][i_phase]) / sample['error'][i_phase] for i_phase in [0,1]]
			if sample['accuracy'] == [0.,0.]:
				sample['accuracy'] = ['n/a','n/a']
			#sample['precision'] = [sample['error'][i_phase] * numpy.sqrt(sample['num_samples']) / abs(sample['analytical_result'][i_phase])
			#		if abs(sample['analytical_result'][i_phase]) != 0. else 0 for i_phase in [0,1]]
			#sample['percentage'] = [100.0*(abs(sample['analytical_result'][i_phase]-sample['result'][i_phase]) / abs(sample['analytical_result'][i_phase]))
			#		if abs(sample['analytical_result'][i_phase]) != 0. else 0 for i_phase in [0,1]]
			sample['percentage'] = [100.0*(abs(sample['analytical_result'][i_phase]-sample['result'][i_phase]) / abs(sample['analytical_result'][i_phase]))
					if abs(sample['analytical_result'][i_phase]) != 0. else 'n/a' for i_phase in [0,1]]
			if sample['percentage'] == [0.,0.]:
				sample['percentage'] = ['n/a','n/a']
			if numpy.sqrt(sample['analytical_result'][0]**2+sample['analytical_result'][1]**2) != 0.:
				sample['abs_error'] = 100.*numpy.sqrt(
					(sample['analytical_result'][0]-sample['result'][0])**2
					+ (sample['analytical_result'][1]-sample['result'][1])**2
					)/numpy.sqrt(sample['analytical_result'][0]**2+sample['analytical_result'][1]**2)
			else:
				sample['abs_error'] = 'n/a'
			if sample['abs_error'] == 0.:
				sample['abs_error'] = 'n/a'
			if sample.get('maximal_overlap') is not None:
				sample['overlap_multiplicity'] = [len(overlap) for overlap in sample['maximal_overlap']]
			else:
				sample['overlap_multiplicity'] = None
			if [sample['analytical_result'][0]-sample['result'][0],sample['analytical_result'][1]-sample['result'][1]] == [0.,0.]:
				sample['analytical_result'] = ['n/a','n/a']
		#print(sample_daty.
		#sample_data=sorted(sample_data,key=lambda x: x['topology'])
		new_sample_data = []
		for sample in sample_data:
			if sample.get('n_cuts') is not None:
				new_sample_data += [sample]
		all_sample_data += new_sample_data

	#print(all_sample_data)

	data_sets = sort_and_split_data(all_sample_data)
	for data_set in data_sets:
		print(set_table(columns,[1 for i in range(len(data_set))],data=data_set))











