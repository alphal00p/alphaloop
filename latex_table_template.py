#generate Latex table template

import json
import os
from uncertainties import ufloat
import numpy
import argparse

file_path = os.path.dirname(os.path.realpath( __file__ ))
#file_path = '/Users/Dario/Desktop/PhD/Code/pynloop/'
pjoin = os.path.join

ps_pts_per_topology = [2,3,2,2,1] #number of ps points per topology
columns = []
columns += [{'header': r'Topology',
			'shared': True,
			'json_link': lambda sample: sample.get('topology') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'PS',
			'shared': True,
			'json_link': lambda sample: '?' if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{c}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_cuts') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{e}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_unique_existing_E_surface') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{s}}}$',
			'shared': True,
			'json_link': lambda sample: sample.get('n_sources') if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{L_{\text{max}}}$',
			'shared': True,
			'json_link': lambda sample: format_overlap(sample.get('maximal_overlap')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{N_{\text{p}} \ [10^6]}$',
			'shared': True,
			'json_link': lambda sample: format_n_points(sample.get('num_samples')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\sfrac{t}{p} \ [\mu \text{s}]}$',
			'shared': True,
			'json_link': lambda sample: format_timing(sample.get('t_per_ps_point_in_s')) if sample is not None else '?',
			'align': 'l'}]
columns += [{'header': r'Phase',
			'shared': False,
			'json_link': [lambda sample: r'$\Re$',
						lambda sample: r'$\Im$'],
			'align': 'l'}]
columns += [{'header': r'Reference',
			'shared': False,
			'json_link': [lambda sample: format_reference(sample.get('analytical_result')[0]) if sample is not None else '?',
							lambda sample: format_reference(sample['analytical_result'][1]) if sample is not None else '?'],
			'align': 'l'}]
columns += [{'header': r'Numerical LTD',
			'shared': False,
			'json_link': [lambda sample: format_ltd(ufloat(sample.get('result')[0], sample.get('error')[0])) if sample is not None else '?',
						lambda sample: format_ltd(ufloat(sample.get('result')[1], sample.get('error')[1])) if sample is not None else '?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\sigma]}$',
			'shared': False,
			'json_link': [lambda sample: '{:.2f}'.format(sample.get('accuracy')[0]) if sample is not None else '?',
							lambda sample: '{:.2f}'.format(sample.get('accuracy')[1]) if sample is not None else '?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%]}$',
			'shared': False,
			'json_link': [lambda sample: '{:.2f}'.format(sample.get('percentage')[0]) if sample is not None else '?',
							lambda sample: '{:.2f}'.format(sample.get('percentage')[1]) if sample is not None else '?'],
			'align': 'l'}]
columns += [{'header': r'$\mathtt{\Delta \ [\%] |\cdot|}$',
			'shared': True,
			'json_link': lambda sample: '{:.2f}'.format(sample.get('abs_error')) if sample is not None else '?',
			'align': 'l'}]

NHEADER = len(columns)

def format_overlap(overlap):
	overlap_str = str(overlap)
	overlap_str = overlap_str.replace(' ','')
	half = 40.
	center_index = 0
	distance = 1000
	for i,s in enumerate(overlap_str[:-1]):
		if s == ']' and overlap_str[i+1]==',':
			if abs(half-(i+1)) < distance:
				distance = abs(half-(i+1))
				center_index = i+1
	if center_index == 0 or len(overlap_str) < half: 
		overlap_str = overlap_str
	else:
		overlap_str = r'&\mathtt{'+overlap_str[:center_index+1]+r'}\\&\mathtt{'+overlap_str[center_index+1:]+r'}'
	overlap_str = overlap_str.replace('[[','(')
	overlap_str = overlap_str.replace('[','(')
	overlap_str = overlap_str.replace(']]',')')
	overlap_str = overlap_str.replace(']',')')
	return r'$\subalign{'+overlap_str+'}$'

def format_n_points(n_points):
	return '{:d}'.format(int(n_points/10**6))

def format_timing(time):
	return '{:.0f}'.format(time*10**6)

def format_reference(reference):
	reference_str = r'{: .5e}'.format(reference)
	reference_str = reference_str.replace(r'e',r'$\mathtt{\cdot 10^{')
	reference_str += '}}$'
	return r'\texttt{'+reference_str+r'}'

def format_ltd(res):
	res_str = '{:ue}'.format(res)
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
	res_str = res_str.replace(r'+/-',r'~'*(6-counter)+r'+/-~')
	if ')e' in res_str:
		res_str = res_str.replace(r')e',r'~'*(5-counter)+r')'+r'$\mathtt{\cdot 10^{')
		res_str += '}}$'
	else:
		res_str += r'~'*(5-counter)+r')'+r'$\mathtt{\phantom{{}\cdot{}10^{-10}}}$'
	return r'\texttt{'+res_str+r'}'


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
				column_str = str(column['json_link'][0](data))
				topology_string +=  r'& ' + column_str	
		topology_string += r'\\'+'\n'
		for column in columns[1:]:
			if column['shared']:
				topology_string += r'& '
			else:
				column_str = str(column['json_link'][1](data))
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


if __name__ == "__main__":

	DEFAULT_PATH = pjoin(file_path,"deformation_paper_results/explore_1loop_3B.json")

	parser = argparse.ArgumentParser(description='Tootl to tabulate json data')
	parser.add_argument('--path', default=DEFAULT_PATH ,help='Specify the path to the json file.')
	args = parser.parse_args()

	PATH = args.path

    #print(set_table(columns,ps_pts_per_topology))

	historical_data = get_history(PATH)
	sample_data = []
	for t, v in historical_data.items():
		samples = v['samples']
		sample_data += samples

	for sample in sample_data:
		sample['accuracy'] = [abs(sample['analytical_result'][i_phase] - sample['result'][i_phase]) / sample['error'][i_phase] for i_phase in [0,1]]
		sample['precision'] = [sample['error'][i_phase] * numpy.sqrt(sample['num_samples']) / abs(sample['analytical_result'][i_phase])
				if abs(sample['analytical_result'][i_phase]) != 0. else 0 for i_phase in [0,1]]
		#sample['percentage'] = [100.0*(abs(sample['analytical_result'][i_phase]-sample['result'][i_phase]) / abs(sample['analytical_result'][i_phase]))
		#		if abs(sample['analytical_result'][i_phase]) != 0. else 0 for i_phase in [0,1]]
		sample['percentage'] = [100.*(abs(sample['error'][i_phase]/sample['result'][i_phase])) for i_phase in [0,1]]
		sample['abs_error'] = 100.*numpy.sqrt(sample['error'][0]**2 + sample['error'][1]**2)/numpy.sqrt(sample['result'][0]**2+sample['result'][1]**2)
	#print(sample_daty.
	sample_data=sorted(sample_data,key=lambda x: x['topology'])
	for sample in sample_data:
		if sample.get('n_cuts') is None:
			sample_data.remove(sample)
	#print(sample_data[0])

	print(set_table(columns,[1 for i in range(len(sample_data))],data=sample_data))











