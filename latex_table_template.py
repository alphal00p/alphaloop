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
			'json_link': lambda sample: sample['topology'] if sample is not None else '?'}]
columns += [{'header': r'PS point',
			'shared': True,
			'json_link': lambda sample: '?' if sample is not None else '?'}] #sample['ps_name']}]
columns += [{'header': r'$N_{\text{e}}$',
			'shared': True,
			'json_link': lambda sample: '?' if sample is not None else '?'}] #sample['n_ellipsoids']}]
columns += [{'header': r'$N_{\text{s}}$',
			'shared': True,
			'json_link': lambda sample: '?' if sample is not None else '?'}] #sample['n_sources']}]
columns += [{'header': r'$L_{\text{max}}$',
			'shared': True,
			'json_link': lambda sample: '?' if sample is not None else '?'}] #sample['n_sources']}]
columns += [{'header': r'$N_{\text{points}}$',
			'shared': True,
			'json_link': lambda sample: sample['num_samples'] if sample is not None else '?'}]
columns += [{'header': r'$\frac{t}{p} [\mu \text{s}]$',
			'shared': True,
			'json_link': lambda sample: '?'}]
columns += [{'header': r'Phase',
			'shared': False,
			'json_link': [lambda sample: r'$\Re$',
						lambda sample: r'$\Im$']}]
columns += [{'header': r'Reference result',
			'shared': False,
			'json_link': [lambda sample: '{:.5e}'.format(sample['analytical_result'][0]) if sample is not None else '?',
							lambda sample: '{:.5e}'.format(sample['analytical_result'][1]) if sample is not None else '?']}]
columns += [{'header': r'numerical LTD',
			'shared': False,
			'json_link': [lambda sample: ufloat(sample['result'][0], sample['error'][0]) if sample is not None else '?',
						lambda sample: ufloat(sample['result'][1], sample['error'][1]) if sample is not None else '?']
			}]
columns += [{'header': r'$\Delta[\sigma]$',
			'shared': False,
			'json_link': [lambda sample: '{:.2f}'.format(sample['accuracy'][0]) if sample is not None else '?',
							lambda sample: '{:.2f}'.format(sample['accuracy'][1]) if sample is not None else '?']}]
columns += [{'header': r'$\Delta[\%]$',
			'shared': False,
			'json_link': [lambda sample: '{:.2f}'.format(sample['percentage'][0]) if sample is not None else '?',
							lambda sample: '{:.2f}'.format(sample['percentage'][1]) if sample is not None else '?']}]
columns += [{'header': r'$\Delta[\%] |\cdot|$',
			'shared': True,
			'json_link': lambda sample: '{:.2e}'.format(sample['abs_error']) if sample is not None else '?'}]

NHEADER = len(columns)

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
	"""
		for column in columns[1:5]:
			name = str(column['json_link'](data))
			topology_string += r'& \multirow{2}{*}{'+name+'}'
		topology_string += '\n'
		topology_string += r'& $\Re$'
		for column in columns[6:]:
			name = str(column['json_link'][0](data))
			topology_string +=  r' & ' + name 
		topology_string += r'\\'+'\n'
		topology_string += r'& & & &' + '\n'
		topology_string += r'& $\Im$'
		for column in columns[6:]:
			name = str(column['json_link'][1](data))
			topology_string += r' & ' + name
		topology_string += r'\\'+'\n'
	return topology_string
	"""

def set_table(columns,ps_pts_per_topology,data=None):
	table_string = r'\begin{table}[tbp]'+'\n'
	table_string += r'\centering'+'\n'
	table_string += r'\begin{tabular}{'+'c'*NHEADER+'}'+'\n'
	table_string += set_header(columns)
	if data is not None:
		for n_ps_points,topo_data in zip(ps_pts_per_topology,data):
			table_string += set_topology(n_ps_points,data=topo_data)
	else:
		for n_ps_points in ps_pts_per_topology:
			table_string += set_topology(n_ps_points,data=None)
	table_string += r'\end{tabular}'+'\n'
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
		sample['percentage'] = [100.0*(abs(sample['error'][i_phase]/sample['result'][i_phase])) for i_phase in [0,1]]
		sample['abs_error'] = numpy.sqrt(sample['error'][0]**2 + sample['error'][1]**2)/numpy.sqrt(sample['result'][0]**2+sample['result'][1]**2)
	#print(sample_daty.
	sample_data=sorted(sample_data,key=lambda x: x['topology'])
	#print(sorted_data)
	print(set_table(columns,[1 for i in range(len(sample_data))],data=sample_data))











