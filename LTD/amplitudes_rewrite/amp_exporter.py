
import os
from pathlib import Path
import sys
import yaml
import json

if True:
    abspath = os.path.abspath
    pjoin = os.path.join
    src_path = os.path.dirname(os.path.realpath(__file__))    
    sys.path.append(src_path)
    sys.path.append(pjoin(src_path,'..'))
    from amplitudes_rewrite.amplitudes import Amplitude, FormProcessorAmp
    from amplitudes_rewrite.topology_generator import SquaredTopologyGeneratorForAmplitudes



class AmpExporter():

    def __init__(self, out_dir=None, alphaloop_dir=None, process_name=None, amplitude_dict=None, amplitude_options=None, form_options=None):
        
        self.out_dir = out_dir
        self.alphaloop_dir = alphaloop_dir
        self.process_name = process_name
        self.amplitude_dict = amplitude_dict
        self.amplitude_options = amplitude_options
        self.form_options = form_options

    @classmethod
    def from_runcard(AmpExporter, amplitude_runcard_path):
        """Initializes and AmpExporter from an amplitude runcard

        Args:
            amplitude_runcard_path (path): Path to an amplitude runcard
        Returns:
            Class: AmpExporter
        """
        filename, file_extension = os.path.splitext(amplitude_runcard_path)

        with open(amplitude_runcard_path) as f:
            if file_extension == '.yaml':
                runcard_options = yaml.load(f, Loader=yaml.FullLoader)
            elif file_extension == '.json':
                runcard_options = json.load(f)

        amplitude = runcard_options.get('amplitude')
        process_name = amplitude.get('name', 'my_amp')
        amplitude_options = runcard_options.get('amplitude_options', '')
        form_options = runcard_options.get('form_options', '')
        out_dir = os.path.abspath(pjoin(runcard_options.get(
            'directories').get('out_dir', './'), process_name))
        alphaloop_dir = os.path.abspath(
            runcard_options.get('directories').get('alphaloop'))
        amp_exporter = AmpExporter(out_dir=out_dir, amplitude_dict=amplitude, process_name=process_name,
                                   form_options=form_options, amplitude_options=amplitude_options, alphaloop_dir=alphaloop_dir)
        return amp_exporter

    def generate_dir_structure(self, add_out_dir='', mode=''):
        """Generates the output structure needed for and amplitude

        Args:
            add_out_dir (str, optional): Allows for a modification of the output directory. Defaults to ''.
            mode (str, optional): If 'full', then all neccessary sub-directories are build. Defaults to ''.

        Returns:
            out_path (str): path to output directory
        """
        out_dir = pjoin(self.out_dir, add_out_dir)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        if not mode == 'full':
            # create only top_dir
            pass
        else:
            Path(pjoin(out_dir, 'Rust_inputs')).mkdir(
                parents=True, exist_ok=True)
            Path(pjoin(out_dir, 'lib')).mkdir(parents=True, exist_ok=True)
            Path(pjoin(out_dir, 'FORM', 'workspace')).mkdir(
                parents=True, exist_ok=True)
        return out_dir

    def export(self):
        """ Produces the output for an amplitude
        """
        alphaloop_dir = self.alphaloop_dir
        out_dir = self.generate_dir_structure()

        amplitude = Amplitude.import_amplitude(
            self.amplitude_dict, self.amplitude_options)
        amplitude_list = amplitude.split_color(
            out_dir=pjoin(out_dir, 'color_computation'), alphaloop_dir=alphaloop_dir)

        
        for i, amp in enumerate(amplitude_list.amplitudes):
            out_dir = self.generate_dir_structure(
                add_out_dir='color_struc_'+str(i), mode='full')
            topo_gen = SquaredTopologyGeneratorForAmplitudes.from_amplitude(
                amp)
            topo_gen.export_yaml(out_dir+'/Rust_inputs')
            form_processor = FormProcessorAmp(
                amp, alphaloop_dir, form_options=self.form_options, form_wrk_space=pjoin(out_dir, 'FORM', 'workspace'))
            form_processor.generate_output()