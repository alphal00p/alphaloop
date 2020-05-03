#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import os
import stat
import logging
import itertools
from pprint import pformat
from math import fmod
import aloha
import shutil

import aloha.create_aloha as create_aloha

plugin_path = os.path.dirname(os.path.realpath( __file__ ))

from pathlib import Path

from madgraph import MadGraph5Error, InvalidCmd, MG5DIR

import madgraph.core.base_objects as base_objects
import madgraph.iolibs.group_subprocs as group_subprocs

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('alphaLoop.Exporter')

pjoin = os.path.join 

class alphaLoopModelConverter(export_v4.UFO_model_to_mg4):

    def __init__(self, *args, **opts):
        """ initialization of the objects """
        super(alphaLoopModelConverter, self).__init__(*args, **opts)

class alphaLoopExporter(export_v4.ProcessExporterFortranSA):
    check = False
    exporter = 'v4'
    output = 'Template'

    default_opt = {
        'clean': False, 'complex_mass':False,
        'export_format':'standalone', 'mp': False,
        'v5_model': True,
        'sa_symmetry' : False,
        'output_options':{}
    }

    def __init__(self, *args, **opts):
        """ initialization of the objects """
        super(alphaLoopExporter, self).__init__(*args, **opts)

    def copy_template(self, model):
        """Additional actions needed for setup of Template
        """

        #First copy the full template tree if dir_path doesn't exit
        if os.path.isdir(self.dir_path):
            return
        
        logger.info('initialize a new standalone directory: %s' % \
                        os.path.basename(self.dir_path))
        temp_dir = pjoin(self.mgme_dir, 'Template/LO')
        
        # Create the directory structure
        Path(self.dir_path).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'Source')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'Source', 'MODEL')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'Source', 'DHELAS')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'SubProcesses')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'bin')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'bin', 'internal')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'lib')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path, 'Cards')).mkdir(parents=True, exist_ok=True)
        
        # Information at top-level
        #Write version info
        shutil.copy(pjoin(temp_dir, 'TemplateVersion.txt'), self.dir_path)
        try:
            shutil.copy(pjoin(self.mgme_dir, 'MGMEVersion.txt'), self.dir_path)
        except IOError:
            MG5_version = misc.get_pkg_info()
            open(pjoin(self.dir_path, 'MGMEVersion.txt'), 'w').write( \
                "5." + MG5_version['version'])
        
        
        # Add file in SubProcesses
        shutil.copy(pjoin(self.mgme_dir, 'madgraph', 'iolibs', 'template_files', 'makefile_sa_f_sp'), 
                    pjoin(self.dir_path, 'SubProcesses', 'makefileP'))
        
        if self.format == 'standalone':
            shutil.copy(pjoin(self.mgme_dir, 'madgraph', 'iolibs', 'template_files', 'check_sa.f'), 
                    pjoin(self.dir_path, 'SubProcesses', 'check_sa.f'))
                        
        # Add file in Source
        shutil.copy(pjoin(temp_dir, 'Source', 'make_opts'), 
                    pjoin(self.dir_path, 'Source'))        
        # add the makefile 
        filename = pjoin(self.dir_path,'Source','makefile')
        self.write_source_makefile(writers.FileWriter(filename))    

    def finalize(self, matrix_elements, history, mg5options, flaglist, *args, **opts):
        """Finalize Standalone MG4 directory by 
           generation proc_card_mg5.dat
           generate a global makefile
        """
        super(alphaLoopExporter, self).finalize(
            matrix_elements, history, mg5options, flaglist, *args, **opts)

    #===========================================================================
    # process exporter fortran switch between group and not grouped
    #===========================================================================
    def export_processes(self, matrix_elements, fortran_model, *args, **opts):
        """Make the switch between grouped and not grouped output"""
        
        calls = 0
        if isinstance(matrix_elements, group_subprocs.SubProcessGroupList):
            for (group_number, me_group) in enumerate(matrix_elements):
                calls = calls + self.generate_subprocess_directory(\
                                          me_group, fortran_model, group_number)
        else:
            for me_number, me in enumerate(matrix_elements.get_matrix_elements()):
                calls = calls + self.generate_subprocess_directory(\
                                                   me, fortran_model, me_number)    
                        
        return calls

    def convert_model(self, model, wanted_lorentz = [], wanted_couplings = []):
        """ Create a full valid MG4 model from a MG5 model (coming from UFO)"""

        # Make sure aloha is in quadruple precision if needed
        old_aloha_mp=aloha.mp_precision
        aloha.mp_precision=self.opt['mp']

        # create the MODEL
        write_dir=pjoin(self.dir_path, 'Source', 'MODEL')
        model_builder = alphaLoopModelConverter(model, write_dir, self.opt + self.proc_characteristic)
        model_builder.build(wanted_couplings)

        # Backup the loop mode, because it can be changed in what follows.
        old_loop_mode = aloha.loop_mode

        # Create the aloha model or use the existing one (for loop exporters
        # this is useful as the aloha model will be used again in the 
        # LoopHelasMatrixElements generated). We do not save the model generated
        # here if it didn't exist already because it would be a waste of
        # memory for tree level applications since aloha is only needed at the
        # time of creating the aloha fortran subroutines.
        if hasattr(self, 'aloha_model'):
            aloha_model = self.aloha_model
        else:
            aloha_model = create_aloha.AbstractALOHAModel(os.path.basename(model.get('modelpath')))
        aloha_model.add_Lorentz_object(model.get('lorentz'))

        # Compute the subroutines
        if wanted_lorentz:
            aloha_model.compute_subset(wanted_lorentz)
        else:
            aloha_model.compute_all(save=False)

        # Write them out
        write_dir=pjoin(self.dir_path, 'Source', 'DHELAS')
        aloha_model.write(write_dir, 'Fortran')

        # Revert the original aloha loop mode
        aloha.loop_mode = old_loop_mode

        #copy Helas Template
        cp(MG5DIR + '/aloha/template_files/Makefile_F', write_dir+'/makefile')
        if any([any(['L' in tag for tag in d[1]]) for d in wanted_lorentz]):
            cp(MG5DIR + '/aloha/template_files/aloha_functions_loop.f', 
                                                 write_dir+'/aloha_functions.f')
            aloha_model.loop_mode = False
        else:
            cp(MG5DIR + '/aloha/template_files/aloha_functions.f', 
                                                 write_dir+'/aloha_functions.f')
        create_aloha.write_aloha_file_inc(write_dir, '.f', '.o')

        # Make final link in the Process
        self.make_model_symbolic_link()
    
        # Re-establish original aloha mode
        aloha.mp_precision=old_aloha_mp

    def generate_subprocess_directory(self, matrix_element, fortran_model, number, *args, **opts):
        """Generate the Pxxxxx directory for a subprocess in MG4 standalone,
        including the necessary matrix.f and nexternal.inc files"""

        super(alphaLoopExporter, self).generate_subprocess_directory(
            matrix_element, fortran_model, number, *args, **opts)

    #===========================================================================
    # write_source_makefile
    #===========================================================================
    def write_source_makefile(self, writer, *args, **opts):
        super(alphaLoopExporter, self).write_source_makefile(writer, *args, **opts)

    #===========================================================================
    # write_matrix_element_v4
    #===========================================================================
    def write_matrix_element_v4(self, writer, matrix_element, fortran_model,
                                write=True, proc_prefix='',**opts):
        super(alphaLoopExporter, self).write_matrix_element_v4(
            writer, matrix_element, fortran_model, write=write, proc_prefix=proc_prefix,**opts)

         