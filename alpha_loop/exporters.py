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
import progressbar
import aloha
import shutil

import aloha.create_aloha as create_aloha

plugin_path = os.path.dirname(os.path.realpath( __file__ ))

from pathlib import Path

from madgraph import MadGraph5Error, InvalidCmd, MG5DIR

import madgraph.core.base_objects as base_objects
import madgraph.iolibs.group_subprocs as group_subprocs
import madgraph.iolibs.drawing_eps as draw

import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc
import madgraph.iolibs.files as files
from madgraph.iolibs.files import cp, ln, mv

logger = logging.getLogger('alphaLoop.Exporter')

import alpha_loop.LTD_squared as LTD_squared
import alpha_loop.helas_call_writers as aL_helas_call_writers
import alpha_loop.utils as utils

import yaml
from yaml import Loader, Dumper

plugin_src_path = os.path.dirname(os.path.realpath( __file__ ))

pjoin = os.path.join 

class alphaLoopExporterError(MadGraph5Error):
    """ Error for the alphaLoop exporter """
    pass

class alphaLoopModelConverter(export_v4.UFO_model_to_mg4):

    def __init__(self, *args, **opts):
        """ initialization of the objects """
        super(alphaLoopModelConverter, self).__init__(*args, **opts)

class alphaLoopExporter(export_v4.ProcessExporterFortranSA):
    check = False
    exporter = 'v4'
    output = 'Template'

    matrix_template = pjoin(plugin_src_path,'Templates','matrix_standalone.inc')

    default_opt = {
        'clean': False, 'complex_mass':False,
        'export_format':'standalone', 'mp': False,
        'v5_model': True,
        'sa_symmetry' : False,
        'output_options':{}
    }

    def __init__(self, *args, **opts):
        """ initialization of the objects """
        
        self.alphaLoop_options = opts.pop('alphaLoop_options',{})
        
        super(alphaLoopExporter, self).__init__(*args, **opts)
        
        # Force prefixing
        self.cmd_options['prefix']='int'
        # Keep track of all ME exported so as to build the C_bindings dispatcher in finalize
        self.all_processes_exported = []
        # Keep track of all supergraphs generated so as to write the overall LTD squared info at the end of the
        # generation in finalize()
        self.all_super_graphs = []

    #===========================================================================
    # write LTD^2 information destined to the rust backend
    #===========================================================================
    def write_LTD_squared_info(self, matrix_element, proc_number, *args, **opts):
        """ Process the matrix element object so as to extract all necessary information
        for the fortran matrix element to be used by the rust_backend processor so as
        to perform LTD^2 cross-section computations. """

        LTD_square_info_writer = writers.FileWriter('LTD_squared_info.yaml')
        
        # Get a diagram_generation.Amplitude instance corresponding to this HelasMatrixElement
        base_amplitude = matrix_element.get('base_amplitude')
        base_diagrams = base_amplitude.get('diagrams')

        # First build an instance "LT2DiagramLst" from the Helas matrix element object
        LTD2_diagram_list = LTD_squared.LTD2DiagramList(matrix_element)
        all_super_graphs = LTD_squared.SuperGraphList(LTD2_diagram_list, proc_number)
        logger.info("%s%s involves %d supergraphs.%s"%(
            utils.bcolors.GREEN,
            matrix_element.get('processes')[0].nice_string(),
            len(all_super_graphs),
            utils.bcolors.ENDC,
        ))

        # Now write out a yaml file for each of these 
        rust_inputs_path = pjoin(self.dir_path, 'Rust_inputs')
        Path(rust_inputs_path).mkdir(parents=True, exist_ok=True)

        max_super_graphs = 0
        if self.alphaLoop_options['n_rust_inputs_to_generate'] == 0:
            logger.warning("%sUser requested to skip the (slow) generation of rust LTD^2 inputs.%s"%(
                utils.bcolors.RED,utils.bcolors.ENDC
            ))
            return

        if self.alphaLoop_options['n_rust_inputs_to_generate'] > 0:
            logger.warning("%sUser requested to only generate LTD^2 inputs for the first %d inputs.%s"%(
                utils.bcolors.RED,self.alphaLoop_options['n_rust_inputs_to_generate'],utils.bcolors.ENDC
            ))

        base_proc_name = matrix_element.get('processes')[0].shell_string().split('_',1)[1]
        with progressbar.ProgressBar(
            prefix = 'Generating rust inputs for {variables.super_graph_name} : ',
            max_value=(len(all_super_graphs) if self.alphaLoop_options['n_rust_inputs_to_generate']<0 else
                        min(len(all_super_graphs), self.alphaLoop_options['n_rust_inputs_to_generate'])),
            variables = {'super_graph_name' : '%s_%d'%(base_proc_name,1)}) as bar:
            for i_super_graph, super_graph in enumerate(all_super_graphs):
                squared_topology_name = '%s_%d'%(base_proc_name,(i_super_graph+1))
                super_graph.set_name(squared_topology_name)
                bar.update(super_graph_name=super_graph.name)
                bar.update(i_super_graph+1)
                # Keep track of all supergraphs generated so as to write the overall LTD squared info at the end of the
                # generation in finalize()
                self.all_super_graphs.append(super_graph)
                file_path = pjoin(rust_inputs_path,'%s.yaml'%squared_topology_name)
                super_graph.generate_yaml_input_file(
                    file_path,
                    matrix_element.get('processes')[0].get('model'),
                    self.alphaLoop_options
                )
                if (self.alphaLoop_options['n_rust_inputs_to_generate']>0) and (
                    (i_super_graph+1)==self.alphaLoop_options['n_rust_inputs_to_generate']):
                    break

    def write_overall_LTD_squared_info(self, *args, **opts):
        """ Write the overall cross-section yaml input file for the rust_backend to
        be able to generate the whole cross-section at once."""

        rust_inputs_path = pjoin(self.dir_path, 'Rust_inputs')
        # And now finally generate the overall cross section yaml input file.
        base_name = os.path.basename(self.dir_path)
        overall_xsec_yaml_file_path = pjoin(rust_inputs_path,'%s.yaml'%base_name)
        overall_xsec_yaml = {
            'name' : base_name,
            'madgraph_numerators_lib_dir' : pjoin(self.dir_path,'lib'),
            'topologies' : [
                { 
                    'name' : super_graph.name,
                    'multiplicity' : 1
                }
                for super_graph in self.all_super_graphs
            ]
        }
        open(overall_xsec_yaml_file_path,'w').write(yaml.dump(overall_xsec_yaml, Dumper=Dumper, default_flow_style=False))

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
           
        # Add file in Source
        shutil.copy(pjoin(temp_dir, 'Source', 'make_opts'), 
                    pjoin(self.dir_path, 'Source'))        
        # add the makefile 
        filename = pjoin(self.dir_path,'Source','makefile')
        self.write_source_makefile(writers.FileWriter(filename))

    def write_hook_makefile(self):
        """ """
        # Add file in SubProcesses
        shutil.copy(pjoin(plugin_path, 'Templates', 'SubProcesses_makefile'), 
                    pjoin(self.dir_path, 'SubProcesses', 'makefile'))

    def finalize(self, matrix_elements, history, mg5options, flaglist, *args, **opts):
        """Finalize Standalone MG4 directory by 
           generation proc_card_mg5.dat
           generate a global makefile
        """
            
        compiler =  {'fortran': mg5options['fortran_compiler'],
                     'cpp': mg5options['cpp_compiler'],
                     'f2py': mg5options['f2py_compiler']}

        self.compiler_choice(compiler)
        self.make()

        # Write command history as proc_card_mg5
        if history and os.path.isdir(pjoin(self.dir_path, 'Cards')):
            output_file = pjoin(self.dir_path, 'Cards', 'proc_card_mg5.dat')
            history.write(output_file)
        
        export_v4.ProcessExporterFortran.finalize(self, matrix_elements, 
                                             history, mg5options, flaglist)

        open(pjoin(self.dir_path,'__init__.py'),'w')
        open(pjoin(self.dir_path,'SubProcesses','__init__.py'),'w')

        if 'mode' in self.opt and self.opt['mode'] == "reweight":
            #add the module to hande the NLO weight
            files.copytree(pjoin(MG5DIR, 'Template', 'RWGTNLO'),
                          pjoin(self.dir_path, 'Source'))
            files.copytree(pjoin(MG5DIR, 'Template', 'NLO', 'Source', 'PDF'),
                           pjoin(self.dir_path, 'Source', 'PDF'))
            self.write_pdf_opendata()
            
        self.write_f2py_splitter()
        self.write_hook_makefile()

        # Add the C_bindings
        filename = pjoin(self.dir_path, 'SubProcesses','C_bindings.f')
        self.write_c_bindings(writers.FortranWriter(filename))

        # Write the overall cross-section yaml file
        self.write_overall_LTD_squared_info()

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

    def write_c_bindings(self, writer):
        """ Write the fortran dispatcher to all processes callable from C/Rust."""
        
        # Keep track of all ME exported so as to build the C_bindings dispatcher in finalize
        processes = self.all_processes_exported

        replace_dict = {}
        replace_dict['binding_prefix'] = 'C_'
        all_proc_numbers = [proc[1] for proc in processes]
        max_proc_number = max(all_proc_numbers)
        if max_proc_number > 1000:
            raise alphaLoopExporterError("The way the process number map is implemented in"+
                " the alphaLoopExporter is not ideal if your process numbers exceed 1000."+
                " So genealise it by using an actual function instead of an array for the mapping.")
        replace_dict['max_proc_number'] = max_proc_number
        proc_number_to_position_map = []
        for proc_number in range(max_proc_number+1):
            try:
                position = all_proc_numbers.index(proc_number)
                proc_number_to_position_map.append('%d'%position)
            except ValueError:
                proc_number_to_position_map.append('-1')
        replace_dict['proc_number_to_position_map'] = ','.join(proc_number_to_position_map)

        # Now gather all possible combination of number of external legs
        all_externals_combination = set([])
        for me, number in processes:
            all_externals_combination.add(len(me.get('processes')[0].get('legs')))
        # And now generate an input momenta configuration for all of them
        truncated_mom_list = []
        for n_external in sorted(list(all_externals_combination)):
            truncated_mom_list.append("double complex P%d(0:3,%d)"%(n_external,n_external))
        replace_dict['truncated_mom_list'] = '\n'.join(truncated_mom_list)

        replace_dict['max_n_external'] = max(all_externals_combination)

        replace_dict['proc_positions_list'] = ','.join('%d'%pos for pos in range(1,len(processes)+1))

        matrix_element_call_dispatch = []
        for i_proc, (me, proc_number) in enumerate(processes):
            matrix_element_call_dispatch.append('%d CONTINUE'%(i_proc+1))
            # Setup the kinematics now
            matrix_element_call_dispatch.append("DO I=0,3")
            n_external = len(me.get('processes')[0].get('legs'))
            matrix_element_call_dispatch.append("DO J=0,%d"%(n_external-1))
            matrix_element_call_dispatch.append("P%d(I,J+1)=DCMPLX(P(J*8+2*I),P(J*8+2*I+1))"%n_external)
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append(
                'CALL M%d_SMATRIXHEL(P%d, -1, SDL, SDR, ANS)'%(i_proc,n_external))
            matrix_element_call_dispatch.append("GOTO 9999")
        replace_dict['matrix_element_call_dispatch'] = '\n'.join(matrix_element_call_dispatch)

        replace_dict['first_proc_prefix'] = 'M%d_'%(processes[0][1])

        template_path = pjoin(plugin_path,'Templates','C_bindings.inc')
        writer.writelines(open(template_path,'r').read()%replace_dict)
        
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
        aloha.loop_mode = True

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
            cp(pjoin(plugin_src_path,'Templates','aloha_functions.f'), 
                                                 write_dir+'/aloha_functions.f')
        create_aloha.write_aloha_file_inc(write_dir, '.f', '.o')

        # Make final link in the Process
        self.make_model_symbolic_link()
    
        # Re-establish original aloha mode
        aloha.mp_precision=old_aloha_mp

    #===========================================================================
    # generate_subprocess_directory
    #===========================================================================
    def generate_subprocess_directory(self, matrix_element, fortran_model, number):
        """Generate the Pxxxxx directory for a subprocess in MG4 standalone,
        including the necessary matrix.f and nexternal.inc files"""
        
        self.all_processes_exported.append((matrix_element,number))

        # Overwrite fortran model, slightly inefficient but tolerable
        fortran_model = aL_helas_call_writers.alphaLoopHelasCallWriter(
                    matrix_element.get('processes')[0].get('model'),
                    use_physical_gluon_helicity_sum = self.alphaLoop_options['use_physical_gluon_helicity_sum']
                    )


        # Create the directory PN_xx_xxxxx in the specified path
        dirpath = pjoin(self.dir_path, 'SubProcesses', \
                       "P%s" % matrix_element.get('processes')[0].shell_string())

        if self.opt['sa_symmetry']:
            # avoid symmetric output
            for i,proc in enumerate(matrix_element.get('processes')):
                   
                tag = proc.get_tag()     
                legs = proc.get('legs')[:]
                leg0 = proc.get('legs')[0]
                leg1 = proc.get('legs')[1]
                if not leg1.get('state'):
                    proc.get('legs')[0] = leg1
                    proc.get('legs')[1] = leg0
                    flegs = proc.get('legs')[2:]
                    for perm in itertools.permutations(flegs):
                        for i,p in enumerate(perm):
                            proc.get('legs')[i+2] = p
                        dirpath2 =  pjoin(self.dir_path, 'SubProcesses', \
                               "P%s" % proc.shell_string())
                        #restore original order
                        proc.get('legs')[2:] = legs[2:]              
                        if os.path.exists(dirpath2):
                            proc.get('legs')[:] = legs
                            return 0
                proc.get('legs')[:] = legs

        try:
            os.mkdir(dirpath)
        except os.error as error:
            logger.warning(error.strerror + " " + dirpath)

        #try:
        #    os.chdir(dirpath)
        #except os.error:
        #    logger.error('Could not cd to directory %s' % dirpath)
        #    return 0

        logger.info('Creating files in directory %s' % dirpath)

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()

        # Create the matrix.f file and the nexternal.inc file
        if self.opt['export_format']=='standalone_msP':
            filename = pjoin(dirpath, 'matrix_prod.f')
        else:
            filename = pjoin(dirpath, 'matrix.f')
        
        proc_prefix = ''
        if 'prefix' in self.cmd_options:
            if self.cmd_options['prefix'] == 'int':
                proc_prefix = 'M%s_' % number
            elif self.cmd_options['prefix'] == 'proc':
                proc_prefix = matrix_element.get('processes')[0].shell_string().split('_',1)[1]
            else:
                raise Exception('--prefix options supports only \'int\' and \'proc\'')
            for proc in matrix_element.get('processes'):
                ids = [l.get('id') for l in proc.get('legs_with_decays')]
                self.prefix_info[tuple(ids)] = [proc_prefix, proc.get_tag()]

        open(pjoin(dirpath, 'check_sa.f'),'w').write(
            open(pjoin(plugin_src_path,'Templates','check_sa.f'),'r').read().replace('PROC_PREFIX_',proc_prefix)
        )

        calls = self.write_matrix_element_v4(
            writers.FortranWriter(filename),
            matrix_element,
            fortran_model,
            proc_prefix=proc_prefix)

        if self.opt['export_format'] == 'standalone_msP':
            filename =  pjoin(dirpath,'configs_production.inc')
            mapconfigs, s_and_t_channels = self.write_configs_file(\
                writers.FortranWriter(filename),
                matrix_element)

            filename =  pjoin(dirpath,'props_production.inc')
            self.write_props_file(writers.FortranWriter(filename),
                             matrix_element,
                             s_and_t_channels)

            filename =  pjoin(dirpath,'nexternal_prod.inc')
            self.write_nexternal_madspin(writers.FortranWriter(filename),
                             nexternal, ninitial)

        if self.opt['export_format']=='standalone_msF':
            filename = pjoin(dirpath, 'helamp.inc')
            ncomb=self.get_helicity_combinations(matrix_element)
            self.write_helamp_madspin(writers.FortranWriter(filename),
                             ncomb)
            
        filename = pjoin(dirpath, 'nexternal.inc')
        self.write_nexternal_file(writers.FortranWriter(filename),
                             nexternal, ninitial)

        filename = pjoin(dirpath, 'pmass.inc')
        self.write_pmass_file(writers.FortranWriter(filename),
                         matrix_element)

        filename = pjoin(dirpath, 'ngraphs.inc')
        self.write_ngraphs_file(writers.FortranWriter(filename),
                           len(matrix_element.get_all_amplitudes()))

        # Generate diagrams
        if not 'noeps' in self.opt['output_options'] or self.opt['output_options']['noeps'] != 'True':
            filename = pjoin(dirpath, "matrix.ps")
            plot = draw.MultiEpsDiagramDrawer(matrix_element.get('base_amplitude').\
                                                 get('diagrams'),
                                              filename,
                                              model=matrix_element.get('processes')[0].\
                                                 get('model'),
                                              amplitude=True)
            logger.info("Generating Feynman diagrams for " + \
                         matrix_element.get('processes')[0].nice_string())
            plot.draw()

        linkfiles = ['coupl.inc']

        for file in linkfiles:
            ln('../%s' % file, cwd=dirpath)
        ln('../makefileP', name='makefile', cwd=dirpath)
        # Return to original PWD
        #os.chdir(cwd)

        # Generate LTD squared info
        self.write_LTD_squared_info(matrix_element, number)

        if not calls:
            calls = 0

        return calls

    def get_helicity_states(self, model, pdg, state, allow_reverse=True):
        """ Modified list of helicities so as to have full control over the polarisation sum
        of our numerator."""

        particle = model.get('particle_dict')[pdg]
        if particle.get('spin')>3:
            raise alphaLoopExporterError("alphaLoop only support particles with spin 0, 1/2 or 1.")
        # Only modify the helicity of *final state* particles that are spin-1 vectors.
        if state==True and particle.get('spin') in [2,3]:
            return [0, 1, 2, 3]
        return particle.get_helicity_states(allow_reverse)

    def get_helicity_matrix(self, matrix_element, allow_reverse=True):
        """Gives the helicity matrix for external wavefunctions"""

        if not matrix_element.get('processes'):
            return None

        process = matrix_element.get('processes')[0]
        model = process.get('model')
        hel_per_part = [ self.get_helicity_states(model, wf.get('pdg_code'), wf.get('leg_state'),
            allow_reverse=allow_reverse) for wf in matrix_element.get_external_wavefunctions()]
        return itertools.product(*hel_per_part)

    def get_helicity_combinations(self, matrix_element):
        return len(list(self.get_helicity_matrix(matrix_element)))

    def get_helicity_lines(self, matrix_element,array_name='NHEL'):
        """Return the Helicity matrix definition lines for this matrix element"""

        helicity_line_list = []
        i = 0
        for helicities in self.get_helicity_matrix(matrix_element):
            i = i + 1
            int_list = [i, len(helicities)]
            int_list.extend(helicities)
            helicity_line_list.append(\
                ("DATA ("+array_name+"(I,%4r),I=1,%d) /" + \
                 ",".join(['%2r'] * len(helicities)) + "/") % tuple(int_list))

        return "\n".join(helicity_line_list)

    #===========================================================================
    # write_matrix_element_v4
    #===========================================================================
    def write_matrix_element_v4(self, writer, matrix_element, fortran_model,
                                write=True, proc_prefix=''):
        """Export a matrix element to a matrix.f file in MG4 standalone format
        if write is on False, just return the replace_dict and not write anything."""

        if not matrix_element.get('processes') or \
               not matrix_element.get('diagrams'):
            return 0

        if writer:
            if not isinstance(writer, writers.FortranWriter):
                raise writers.FortranWriter.FortranWriterError(\
                "writer not FortranWriter but %s" % type(writer))
            # Set lowercase/uppercase Fortran code
            writers.FortranWriter.downcase = False

        if 'sa_symmetry' not in self.opt:
            self.opt['sa_symmetry']=False

        # The proc_id is for MadEvent grouping which is never used in SA.
        replace_dict = {'global_variable':'', 'amp2_lines':'',
                                       'proc_prefix':proc_prefix, 'proc_id':''}

        # Extract helas calls
        helas_calls = fortran_model.get_matrix_element_calls(\
                    matrix_element)

        replace_dict['helas_calls'] = "\n".join(helas_calls)

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(matrix_element)
        replace_dict['process_lines'] = process_lines

        # List of diagram indices
        replace_dict['diagram_indices'] = ','.join(
            '%d'%diag.get('number') for diag in matrix_element.get('diagrams'))

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
        replace_dict['nexternal'] = nexternal
        replace_dict['nincoming'] = ninitial

        # Extract ncomb
        ncomb = self.get_helicity_combinations(matrix_element)
        replace_dict['ncomb'] = ncomb

        # Extract helicity lines
        helicity_lines = self.get_helicity_lines(matrix_element)
        replace_dict['helicity_lines'] = helicity_lines

        # Extract overall denominator
        # Averaging initial state color, spin, and identical FS particles
        replace_dict['den_factor_line'] = self.get_den_factor_line(matrix_element)

        # Extract ngraphs
        ngraphs = matrix_element.get_number_of_amplitudes()
        replace_dict['ngraphs'] = ngraphs

        # Extract nwavefuncs
        nwavefuncs = matrix_element.get_number_of_wavefunctions()
        replace_dict['nwavefuncs'] = nwavefuncs

        # Extract ncolor
        ncolor = max(1, len(matrix_element.get('color_basis')))
        replace_dict['ncolor'] = ncolor

        replace_dict['hel_avg_factor'] = matrix_element.get_hel_avg_factor()
        replace_dict['beamone_helavgfactor'], replace_dict['beamtwo_helavgfactor'] =\
                                       matrix_element.get_beams_hel_avg_factor()

        # Extract color data lines
        color_data_lines = self.get_color_data_lines(matrix_element)
        replace_dict['color_data_lines'] = "\n".join(color_data_lines)

        if self.opt['export_format']=='standalone_msP':
        # For MadSpin need to return the AMP2
            amp2_lines = self.get_amp2_lines(matrix_element, [] )
            replace_dict['amp2_lines'] = '\n'.join(amp2_lines)
            replace_dict['global_variable'] = \
         "       Double Precision amp2(NGRAPHS)\n       common/to_amps/  amp2\n"

        # JAMP definition, depends on the number of independent split orders
        split_orders=matrix_element.get('processes')[0].get('split_orders')

        if len(split_orders)==0:
            replace_dict['nSplitOrders']=''
            # Extract JAMP lines
            jamp_lines = self.get_JAMP_lines(matrix_element,JAMP_format="JAMP(K,%s)")
            # Consider the output of a dummy order 'ALL_ORDERS' for which we
            # set all amplitude order to weight 1 and only one squared order
            # contribution which is of course ALL_ORDERS=2.
            squared_orders = [(2,),]
            amp_orders = [((1,),tuple(range(1,ngraphs+1)))]
            replace_dict['chosen_so_configs'] = '.TRUE.'
            replace_dict['nSqAmpSplitOrders']=1
            replace_dict['split_order_str_list']=''
        else:
            squared_orders, amp_orders = matrix_element.get_split_orders_mapping()
            replace_dict['nAmpSplitOrders']=len(amp_orders)
            replace_dict['nSqAmpSplitOrders']=len(squared_orders)
            replace_dict['nSplitOrders']=len(split_orders)
            replace_dict['split_order_str_list']=str(split_orders)
            amp_so = self.get_split_orders_lines(
                    [amp_order[0] for amp_order in amp_orders],'AMPSPLITORDERS')
            sqamp_so = self.get_split_orders_lines(squared_orders,'SQSPLITORDERS')
            replace_dict['ampsplitorders']='\n'.join(amp_so)
            replace_dict['sqsplitorders']='\n'.join(sqamp_so)           
            jamp_lines = self.get_JAMP_lines_split_order(\
                       matrix_element,amp_orders,
                       JAMP_format="JAMP(K,%s)",
                       split_order_names=split_orders)
            
            # Now setup the array specifying what squared split order is chosen
            replace_dict['chosen_so_configs']=self.set_chosen_SO_index(
                              matrix_element.get('processes')[0],squared_orders)
            
            # For convenience we also write the driver check_sa_splitOrders.f
            # that explicitely writes out the contribution from each squared order.
            # The original driver still works and is compiled with 'make' while
            # the splitOrders one is compiled with 'make check_sa_born_splitOrders'
            check_sa_writer=writers.FortranWriter('check_sa_born_splitOrders.f')
            self.write_check_sa_splitOrders(squared_orders,split_orders,
              nexternal,ninitial,proc_prefix,check_sa_writer)

        if write:
            writers.FortranWriter('nsqso_born.inc').writelines(
                """INTEGER NSQSO_BORN
                   PARAMETER (NSQSO_BORN=%d)"""%replace_dict['nSqAmpSplitOrders'])

        # We must make sure to replace IMAG by its appropriate sign since we are doing
        # "manual" complex-conjugation.
        replace_dict['jamp_lines'] = ('\n'.join(jamp_lines)).upper().replace('IMAG1','IMAG1(K)')

        matrix_template = self.matrix_template

        if len(split_orders)>0:
            raise alphaLoopExporterError("Split orders currently not supported by alphaLoop.")
            if self.opt['export_format'] in ['standalone_msP', 'standalone_msF']:
                logger.debug("Warning: The export format %s is not "+\
                  " available for individual ME evaluation of given coupl. orders."+\
                  " Only the total ME will be computed.", self.opt['export_format'])
            else:
                matrix_template = "matrix_standalone_splitOrders_v4.inc"

        replace_dict['template_file'] = pjoin(MG5DIR,'madgraph', 'iolibs', 'template_files', matrix_template)
        replace_dict['template_file2'] = pjoin(MG5DIR,'madgraph', \
                                   'iolibs/template_files/split_orders_helping_functions.inc')

        if write and writer:
            path = replace_dict['template_file']
            content = open(path).read()
            content = content % replace_dict
            # Write the file
            writer.writelines(content)
            # Add the helper functions.
            if len(split_orders)>0:
                content = '\n' + open(replace_dict['template_file2'])\
                                   .read()%replace_dict
                writer.writelines(content)
            return len([call for call in helas_calls if call.find('#') != 0])
        else:
            replace_dict['return_value'] = len([call for call in helas_calls if call.find('#') != 0])
            return replace_dict # for subclass update

    #===========================================================================
    # write_source_makefile
    #===========================================================================
    def write_source_makefile(self, writer, *args, **opts):
        super(alphaLoopExporter, self).write_source_makefile(writer, *args, **opts)

         