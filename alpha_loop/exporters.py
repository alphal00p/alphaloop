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
import copy
import subprocess
import pickle
import time

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
import madgraph.interface.common_run_interface as common_run_interface
import models.check_param_card as check_param_card
import models.model_reader as model_reader

logger = logging.getLogger('alphaLoop.Exporter')

import alpha_loop.LTD_squared as LTD_squared
import alpha_loop.helas_call_writers as aL_helas_call_writers
import alpha_loop.utils as utils
import alpha_loop.FORM_processing as FORM_processing
import alpha_loop.amplitudes as AMP_TOOLS
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

    def copy_standard_file(self):
        """ Specialize this function so as to overwrite the makefile of the MODEL 
        with our customized one."""
        super(alphaLoopModelConverter,self).copy_standard_file()
        shutil.copy(pjoin(plugin_path, 'Templates', 'model_makefile'), 
                    pjoin(self.dir_path, 'makefile'))

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
        self.MG5aMC_options = opts.pop('MG5aMC_options',{})

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
        
        # Get a diagram_generation.Amplitude instance corresponding to this HelasMatrixElement
        process = matrix_element.get('processes')[0]
        # This process is meant to accomodate self-energies if any of its external legs contain
        # a self-energy "anchor" particle
        is_self_energy = any(abs(l.get('id')) >= LTD_squared.self_energy_global_id_offset for l in process.get('legs'))
        # Choose the appropriate class for the kind of super-graph handled here.
        # We could use a class factory, but that would be overkill for the present case.
        LTD2_diagram_list_class = (LTD_squared.LTD2DiagramList if not is_self_energy 
                                                    else LTD_squared.SelfEnergyLTD2DiagramList)
        super_graph_list_class = (LTD_squared.SuperGraphList if not is_self_energy 
                                                    else LTD_squared.SelfEnergySuperGraphList)

        base_amplitude = matrix_element.get('base_amplitude')
        base_diagrams = base_amplitude.get('diagrams')

        # First build an instance "LT2DiagramLst" from the Helas matrix element object
        LTD2_diagram_list = LTD2_diagram_list_class(matrix_element, self.alphaLoop_options)
        all_super_graphs = super_graph_list_class(LTD2_diagram_list, proc_number,self.alphaLoop_options['supergraphs_selection'])
        logger.info("%s%s involves %d %ssupergraphs.%s"%(
            utils.bcolors.GREEN,
            matrix_element.get('processes')[0].nice_string(),
            len(all_super_graphs),
            '' if not is_self_energy else 'self-energy ',
            utils.bcolors.ENDC,
        ))

        # Keep track of all supergraphs generated so as to write the overall LTD squared info at the end of the
        # generation in finalize()
        self.all_super_graphs.append((all_super_graphs,matrix_element))

    def write_overall_LTD_squared_info(self, *args, **opts):
        """ Write the overall cross-section yaml input file for the rust_backend to
        be able to generate the whole cross-section at once."""

        if len(self.all_super_graphs)==0:
            raise alphaLoopExporterError("No supergraph generated.")

        representative_process = self.all_super_graphs[0][1].get('processes')[0]
        model = representative_process.get('model')
        computed_model = model_reader.ModelReader(model)
        computed_model.set_parameters_and_couplings(pjoin(self.dir_path,'Cards','param_card.dat'))      

        # Now write out a yaml file for each of these 
        rust_inputs_path = pjoin(self.dir_path, 'Rust_inputs')
        Path(rust_inputs_path).mkdir(parents=True, exist_ok=True)

        # Compute final state symmetry factor
        for super_graph_list, ME in self.all_super_graphs:
            final_ids = [l.get('id') for l in ME.get('processes')[0].get('legs') if l.get('state') is True]
            fs_symm_factor = 1
            for pdg in set(final_ids):
                fs_symm_factor *= final_ids.count(pdg)
            for sg in super_graph_list:
                sg.fs_symm_factor = fs_symm_factor
                
        # Filter once again for isomorphism across all supergraphs generated
        if self.alphaLoop_options['apply_graph_isomorphisms']:
            overall_basis = []

            def find_isomorphic_partner(sg):
                for partner_candidate in overall_basis:
                    if sg.is_isomorphic_to(partner_candidate):
                        return partner_candidate
                return None

            n_tot_supergraphs = sum(len(sg[0]) for sg in self.all_super_graphs)
            logger.info("Detecting isomorphisms between all %d 'unique' super-graphs generated across all processes..."%n_tot_supergraphs)
            i_graph = 0
            with progressbar.ProgressBar(
                prefix = 'Filtering supergraphs across processes ({variables.n_unique_super_graphs_sofar} unique found so far) : ',
                max_value=n_tot_supergraphs,
                variables = {'n_unique_super_graphs_sofar' : '0'}
                ) as bar:
                for super_graph_list, ME in self.all_super_graphs:
                    new_selection_for_this_process = []
                    for super_graph in list(super_graph_list):
                        i_graph += 1
                        partner_sg = find_isomorphic_partner(super_graph)
                        if partner_sg is None:
                            overall_basis.append(super_graph)
                            new_selection_for_this_process.append(super_graph)
                        else:
                            partner_sg.MG_LO_cuts_corresponding_to_this_supergraph.extend(
                                super_graph.MG_LO_cuts_corresponding_to_this_supergraph)
                        bar.update(n_unique_super_graphs_sofar='%d'%len(overall_basis))
                        bar.update(i_graph)
                    super_graph_list[:] = new_selection_for_this_process
            logger.info("%sA grand total of %d unique super-graphs have been generated across all processes.%s"%(
                utils.bcolors.GREEN,len(overall_basis),utils.bcolors.ENDC
            ))

        # Now output the Rust inputs for the remaining supergraphs.
        is_LO = (not self.alphaLoop_options['perturbative_orders']) or sum(self.alphaLoop_options['perturbative_orders'].values())==0
        # Adjust below to `and False` in order to disable debugging.
        do_debug = is_LO
        if do_debug:
            if len(self.all_super_graphs)==1:
                log_func = logger.critical
                err_text = 'ERROR'
                color = utils.bcolors.RED
            else:
                log_func = logger.warning
                err_text = 'WARNING'
                color = utils.bcolors.BLUE
        else:
            log_func = lambda el: None
            err_text = 'WARNING'
            color = utils.bcolors.BLUE

        _HAS_WARNED_ABOUT_SYMMETRY_FACTORS = False 
        FORM_id=0
        for all_super_graphs, matrix_element in self.all_super_graphs:

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
            n_tot_cutkosky_cuts = 0
            with progressbar.ProgressBar(
                prefix = 'Generating rust inputs (from MadGraph) for {variables.super_graph_name} : ',
                max_value=(len(all_super_graphs) if self.alphaLoop_options['n_rust_inputs_to_generate']<0 else
                            min(len(all_super_graphs), self.alphaLoop_options['n_rust_inputs_to_generate'])),
                variables = {'super_graph_name' : '%s_%d'%(base_proc_name,1)}) as bar:
                for i_super_graph, super_graph in enumerate(all_super_graphs):
                    squared_topology_name = '%s_%d'%(base_proc_name,(i_super_graph+1))
                    super_graph.set_name(squared_topology_name)
                    bar.update(super_graph_name=super_graph.name)
                    bar.update(i_super_graph+1)
                    file_path = pjoin(rust_inputs_path,'%s.yaml'%squared_topology_name)
                    if self.alphaLoop_options['FORM_processing_output_format']:
                        FORM_id_to_supply = FORM_id
                    else:
                        FORM_id_to_supply = None
                    try:
                        super_graph.generate_yaml_input_file(
                            file_path,
                            computed_model,
                            self.alphaLoop_options,
                            FORM_id=FORM_id_to_supply
                        )
                    except Exception as e:
                        logger.exception("%sRust input yaml file generation from MadGraph for supergraph '%s' failed. Proceeding nonetheless.%s"%(
                                        utils.bcolors.RED, super_graph.get_name(), utils.bcolors.ENDC))
                        super_graph.cutkosky_cuts_generated = []             

                    if len(super_graph.cutkosky_cuts_generated)==0:
                        log_func("%s%s: Supergraph '%s' yielded no Cutkosky cuts. That should not happen at LO.%s"%(
                                color,err_text,
                                squared_topology_name,
                                utils.bcolors.ENDC
                            ) )
                        super_graph.symmetry_factor = 1
                    else:
                        fs_symm_factor = super_graph.fs_symm_factor
                        if is_LO:
                            if len(super_graph.MG_LO_cuts_corresponding_to_this_supergraph)%len(super_graph.cutkosky_cuts_generated)!=0:
                                log_func("%s%s: Incorrect LO cutkosky cut generation. %d not divisible by %d.%s"%(
                                    color,err_text,
                                    len(super_graph.MG_LO_cuts_corresponding_to_this_supergraph),
                                    len(super_graph.cutkosky_cuts_generated),
                                    utils.bcolors.ENDC
                                ) )
                            reconstructed_symmetry_factor = int(len(super_graph.MG_LO_cuts_corresponding_to_this_supergraph)/len(super_graph.cutkosky_cuts_generated))
                            # Actually this can happen, like in the supergraph 7x8 of e+ e- > a > t t~ g g featuring a purely gluonic bubble that must be associated
                            # with a factor 1/2.
#                            if reconstructed_symmetry_factor%fs_symm_factor!=0:
#                                log_func("%s%s: Incorrect LO cutkosky cut generation. %d not divisible by %d.%s"%(
#                                    color,err_text,
#                                    reconstructed_symmetry_factor,fs_symm_factor,
#                                    utils.bcolors.ENDC
#                                ) )

                        if is_LO:
                            super_graph.symmetry_factor = reconstructed_symmetry_factor/float(fs_symm_factor)
                        else:
                            if not _HAS_WARNED_ABOUT_SYMMETRY_FACTORS:
                                logger.critical("%sNo strategy yet for determining symmetry factors beyond LO contributions. Setting it to 1 for now, but this is incorrect.%s"%(
                                                                                                                                        utils.bcolors.RED,utils.bcolors.ENDC))
                                _HAS_WARNED_ABOUT_SYMMETRY_FACTORS = True

                            super_graph.symmetry_factor = 1

                        n_tot_cutkosky_cuts += fs_symm_factor * super_graph.symmetry_factor * max(len(super_graph.cutkosky_cuts_generated),0)
                        # For now the debug/check below always pass by constribution, but in the future we could think of computing super_graph.symmetry_factor independently
                        # from the number of MG_LO cuts.
                        if do_debug:
                            if fs_symm_factor*super_graph.symmetry_factor*len(super_graph.cutkosky_cuts_generated) != len(super_graph.MG_LO_cuts_corresponding_to_this_supergraph):
                                log_func(("%s%s: alphaLoop only found %d (fs_sym) x %d (graph_sym) x %d (n_cuts) LO cutkosky cuts in diagram '%s'%s, whereas the "+
                                        "following %d were expected from the following amplitude diagram interferences:\n%s%s")%(
                                    color,
                                    err_text,
                                    fs_symm_factor,
                                    super_graph.symmetry_factor,
                                    len(super_graph.cutkosky_cuts_generated),
                                    squared_topology_name,
                                    ' (FORM_ID: %d)'%FORM_id_to_supply if FORM_id_to_supply is not None else '',
                                    len(super_graph.MG_LO_cuts_corresponding_to_this_supergraph),
                                    ', '.join('%dx%d'%interf for interf in super_graph.MG_LO_cuts_corresponding_to_this_supergraph),
                                    utils.bcolors.ENDC
                                ))

                    FORM_id+=1
                    if (self.alphaLoop_options['n_rust_inputs_to_generate']>0) and (
                        (i_super_graph+1)==self.alphaLoop_options['n_rust_inputs_to_generate']):
                        break

            if do_debug:
                n_diag = len(matrix_element.get('diagrams'))
                target_n_tot_cutkosky_cuts = n_diag**2
                if n_tot_cutkosky_cuts != target_n_tot_cutkosky_cuts:
                    log_func("%s%s: alphaLoop only found %d LO cutkosky cuts whereas %d were expected from the %d amplitude diagrams generated!%s"%(
                        color,
                        err_text,
                        n_tot_cutkosky_cuts,target_n_tot_cutkosky_cuts,n_diag,
                        utils.bcolors.ENDC
                    ))

        if self.alphaLoop_options['FORM_processing_output_format']:

            FORM_super_graph_list = FORM_processing.FORMSuperGraphList(
                sum([sg_list for sg_list, ME in  self.all_super_graphs],[])
            )
            FORM_output_path = pjoin(self.dir_path, 'FORM')
            Path(FORM_output_path).mkdir(parents=True, exist_ok=True)
            Path(pjoin(self.dir_path, 'FORM', 'Rust_inputs')).mkdir(parents=True, exist_ok=True)

            FORM_workspace = pjoin(self.dir_path, 'FORM', 'workspace')
            Path(FORM_workspace).mkdir(parents=True, exist_ok=True)
            shutil.copy(pjoin(plugin_path, 'Templates', 'FORM_output_makefile'), 
                    pjoin(FORM_output_path, 'Makefile'))
            # Also output a copy of the input super-graph sent to the FORM processor.
            graph_dumps = []
            graph_dumps.append("graphs=[]")
            graph_dumps.append("graph_names=[]")
            graph_dumps.append("")
            for FORM_super_graph in FORM_super_graph_list:
                graph_dumps.append('graph_names.append("%s")'%FORM_super_graph[0].name)
                graph_dumps.append('graphs.append(\n%s\n)'%pformat(FORM_super_graph[0].to_dict()))
            open(pjoin(FORM_output_path,'all_MG_supergraphs.py'),'w').write(
                '\n'.join(graph_dumps)
            )
            computed_model = model_reader.ModelReader(model)
            computed_model.set_parameters_and_couplings(pjoin(self.dir_path,'Cards','param_card.dat')) 
            characteristic_process_definition = self.all_super_graphs[0][1].get('processes')[0]
            logger.info("Numerators processing with FORM...")
            FORM_processor = FORM_processing.FORMProcessor(
                FORM_processing.FORMSuperGraphList.from_dict(
                    pjoin(FORM_output_path,'all_MG_supergraphs.py'), 
                first=None,  merge_isomorphic_graphs=False),
#                FORM_super_graph_list,
                computed_model, characteristic_process_definition
            )

            drawings_output_path = pjoin(self.dir_path, 'Drawings')
            Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
            shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                        pjoin(drawings_output_path,'Makefile'))
            FORM_processor.draw(drawings_output_path)

            if self.alphaLoop_options['n_jets'] is not None:
                n_jets = self.alphaLoop_options['n_jets']
            else:
                n_jets = len([1 for leg in representative_process.get('legs') if 
                            leg.get('state')==True and leg.get('id') in self.alphaLoop_options['_jet_PDGs']])
                # Relax n_jets constraint based on the perturbative chosen
                if len(self.alphaLoop_options['perturbative_orders'])!=0:
                    n_jets -= sum(self.alphaLoop_options['perturbative_orders'].values())//2
            if self.alphaLoop_options['final_state_pdgs'] is not None:
                final_state_particle_ids = self.alphaLoop_options['final_state_pdgs']
            else:
                final_state_particle_ids = tuple([leg.get('id') for leg in representative_process.get('legs') if 
                            leg.get('state')==True and leg.get('id') not in self.alphaLoop_options['_jet_PDGs']])
            logger.info("Generating yaml input files with FORM processing...")
            FORM_processor.generate_squared_topology_files(
                pjoin(self.dir_path,'FORM','Rust_inputs'), n_jets, jet_ids=self.alphaLoop_options['_jet_PDGs'],
                final_state_particle_ids=final_state_particle_ids,
                filter_non_contributing_graphs=False,
                workspace=FORM_workspace,
                integrand_type=self.alphaLoop_options['FORM_integrand_type']
            )

            logger.info("Generating FORM numerators...")
            FORM_processor.generate_numerator_functions(FORM_output_path, 
                        output_format=self.alphaLoop_options['FORM_processing_output_format'],
                        workspace=FORM_workspace,
                        integrand_type=self.alphaLoop_options['FORM_integrand_type']
            )

        # And now finally generate the overall cross section yaml input file.
        base_name = os.path.basename(self.dir_path)
        overall_xsec_yaml_file_path = pjoin(rust_inputs_path,'%s.yaml'%base_name)
        overall_xsec_yaml = {
            'name' : base_name,
            'topologies' : [
                { 
                    'name' : super_graph.name,
                    'multiplicity' : super_graph.symmetry_factor
                }
                for all_super_graph, matrix_element in self.all_super_graphs for super_graph in all_super_graph
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
        shutil.copy(pjoin(plugin_path, 'Templates', 'Source_make_opts'), 
                    pjoin(self.dir_path, 'Source','make_opts'))

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

        # Fix makefiles compiler definition
        self.replace_make_opt_c_compiler(compiler['cpp'])
#       Uncomment if necessary but it needs to be fixed as for now it seems to corrupt make_opts
#        self.replace_make_opt_f_compiler(compiler)

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

        # Write the overall cross-section yaml file
        self.write_overall_LTD_squared_info()

        # Add the C_bindings
        filename = pjoin(self.dir_path, 'SubProcesses','C_bindings.f')
        C_binding_replace_dict = self.write_c_bindings(writers.FortranWriter(filename))

        # Add the IO_bindings
        filename = pjoin(self.dir_path, 'SubProcesses','IO_bindings.f')
        self.write_IO_bindings(writers.FortranWriter(filename),C_binding_replace_dict)

        # Compile the C_bindings
        logger.info("Compiling C_bindings for MG numerators...")
        try:
            misc.compile(arg=[],cwd=pjoin(self.dir_path, 'SubProcesses'), mode='fortran')
        except MadGraph5Error as e:
            logger.critical("%sMG numerator failed to compile:\n%s\nProceeding nonetheless.%s"%(
                utils.bcolors.RED, str(e), utils.bcolors.ENDC))
        
        # Compile FORM output if present
        FORM_output_dir = pjoin(self.dir_path,"FORM")
        if os.path.exists(FORM_output_dir):
            FORM_processing.FORMProcessor.compile(FORM_output_dir)
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
                proc_number_to_position_map.append('%d'%(position+1))
            except ValueError:
                proc_number_to_position_map.append('-1')
        replace_dict['proc_number_to_position_map'] = ','.join(proc_number_to_position_map)

        # Now gather all possible combination of number of external legs
        all_externals_combination = set([])
        for me, number in processes:
            all_externals_combination.add(len(me.get('processes')[0].get('legs')))
        
        max_n_diags = 0
        for super_graph_list, ME in self.all_super_graphs:
            for sg in super_graph_list:
                max_n_diags=max(max_n_diags,sg.call_signature['left_diagram_id'],sg.call_signature['right_diagram_id'])

        replace_dict['max_n_diags'] = max_n_diags

        max_n_indep_mom = max(all_externals_combination)

        mom_sign_flips = []
        for i_proc, (super_graph_list, ME) in enumerate(self.all_super_graphs):
            (nexternal, ninitial) = ME.get_nexternal_ninitial()
            for sg in super_graph_list:
                mom_signs= (['1',]*ninitial)+(['1',]*(len(sg.cuts)-1) if sg.MG_external_momenta_sign_flips is None else \
                        ['-1' if is_flipped else '1' for is_flipped in sg.MG_external_momenta_sign_flips])
                if len(mom_signs) > max_n_indep_mom:
                    logger.critical("The number of momenta sign in C_bindings.f is larger than the allowed maximum.")
                    logger.critical("The length is now adjusted so as to make sure C_bindings.f compiles, but this is indicative of a bug.")
                    max_n_indep_mom = len(mom_signs)
                elif max_n_indep_mom > len(mom_signs):
                    mom_signs += ['1',]*(max_n_indep_mom-len(mom_signs))

                mom_sign_flips.append("DATA (MOM_SIGNS(%d,%d,%d,I),I=1,%d)/%s/"%(
                    i_proc, 
                    sg.call_signature['left_diagram_id'],
                    sg.call_signature['right_diagram_id'],
                    len(mom_signs),
                    ','.join(mom_signs)
                ))

        replace_dict['max_n_indep_mom'] = max_n_indep_mom
        replace_dict['mom_sign_flips_def'] = '\n'.join(mom_sign_flips)

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
            n_incoming = len([l for l in me.get('processes')[0].get('legs') if l.get('state')==False])
            matrix_element_call_dispatch.append("P%d(I,%d)=DCMPLX(0.0d0,0.0d0)"%(n_external,n_external))
            matrix_element_call_dispatch.append("DO J=0,%d-1"%(n_external-1))
            matrix_element_call_dispatch.append("P%d(I,J+1)=MOM_SIGNS(PROC_ID,SDL,SDR,J+1)*DCMPLX(P(J*8+2*I+1),P(J*8+2*I+2))"%n_external)
            matrix_element_call_dispatch.append("IF ((J+1).le.%d) THEN"%n_incoming)
            matrix_element_call_dispatch.append("P%d(I,%d)=P%d(I,%d)+P%d(I,J+1)"%(n_external,n_external,n_external,n_external,n_external))
            matrix_element_call_dispatch.append("ELSE")
            matrix_element_call_dispatch.append("P%d(I,%d)=P%d(I,%d)-P%d(I,J+1)"%(n_external,n_external,n_external,n_external,n_external))
            matrix_element_call_dispatch.append("ENDIF")
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append(
                'CALL M%d_SMATRIXHEL(P%d, -1, SDL, SDR, ANS)'%(i_proc,n_external))
            overall_phase = complex(1,0)
            number_of_legs_for_each_SE={}
            for leg in me.get('processes')[0].get('legs'):
                # Consider final state legs only
                if leg.get('state')==False:
                    continue
                if abs(leg.get('id'))<LTD_squared.self_energy_global_id_offset:
                    overall_phase *= complex(0,1)
                SE_number = (abs(leg.get('id'))//LTD_squared.self_energy_global_id_offset)%10
                if SE_number!=0:
                    if SE_number in number_of_legs_for_each_SE:
                        number_of_legs_for_each_SE[SE_number] += 1
                    else:
                        number_of_legs_for_each_SE[SE_number] = 1
            # Then accumulate a phase for each additional legs in self-energies:
            if number_of_legs_for_each_SE:
                overall_phase *= (complex(0,-1)**(sum(v-3 for v in number_of_legs_for_each_SE.values())))
            matrix_element_call_dispatch.append("ANS=ANS*DCMPLX(%d.0d0,%d.0d0)"%(int(overall_phase.real),int(overall_phase.imag)))
            matrix_element_call_dispatch.append("GOTO 9999")
        replace_dict['matrix_element_call_dispatch'] = '\n'.join(matrix_element_call_dispatch)

        replace_dict['first_proc_prefix'] = 'M%d_'%(processes[0][1])

        template_path = pjoin(plugin_path,'Templates','C_bindings.inc')
        writer.writelines(open(template_path,'r').read()%replace_dict)

        return replace_dict

    def write_IO_bindings(self, writer, replace_dict):
        """ Write the fortran dispatcher to all processes callable through std I/O."""
        
        # Keep track of all ME exported so as to build the IO_binding dispatcher in finalize
        processes = self.all_processes_exported

        matrix_element_call_dispatch = []
        for i_proc, (me, proc_number) in enumerate(processes):
            matrix_element_call_dispatch.append('%d CONTINUE'%(i_proc+1))
            # Setup the kinematics now
            n_external = len(me.get('processes')[0].get('legs'))
            matrix_element_call_dispatch.append("DO J=1,%d"%n_external)
            matrix_element_call_dispatch.append("READ(*,*) (P(I),I=0,7)")
            matrix_element_call_dispatch.append("DO I=0,3")
            matrix_element_call_dispatch.append("P%d(I,J)=DCMPLX(P(2*I),P(2*I+1))"%n_external)
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append("ENDDO")
            matrix_element_call_dispatch.append(
                'CALL M%d_SMATRIXHEL(P%d, -1, SDL, SDR, ANS)'%(i_proc,n_external))
            matrix_element_call_dispatch.append("GOTO 9999")
        replace_dict['matrix_element_call_dispatch'] = '\n'.join(matrix_element_call_dispatch)

        template_path = pjoin(plugin_path,'Templates','IO_bindings.inc')
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

        if 'use_physical_gluon_helicity_sum' not in self.alphaLoop_options:
            raise alphaLoopExporterError("You probably forgot to specify '--mode=alphaloop' when running ./bin/mg5_aMC.")

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


        replace_dict['ndiags_in_this_proc'] = len(matrix_element.get('diagrams'))

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


class QGRAFExporter(object):
    
    qgraf_style = pjoin(plugin_src_path,'Templates','qgraf','orientedGraphPython.sty')
    qgraf_root_path = pjoin(plugin_src_path,'Templates','qgraf')
    qgraf_template_epem = pjoin(plugin_src_path,'Templates','qgraf','qgraf_epem.dat')
    qgraf_template_no_s = pjoin(plugin_src_path,'Templates','qgraf','qgraf_no_s.dat')

    qgraf_templates = {
        'epem':{'path': qgraf_template_epem, 'example': 'e+ e- > a > ...'},
        'no_s':{'paht': qgraf_template_no_s, 'example': 'a a > h'}
    }
    
    qgraf_field_replace = {}
    for field in ['d', 'u', 's', 'c', 'b', 't']:
        qgraf_field_replace[field] = field
        qgraf_field_replace[field+'~'] = field+'bar'
    qgraf_field_replace['gh'] = 'ghost'
    qgraf_field_replace['gh~'] = 'ghostbar'
    qgraf_field_replace['h'] = 'higgs'
    qgraf_field_replace['e-'] = 'eminus'
    qgraf_field_replace['e+'] = 'eplus'
    qgraf_field_replace['a'] = 'photon'
    qgraf_field_replace['g'] = 'gluon'
 

    default_opt = {'vetos': []}

    def __init__(self, process_definition, model, **opts):
        self.model = model
        self.proc_def = process_definition
        self.QGRAF_path = pjoin(plugin_path,os.path.pardir,'libraries','QGRAF','qgraf')
        # TODO: Ensure no multiparticles in the initial states
        self.initial_states = [leg.get('ids')[0] for leg in process_definition.get('legs') if not leg.get('state')]
        self.final_states = [leg.get('ids') for leg in process_definition.get('legs') if leg.get('state')]

    def generate(self):
        pass

    def output(self):
        pass

    
    def replace_make_opt_c_compiler(self, compiler, root_dir = ""):
        """Set CXX=compiler in Source/make_opts.
        The version is also checked, in order to set some extra flags
        if the compiler is clang (on MACOS)"""
       
        is_clang = misc.detect_if_cpp_compiler_is_clang(compiler)
        is_lc    = misc.detect_cpp_std_lib_dependence(compiler) == '-lc++'


        # list of the variable to set in the make_opts file
        for_update= {'DEFAULT_CPP_COMPILER':compiler,
                     'MACFLAG':'-mmacosx-version-min=10.7' if is_clang and is_lc else '',
                     'STDLIB': '-lc++' if is_lc else '-lstdc++',
                     'STDLIB_FLAG': '-stdlib=libc++' if is_lc and is_clang else ''
                     }

        # for MOJAVE remove the MACFLAG:
        if is_clang:
            import platform
            version, _, _ = platform.mac_ver()
            if not version:# not linux 
                version = 14 # set version to remove MACFLAG
            else:
                version = int(version.split('.')[1])
            if version >= 14:
                for_update['MACFLAG'] = '-mmacosx-version-min=10.8' if is_lc else ''

        if not root_dir:
            root_dir = self.dir_path
        make_opts = pjoin(root_dir, 'Source', 'make_opts')

        try:
            common_run_interface.CommonRunCmd.update_make_opts_full(
                            make_opts, for_update)
        except IOError:
            if root_dir == self.dir_path:
                logger.info('Fail to set compiler. Trying to continue anyway.')  
    
        return


class HardCodedAmpExporter():
        

    def __init__(self, process_definition, model,allow_mom_remap=False, **opts):
        
        self.MG5aMC_options = opts.pop('MG5aMC_options',{})
        self.dir_path = None
        self.amp_runcard = None
        self.alphaLoop_options = opts.pop('alphaLoop_options',{})
        self.model = model
        self.proc_def = process_definition
        self.color_dict = []
        # TODO: Ensure no multiparticles in the initial states
        self.initial_states = [leg.get('ids')[0] for leg in process_definition.get('legs') if not leg.get('state')]
        self.final_states = [leg.get('ids') for leg in process_definition.get('legs') if leg.get('state')]

        

    def output(self, output_path):
        
        self.amp_input = self.alphaLoop_options['AMPLITUDE_runcard_path']
               
        # Create amplitude with fake interference 
        outpath = os.path.abspath(output_path)
        amp = AMP_TOOLS.amplitude(runcard=self.amp_input, output=os.path.abspath(output_path),src_dir=plugin_path)
        color_dicts = amp.initialize_process(clean_up=True)
        external_data = amp.external_data

        # Overwrite default option of FORM processor to forbid UV-treatment
        FORM_processing.FORM_processing_options['generate_integrated_UV_CTs'] = False
        FORM_processing.FORM_processing_options['generate_renormalisation_graphs'] = False

        # loop over color
        for colstruc,dictpath in enumerate(color_dicts):
            self.dir_path = os.path.abspath(os.path.join(outpath,"color_struc_"+str(colstruc)))
            logger.info("Writing output of hardcoded amplitdue process to '%s'"%self.dir_path)    
            
            # Build output system
            self.build_output_directory()
            write_dir=pjoin(self.dir_path, 'Source', 'MODEL')
            model_builder = alphaLoopModelConverter(self.model, write_dir)
            model_builder.build([])
            shutil.copy(
                pjoin(self.dir_path,'Source','MODEL','param_card.dat'),
                pjoin(self.dir_path,'Cards','param_card.dat'),
            )
            shutil.copy(
                pjoin(self.dir_path,'Source','MODEL','ident_card.dat'),
                pjoin(self.dir_path,'Cards','ident_card.dat'),
            )

            computed_model = model_reader.ModelReader(self.model)
            computed_model.set_parameters_and_couplings(
                                                pjoin(self.dir_path,'Source','MODEL','param_card.dat'))   
            
            
            FORM_workspace = pjoin(self.dir_path, 'FORM', 'workspace')
            Path(FORM_workspace).mkdir(parents=True, exist_ok=True)

            # LOAD SGS
            super_graph_list = FORM_processing.FORMSuperGraphList.from_dict_amp(
                dictpath, verbose=True, model=self.model,external_data=external_data)       
            form_processor = FORM_processing.FORMProcessor(super_graph_list, computed_model, self.proc_def)
            shutil.copy(pjoin(plugin_path, 'Templates', 'FORM_output_makefile'), 
                    pjoin(self.dir_path,'FORM', 'Makefile'))

            if self.alphaLoop_options['n_rust_inputs_to_generate']<0:
                if self.alphaLoop_options['n_jets'] is None:
                    n_jets = len([1 for leg in self.proc_def.get('legs') if 
                                leg.get('state')==True and leg.get('ids')[0] in self.alphaLoop_options['_jet_PDGs']])
                else:
                    n_jets = self.alphaLoop_options['n_jets']
                if self.alphaLoop_options['final_state_pdgs'] is None:
                    final_state_particle_ids = tuple([leg.get('ids')[0] for leg in self.proc_def.get('legs') if 
                                leg.get('state')==True and leg.get('ids')[0] not in self.alphaLoop_options['_jet_PDGs']])
                else:
                    final_state_particle_ids = self.alphaLoop_options['final_state_pdgs']
                

                form_processor.generate_squared_topology_files(
                    pjoin(self.dir_path,'Rust_inputs'), n_jets, 
                    final_state_particle_ids=final_state_particle_ids,
                    jet_ids=self.alphaLoop_options['_jet_PDGs'],
                    # Remove non-contributing graphs from the list stored in the form_processor
                    filter_non_contributing_graphs=True,
                    workspace=FORM_workspace, 
                    integrand_type=self.alphaLoop_options['FORM_integrand_type']
                )

            # Draw
            drawings_output_path = pjoin(self.dir_path, 'Drawings')
            Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
            shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                        pjoin(drawings_output_path,'Makefile'))
            form_processor.draw(drawings_output_path)

            if self.alphaLoop_options['n_rust_inputs_to_generate']<0:

                form_processor.generate_numerator_functions(pjoin(self.dir_path,'FORM'), 
                    output_format=self.alphaLoop_options['FORM_processing_output_format'],
                    workspace=FORM_workspace,
                    integrand_type=self.alphaLoop_options['FORM_integrand_type'],
                    force_overall_factor=1
                )

            form_processor.compile(pjoin(self.dir_path,'FORM'))
    
    def build_output_directory(self):

        Path(self.dir_path).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Rust_inputs')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Cards')).mkdir(parents=True, exist_ok=True)
        # TODO add param_card in Source

        Path(pjoin(self.dir_path,'lib')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'FORM')).mkdir(parents=True, exist_ok=True)

        Path(pjoin(self.dir_path,'Source')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Source','MODEL')).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates', 'Source_make_opts'), 
                    pjoin(self.dir_path, 'Source','make_opts'))

        # Fix makefiles compiler definition
        self.replace_make_opt_c_compiler(self.MG5aMC_options['cpp_compiler'])

    def replace_make_opt_c_compiler(self, compiler, root_dir = ""):
        """Set CXX=compiler in Source/make_opts.
        The version is also checked, in order to set some extra flags
        if the compiler is clang (on MACOS)"""
       
        is_clang = misc.detect_if_cpp_compiler_is_clang(compiler)
        is_lc    = misc.detect_cpp_std_lib_dependence(compiler) == '-lc++'


        # list of the variable to set in the make_opts file
        for_update= {'DEFAULT_CPP_COMPILER':compiler,
                     'MACFLAG':'-mmacosx-version-min=10.7' if is_clang and is_lc else '',
                     'STDLIB': '-lc++' if is_lc else '-lstdc++',
                     'STDLIB_FLAG': '-stdlib=libc++' if is_lc and is_clang else ''
                     }

        # for MOJAVE remove the MACFLAG:
        if is_clang:
            import platform
            version, _, _ = platform.mac_ver()
            if not version:# not linux 
                version = 14 # set version to remove MACFLAG
            else:
                version = int(version.split('.')[1])
            if version >= 14:
                for_update['MACFLAG'] = '-mmacosx-version-min=10.8' if is_lc else ''

        if not root_dir:
            root_dir = self.dir_path
        make_opts = pjoin(root_dir, 'Source', 'make_opts')

        try:
            common_run_interface.CommonRunCmd.update_make_opts_full(
                            make_opts, for_update)
        except IOError:
            if root_dir == self.dir_path:
                logger.info('Fail to set compiler. Trying to continue anyway.')  
    
        return







class HardCodedQGRAFExporter(QGRAFExporter):

    def __init__(self, process_definition, model, **opts):
        self.alphaLoop_options = opts.pop('alphaLoop_options',{})
        self.MG5aMC_options = opts.pop('MG5aMC_options',{})
        self.dir_path = None
        self.qgraf_raw_output = None
        self.qgraf_output = None
        self.qgraf_model = pjoin(self.qgraf_root_path, '%s.dat'%self.alphaLoop_options['qgraf_model'])

        if isinstance(process_definition, base_objects.ProcessDefinition):
            all_processes = list(proc for proc in process_definition)
            representative_proc = all_processes[0]
        else:
            representative_proc = process_definition
        
        super(HardCodedQGRAFExporter,self).__init__(process_definition, model, **opts)
        
        #print(self.final_states)
        n_cuts = len(self.final_states)
        virtual_loops = len(representative_proc['perturbation_couplings'])

        overall_phase = -complex(-1.0,0.0)**(n_cuts + virtual_loops)
        if overall_phase.imag != 0:
            raise alphaLoopExporterError("No support for overall complex phase yet (Ben: how do we put a complex number in FORM? ^^)")
        else:
            self.overall_phase = '%d'%int(overall_phase.real)



    def output(self, output_path):
        self.dir_path = os.path.abspath(output_path)
        self.qgraf_raw_output = pjoin(self.dir_path,'qgraf','output.py')
        self.qgraf_output = pjoin(self.dir_path,'qgraf','all_QG_supergraphs.py')

        # Checkpoints 
        self.checkpoints_path = pjoin(self.dir_path,'checkpoints')
        self.checkpoint = [[pjoin(self.checkpoints_path, name) for name in names]
                             for names in [['FORM_iso_sq_list_1'],['FORM_iso_sq_list_2', 'FORM_processor_2']]] 

        # Preparing Drawing Folder
        drawings_output_path = pjoin(self.dir_path, 'Drawings')
        Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                    pjoin(drawings_output_path,'Makefile'))
        
        logger.info("Writing output of hardcoded QGRAF generation to '%s'"%self.dir_path)

        # Extract process information
        if isinstance(self.proc_def, base_objects.ProcessDefinition):
             representative_proc = next(proc for proc in self.proc_def)
        else:
            representative_proc = self.proc_def

        # Generate all supergraphs using QGRAF
        if self.alphaLoop_options['checkpoint_lvl'] < 1 or not os.path.isfile(self.qgraf_output):
            self.build_output_directory()
            getattr(self,"build_qgraf_%s"%self.alphaLoop_options['qgraf_template_model'])(representative_proc)
            self.standalone_qgraf_file()

        # Process supergraph numerators with FORM and output result in the process output.
        write_dir=pjoin(self.dir_path, 'Source', 'MODEL')
        model_builder = alphaLoopModelConverter(self.model, write_dir)
        model_builder.build([])
        shutil.copy(
            pjoin(self.dir_path,'Source','MODEL','param_card.dat'),
            pjoin(self.dir_path,'Cards','param_card.dat'),
        )
        shutil.copy(
            pjoin(self.dir_path,'Source','MODEL','ident_card.dat'),
            pjoin(self.dir_path,'Cards','ident_card.dat'),
        )

        # Example of how to access information directly from a param_card.dat
        #        param_card = check_param_card.ParamCard(pjoin(self.dir_path,'Source','MODEL','param_card.dat'))
        #        misc.sprint(param_card.get_value('mass',25))
        #        misc.sprint(param_card.get_value('width',25))
        #        misc.sprint(param_card.get_value('yukawa',6))
        #        misc.sprint(param_card.get_value('sminputs',1)) # aEWM1
        #        misc.sprint(param_card.get_value('sminputs',2)) # Gf
        #        misc.sprint(param_card.get_value('sminputs',3)) # aS

        computed_model = model_reader.ModelReader(self.model)
        computed_model.set_parameters_and_couplings(
                                            pjoin(self.dir_path,'Source','MODEL','param_card.dat'))        

        # Assign the name of the model to be qgraf name, as it is useful for deciding between SM and HEFT renormalisation
        computed_model.set('name',self.alphaLoop_options['qgraf_model'])

        final_state_particle_ids = [ abs(pdg) for pdg in self.alphaLoop_options['_jet_PDGs']]
        if self.alphaLoop_options['qgraf_cut_filter']:
            cuts = self.get_cuts(representative_proc, final_state_particle_ids=final_state_particle_ids+[None])
        else:
            cuts = None

        FORM_workspace = pjoin(self.dir_path, 'FORM', 'workspace')
        
        if self.alphaLoop_options['checkpoint_lvl'] < 1 or not all(os.path.isfile(cp) for cp in self.checkpoint[0]):
            Path(FORM_workspace).mkdir(parents=True, exist_ok=True)
            super_graph_list = FORM_processing.FORMSuperGraphList.from_dict(
                self.qgraf_output, merge_isomorphic_graphs=True, 
                model=self.model, workspace=FORM_workspace, cuts=cuts)
            pickle.dump(super_graph_list, open(self.checkpoint[0][0],'wb'))
            logger.info("\033[1:32mFirst Checkpoint!\033[0m")
        else:
            logger.info("\033[1:32mRecovering from first checkpoint!\033[0m")
            super_graph_list = pickle.load(open(self.checkpoint[0][0],'rb'))

        # Add phase
        for graphs in super_graph_list:
            for g in graphs:
                g.overall_factor = "({})*({})"\
                    .format(self.overall_phase, g.overall_factor)

        if self.alphaLoop_options['checkpoint_lvl'] < 2 or not all(os.path.isfile(cp) for cp in self.checkpoint[1]):
            # Create FORM_processing instance and generate the squared topologies
            form_processor = FORM_processing.FORMProcessor(super_graph_list, computed_model, self.proc_def)
            shutil.copy(pjoin(plugin_path, 'Templates', 'FORM_output_makefile'), 
                    pjoin(self.dir_path,'FORM', 'Makefile'))

            if self.alphaLoop_options['n_rust_inputs_to_generate']<0:
                if self.alphaLoop_options['n_jets'] is None:
                    n_jets = len([1 for leg in representative_proc.get('legs') if 
                                leg.get('state')==True and leg.get('id') in self.alphaLoop_options['_jet_PDGs']])
                else:
                    n_jets = self.alphaLoop_options['n_jets']
                if self.alphaLoop_options['final_state_pdgs'] is None:
                    final_state_particle_ids = tuple([leg.get('id') for leg in representative_proc.get('legs') if 
                                leg.get('state')==True and leg.get('id') not in self.alphaLoop_options['_jet_PDGs']])
                else:
                    final_state_particle_ids = self.alphaLoop_options['final_state_pdgs']

                #print(n_jets)
                #print(final_state_particle_ids)
                #print(self.alphaLoop_options['_jet_PDGs'])
                form_processor.generate_squared_topology_files(
                    pjoin(self.dir_path,'Rust_inputs'), n_jets, 
                    final_state_particle_ids=final_state_particle_ids,
                    jet_ids=self.alphaLoop_options['_jet_PDGs'],
                    # Remove non-contributing graphs from the list stored in the form_processor
                    filter_non_contributing_graphs=True,
                    workspace=FORM_workspace, 
                    integrand_type=self.alphaLoop_options['FORM_integrand_type']
                )

            # Draw all the graphs contributing to the process
            form_processor.draw(drawings_output_path)
            
            # Checkpoint 2
            pickle.dump(super_graph_list, open(self.checkpoint[1][0],'wb'))
            pickle.dump(form_processor, open(self.checkpoint[1][1],'wb'))
            logger.info("\033[1:32mSecond Checkpoint!\033[0m")
        else:
            logger.info("\033[1:32mRecovering from second checkpoint!\033[0m")
            super_graph_list = pickle.load(open(self.checkpoint[1][0],'rb'))
            form_processor = pickle.load(open(self.checkpoint[1][1],'rb'))


        if self.alphaLoop_options['n_rust_inputs_to_generate']<0:

            form_processor.generate_numerator_functions(pjoin(self.dir_path,'FORM'), 
                output_format=self.alphaLoop_options['FORM_processing_output_format'],
                workspace=FORM_workspace,
                integrand_type=self.alphaLoop_options['FORM_integrand_type']
            )

        form_processor.compile(pjoin(self.dir_path,'FORM'))

    def get_cuts(self, representative_process, final_state_particle_ids='any'):
        cuts=[]
        additional_loops = len(representative_process['perturbation_couplings'])
        for pdgs in self.final_states:
            cuts += [([abs(pdg) for pdg in pdgs],1)]
        #cuts += [('any', additional_loops )]
        cuts += [(final_state_particle_ids, additional_loops )]
        return cuts

    def build_qgraf_epem(self, representative_process):
        # Check if e+ e- > a > ...
        check_initial_states = [leg.get('id') for leg in representative_process.get('legs') if leg.get('state')==False]==[-11,11]
        check_s_channel = representative_process.get('required_s_channels') == [[22]]
        if not (check_initial_states and check_s_channel):
            raise alphaLoopExporterError("The command 'qgraf_generate' with model 'epem' does not support this process.\nAvailable models (: example)\n{}"\
                .format("\n".join(["\t%s: %s"%(k,v['example']) for k,v in self.qgraf_templates.items()]))
            )

        if self.alphaLoop_options['loop_induced']:
            virtual_loops = len(representative_process['perturbation_couplings']) + 2
        else:
            virtual_loops = len(representative_process['perturbation_couplings'])
        n_final_states = len([leg.get('id') for leg in representative_process.get('legs') if leg.get('state')==True])

        pdg_model_map = self.model['particles'].generate_dict()

        dict_replace = {}
        dict_replace['n_loops'] = n_final_states - 1 + virtual_loops 

        # Veto particles that are forbidden
        dict_replace['vetos'] = ''
        for particle in representative_process['forbidden_particles']:
            dict_replace['vetos'] += 'true=iprop[%s,0,0];\n'%self.qgraf_field_replace[pdg_model_map[particle]['name']]

        # Enforce coupling order
        coupling_order_id = 0
        for sq_coupling, order in representative_process['squared_orders'].items():
            if sq_coupling == 'QED':
                coupling_order_id += 2*order
            elif sq_coupling == 'QCD':
                coupling_order_id += 200*order
            else:
                raise alphaLoopExporterError("Unsupported Coupling for QGRAF: %s"%sq_coupling)
        dict_replace['coupling_order'] = 'true=vsum[QCD_QED,{0},{0}];\n'.format(coupling_order_id)

        # Ensure final state
        dict_replace['final_states'] = ''
        required_final_states = [fields[0] for fields in self.final_states if len(fields) == 1]
        for field in set(required_final_states):
            if field >= 0:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['name']]
            else:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['antiname']]
            dict_replace['final_states'] += 'true=iprop[{},{},10];\n'.format(field_name, required_final_states.count(field))

        shutil.copy(self.qgraf_model, pjoin(self.dir_path, 'qgraf','model.dat'))
        shutil.copy(self.qgraf_style, pjoin(self.dir_path, 'qgraf','orientedGraphPython.sty'))
        
        qgraf_output = pjoin(self.dir_path, 'qgraf/qgraf.dat')
        qgraf_folder = pjoin(self.dir_path, 'qgraf')

        with open(self.qgraf_template_epem, 'r') as stream:
            with open(qgraf_output, 'w') as f:
                f.write(stream.read() % dict_replace)

            subprocess.run(['rm', 'output.py'], cwd=qgraf_folder)
            r = subprocess.run([self.QGRAF_path, ],
                               cwd=qgraf_folder,
                               capture_output=True)
            if r.returncode != 0 or not os.path.exists(pjoin(qgraf_folder, 'output.py')):
                raise print("QGRAF generation failed with error:\n%s" %
                            (r.stdout.decode('UTF-8')))

    def build_qgraf_no_s(self, representative_process, vetos=None):
        # Check if a1 a2 a3 ... > a4 a5 ...
        check_s_channel = representative_process.get('required_s_channels') == []
        if not check_s_channel:
            raise alphaLoopExporterError("The command 'qgraf_generate' with model 'no_s' does not support this process.\nAvailable models (: example)\n{}"\
                .format("\n".join(["\t%s: %s"%(k,v['example']) for k,v in self.qgraf_templates.items()]))
            )
        if self.alphaLoop_options['loop_induced']:
            virtual_loops = len(representative_process['perturbation_couplings']) + 2 
        else:
            virtual_loops = len(representative_process['perturbation_couplings'])
        n_final_states = len([leg.get('id') for leg in representative_process.get('legs') if leg.get('state')==True])

        pdg_model_map = self.model['particles'].generate_dict()
        
        dict_replace = {}
        dict_replace['n_loops'] = n_final_states - 1 + virtual_loops
        
        # Veto particles that are forbidden
        dict_replace['vetos'] = ''
        for particle in representative_process['forbidden_particles']:
            dict_replace['vetos'] += 'true=iprop[%s,0,0];\n'%self.qgraf_field_replace[pdg_model_map[particle]['name']]

        # Enforce coupling order
        coupling_order_id = 0
        for sq_coupling, order in representative_process['squared_orders'].items():
            if sq_coupling == 'QED':
                coupling_order_id += 2*order
            elif sq_coupling == 'QCD':
                coupling_order_id += 200*order
            else:
                raise alphaLoopExporterError("Unsupported Coupling for QGRAF: %s"%sq_coupling)
        dict_replace['coupling_order'] = 'true=vsum[QCD_QED,{0},{0}];\n'.format(coupling_order_id)

        # Ensure final state
        dict_replace['final_states'] = ''
        required_final_states = [fields[0] for fields in self.final_states if len(fields) == 1]
        for field in set(required_final_states):
            if field >= 0:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['name']]
            else:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['antiname']]
            dict_replace['final_states'] += 'true=iprop[{},{},10];\n'.format(field_name, required_final_states.count(field))

        shutil.copy(self.qgraf_model, pjoin(self.dir_path, 'qgraf','model.dat'))
        shutil.copy(self.qgraf_style, pjoin(self.dir_path, 'qgraf','orientedGraphPython.sty'))
        
        qgraf_output = pjoin(self.dir_path, 'qgraf/qgraf.dat')
        qgraf_folder = pjoin(self.dir_path, 'qgraf')

        with open(self.qgraf_template_epem, 'r') as stream:
            with open(qgraf_output, 'w') as f:
                f.write(stream.read() % dict_replace)

            subprocess.run(['rm', 'output.py'], cwd=qgraf_folder)
            r = subprocess.run([self.QGRAF_path, ],
                               cwd=qgraf_folder,
                               capture_output=True)
            if r.returncode != 0 or not os.path.exists(pjoin(qgraf_folder, 'output.py')):
                raise print("QGRAF generation failed with error:\n%s" %
                            (r.stdout.decode('UTF-8')))


        # Veto particles that are forbidden
        dict_replace['vetos'] = ''
        for particle in representative_process['forbidden_particles']:
            dict_replace['vetos'] += 'true=iprop[%s,0,0];\n'%self.qgraf_field_replace[pdg_model_map[particle]['name']]

        # Enforce coupling order
        coupling_order_id = 0
        for sq_coupling, order in representative_process['squared_orders'].items():
            if sq_coupling == 'QED':
                coupling_order_id += 2*order
            elif sq_coupling == 'QCD':
                coupling_order_id += 200*order
            else:
                raise alphaLoopExporterError("Unsupported Coupling for QGRAF: %s"%sq_coupling)
        dict_replace['coupling_order'] = 'true=vsum[QCD_QED,{0},{0}];\n'.format(coupling_order_id)

        # Set externals and enforce final states
        dict_replace['incoming_states'] = 'in = '
        dict_replace['outgoing_states'] = 'out = '
        incoming_fields = []
        pn = 1
        for field in self.initial_states:
            if field >= 0:
                incoming_fields += [(self.qgraf_field_replace[self.model.get_particle(abs(field))['name']],'p%d'%pn)]
            else:
                incoming_fields += [(self.qgraf_field_replace[self.model.get_particle(abs(field))['antiname']],'p%d'%pn)]
            pn += 1
        dict_replace['incoming_states'] += ", ".join(["{}[{}]".format(*field) for field in incoming_fields])+";"
        dict_replace['outgoing_states'] += ", ".join(["{}[{}]".format(*field) for field in incoming_fields])+";"

        dict_replace['final_states'] = ''
        required_final_states = [fields[0] for fields in self.final_states if len(fields) == 1]
        for field in set(required_final_states):
            if field >= 0:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['name']]
            else:
                field_name = self.qgraf_field_replace[self.model.get_particle(abs(field))['antiname']]
            dict_replace['final_states'] += 'true=iprop[{},{},10];\n'.format(field_name, required_final_states.count(field))

        shutil.copy(self.qgraf_model, pjoin(self.dir_path, 'qgraf','model.dat'))
        shutil.copy(self.qgraf_style, pjoin(self.dir_path, 'qgraf','orientedGraphPython.sty'))
        
        qgraf_output = pjoin(self.dir_path, 'qgraf/qgraf.dat')
        qgraf_folder = pjoin(self.dir_path, 'qgraf')

        with open(self.qgraf_template_no_s, 'r') as stream:
            with open(qgraf_output, 'w') as f:
                f.write(stream.read() % dict_replace)

            subprocess.run(['rm', 'output.py'], cwd=qgraf_folder)
            r = subprocess.run([self.QGRAF_path, ],
                               cwd=qgraf_folder,
                               capture_output=True)
            if r.returncode != 0 or not os.path.exists(pjoin(qgraf_folder, 'output.py')):
                raise print("QGRAF generation failed with error:\n%s" %
                            (r.stdout.decode('UTF-8')))

    def build_output_directory(self):

        Path(self.dir_path).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'checkpoints')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Rust_inputs')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Cards')).mkdir(parents=True, exist_ok=True)
        # TODO add param_card in Source

        Path(pjoin(self.dir_path,'lib')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'FORM')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'qgraf')).mkdir(parents=True, exist_ok=True)

        Path(pjoin(self.dir_path,'Source')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Source','MODEL')).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates', 'Source_make_opts'), 
                    pjoin(self.dir_path, 'Source','make_opts'))

        # Fix makefiles compiler definition
        self.replace_make_opt_c_compiler(self.MG5aMC_options['cpp_compiler'])


    def standalone_qgraf_file(self):
        logger.info('Create standalone QGRAF from: {}'.format(self.qgraf_output))

        with open(self.qgraf_output, 'w') as f:
            f.write("graphs=[]\n")
            f.write("graph_names=[]\n")
            with open(self.qgraf_raw_output, 'r') as stream:
                n_graph = 0
                for line in stream.read().splitlines():
                    if 'graphs=[]' in line or 'graph_names=[]' in line:
                        continue
                    if 'graphs.append' in line:
                        f.write('graph_names+=[\"SG_QG{}\"]\n'.format(n_graph))
                        n_graph += 1
                    f.write(line+'\n')

class LUScalarTopologyExporter(QGRAFExporter):

    def __init__(self, cli, output_path, topology, externals, name, lmb, model, benchmark_result=0.0, **opts):
        
        self.alphaLoop_options = opts.pop('alphaLoop_options',{})
        self.MG5aMC_options = opts.pop('MG5aMC_options',{})
        self.cli = cli
        self.dir_path = os.path.abspath(output_path)
        self.topology = topology
        self.externals = externals
        self.name = name
        self.lmb = lmb
        self.benchmark_result = benchmark_result

        self.model = copy.deepcopy(model)
        self.model['particles'].append(base_objects.Particle({
            'name': 'psi',
            'antiname': 'psi',
            'spin' : 1,
            'color': 1,
            # We can think of generalising the pipeline for massive scalars too.
            'mass': 'ZERO',
            'width': 'ZERO', 
            'pdg_code': 1337,
            'line': 'dashed',
            'is_part':True,
            'self_antipart': True
        }))
        # We must add the scalar particle to it!
        self.model.reset_dictionaries()

    def build_output_directory(self):

        Path(self.dir_path).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Rust_inputs')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Cards')).mkdir(parents=True, exist_ok=True)
        # TODO add param_card in Source

        Path(pjoin(self.dir_path,'lib')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'FORM')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'FORM','workspace')).mkdir(parents=True, exist_ok=True)

        Path(pjoin(self.dir_path,'Source')).mkdir(parents=True, exist_ok=True)
        Path(pjoin(self.dir_path,'Source','MODEL')).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates', 'Source_make_opts'), 
                    pjoin(self.dir_path, 'Source','make_opts'))

        # Fix makefiles compiler definition
        self.replace_make_opt_c_compiler(self.MG5aMC_options['cpp_compiler'])

    def get_sg_list(self,computed_model):

        return FORM_processing.FORMSuperGraphList.from_squared_topology(
            self.topology, self.name, self.externals, computed_model, loop_momenta_names=self.lmb,
            benchmark_result=self.benchmark_result
        )
    
    def get_overall_factor(self):
        return '-1' # there is an overall factor of -1 wrt Forcer, presumable due to the Wick rotation?

    def output(self):

        # Generate all supergraphs using QGRAF
        self.build_output_directory()

        # Process supergraph numerators with FORM and output result in the process output.
        write_dir=pjoin(self.dir_path, 'Source', 'MODEL')
        model_builder = alphaLoopModelConverter(self.model, write_dir)
        model_builder.build([])
        shutil.copy(
            pjoin(self.dir_path,'Source','MODEL','param_card.dat'),
            pjoin(self.dir_path,'Cards','param_card.dat'),
        )
        shutil.copy(
            pjoin(self.dir_path,'Source','MODEL','ident_card.dat'),
            pjoin(self.dir_path,'Cards','ident_card.dat'),
        )

        # Example of how to access information directly from a param_card.dat
        #        param_card = check_param_card.ParamCard(pjoin(self.dir_path,'Source','MODEL','param_card.dat'))
        #        misc.sprint(param_card.get_value('mass',25))
        #        misc.sprint(param_card.get_value('width',25))
        #        misc.sprint(param_card.get_value('yukawa',6))
        #        misc.sprint(param_card.get_value('sminputs',1)) # aEWM1
        #        misc.sprint(param_card.get_value('sminputs',2)) # Gf
        #        misc.sprint(param_card.get_value('sminputs',3)) # aS

        computed_model = model_reader.ModelReader(self.model)
        computed_model.set_parameters_and_couplings(
                                            pjoin(self.dir_path,'Source','MODEL','param_card.dat'))

        # Dummy process definition is good enough in this case.
        process_definition=self.cli.extract_process('e+ e- > a > d d~', proc_number=0)

        form_sg_list = self.get_sg_list(computed_model)

        form_processor = FORM_processing.FORMProcessor(form_sg_list, computed_model, process_definition)

        FORM_output_path = pjoin(self.dir_path,'FORM')
        FORM_workspace = pjoin(FORM_output_path,'workspace')
        form_processor.generate_squared_topology_files(
            pjoin(self.dir_path,'Rust_inputs'), 0, final_state_particle_ids=(), jet_ids=None, 
            filter_non_contributing_graphs=False, workspace=FORM_workspace,
            integrand_type=self.alphaLoop_options['FORM_integrand_type']
        )

        form_processor.generate_numerator_functions(
            FORM_output_path, output_format=self.alphaLoop_options['FORM_processing_output_format'],
            workspace=FORM_workspace, header="", integrand_type=self.alphaLoop_options['FORM_integrand_type'],
            force_overall_factor=self.get_overall_factor()
        )

        # Draw
        drawings_output_path = pjoin(self.dir_path, 'Drawings')
        Path(drawings_output_path).mkdir(parents=True, exist_ok=True)
        shutil.copy(pjoin(plugin_path, 'Templates','Drawings_makefile'),
                    pjoin(drawings_output_path,'Makefile'))
        form_processor.draw(drawings_output_path)

        # Compile
        shutil.copy(pjoin(plugin_path, 'Templates', 'FORM_output_makefile'), 
                    pjoin(FORM_output_path, 'Makefile'))
        form_processor.compile(FORM_output_path)

class ScalarIntegralTopologyExporter(LUScalarTopologyExporter):

    def __init__(self, cli, output_path, topology, externals, default_kinematics, masses, name, lmb, model, benchmark_result=0.0,  numerator=None, **opts):
        
        self.default_kinematics = default_kinematics
        # TODO add support for masses and numerators
        if masses is not None:
            raise alphaLoopExporterError("No support yet for masses in ScalarIntegralTopologyExporter.")
        self.masses = masses

        self.numerator = numerator

        super(ScalarIntegralTopologyExporter, self).__init__(cli, output_path, topology, externals, name, lmb, model, benchmark_result=benchmark_result, **opts)

        # We now need to create a squared topology from the provided scalar integral topology
        self.squared_topology_info = self.create_squared_topology()

    def create_squared_topology(self):
        """ Creates and LU squared effective topology from a scalar integral topology."""

        # This function must return the following info
        squared_topology_info = {
            'squared_topology': None,
            'lmb': None,
            'externals': None,
            'effective_vertex_id': None,
            'default_kinematics': None
        }

        # Create a map of the edge name to nodes connected
        edges_map = {e: (u, v) for e,u,v in self.topology}
        external_edges = {}
        internal_edges = []
        for e, nodes in edges_map.items():
            lone_nodes = [ n for n in nodes if all(n not in other_nodes for e_i, other_nodes in edges_map.items() if e_i!=e) ]
            if len(lone_nodes)==2:
                raise alphaLoopExporterError("Edge '%s' of the topology is completely disconnected."%e)
            elif len(lone_nodes)==1:
                external_edges[e] = lone_nodes[0]
            else:           
                internal_edges.append(e)
    
        # First scan the externals of the scalar integral topology and make sure they are properly oriented.
        incoming_nodes = []
        for incoming_edge in self.externals[0]:
            if incoming_edge not in external_edges:
                raise alphaLoopExporterError("Specified external incoming edge '%s' could not be found."%incoming_edge)
            if external_edges[incoming_edge] != edges_map[incoming_edge][0]:
                raise alphaLoopExporterError("Incorrect direction for incoming external edge '%s'."%incoming_edge)
            incoming_nodes.append(external_edges[incoming_edge])

        outgoing_nodes = []
        for outgoing_edge in self.externals[1]:
            if outgoing_edge not in external_edges:
                raise alphaLoopExporterError("Specified external outgoing edge '%s' could not be found."%outgoing_edge)
            if external_edges[outgoing_edge] != edges_map[outgoing_edge][1]:
                raise alphaLoopExporterError("Incorrect direction for outgoing external edge '%s'."%outgoing_edge)
            outgoing_nodes.append(external_edges[outgoing_edge])

        # Now sew external outgoing edges with an effective vertex.
        # First copy edges that will remain in the squared topology
        squared_topology_info['squared_topology'] = [(e, edges_map[e][0], edges_map[e][1]) for e in internal_edges+list(self.externals[0])]
        
        # Create the effectie node ID
        max_node_id = 0
        for u,v in edges_map.values():
            max_node_id = max(max_node_id,u,v)
        max_node_id += 1
        effective_node_id = max_node_id
        squared_topology_info['effective_vertex_id']=effective_node_id

        # Connect the outgoing externals to the effective node_id
        for outgoing_edge in self.externals[1]:
            # Overwrite the external node to be that of the effective node id instead
            squared_topology_info['squared_topology'].append( (outgoing_edge, edges_map[outgoing_edge][0], effective_node_id) )

        # The externals are simply the specified incoming edges
        squared_topology_info['externals'] = list(self.externals[0])
        
        # Create copies of the incoming edges
        external_incoming_number = 1
        while 'q%i'%external_incoming_number in edges_map:
            external_incoming_number +=1
        for _in_node in incoming_nodes:
            max_node_id += 1
            squared_topology_info['squared_topology'].append( ('q%i'%external_incoming_number, effective_node_id, max_node_id) )
            external_incoming_number += 1

        # Combine the user's specified LMB with the LMB from the externals
        if self.lmb is not None:
            squared_topology_info['lmb'] = list(self.lmb) + list(self.externals[1])[:-1]
        else:
            squared_topology_info['lmb'] = None

        # Fill in the specified default kinematics as a two-tuple of list of momenta
        squared_topology_info['default_kinematics'] = [
            [ self.default_kinematics[incoming_edge] for incoming_edge in self.externals[0] ],
            [ self.default_kinematics[outgoing_edge] for outgoing_edge in list(self.externals[1])[:-1] ],
        ]

        # We must divide the (-i)^N of the N-point effective mock-up scalar vertex on the right of the cutkosky cut.
        squared_topology_info['overall_factor'] = -1./(complex(0,-1.)**(len(self.externals[0])*2-2)) # extra minus one of unknown origin

        return squared_topology_info

    def get_sg_list(self,computed_model):

        return FORM_processing.FORMSuperGraphList.from_squared_topology(
            self.squared_topology_info['squared_topology'], self.name, 
            self.squared_topology_info['externals'], computed_model, 
            loop_momenta_names=self.squared_topology_info['lmb'],
            benchmark_result=self.benchmark_result,
            default_kinematics=self.squared_topology_info['default_kinematics'],
            effective_vertex_id=self.squared_topology_info['effective_vertex_id'],
        )

    def get_overall_factor(self):
        # Dress the overall factor with the numerator
        return '((%d)+(%d*i_))*(%s)'%(
            self.squared_topology_info['overall_factor'].real,self.squared_topology_info['overall_factor'].imag,
            '1' if self.numerator is None else self.numerator
        )
