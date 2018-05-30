#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruij               #
#                                                   #
#####################################################

import os
import logging
import sys
import shutil
import re
import random

plugin_path = os.path.dirname(os.path.realpath( __file__ ))
import madgraph

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite, MPI_RANK, MPI_SIZE, MPI_ACTIVE
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
import madgraph.various.cluster as cluster

from madgraph.iolibs.files import cp, ln, mv

import madgraph.integrator.phase_space_generators as phase_space_generators
from madgraph.integrator.vegas3_integrator import Vegas3Integrator
# It will be nice to experiment with Thomas Hahn's integrators too later
# import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator

import nloop_integrands


logger = logging.getLogger('pyNLoop.Interface')

pjoin = os.path.join
template_dir = pjoin(plugin_path, 'Templates')

class pyNLoopInterfaceError(MadGraph5Error):
    """ Error for the pyNLoop plugin """

class pyNLoopInvalidCmd(InvalidCmd):
    """ Invalid command issued to the pyNLoop interface. """

class pyNLoopInterface(madgraph_interface.MadGraphCmd, cmd.CmdShell):
    """ Interface for steering the generation/output of pyNLoop.
    We make it inherit from CmdShell so that launch_ext_prog does not attempt to start in WebMode."""
    
    _hardcoded_topologies = { 
        'dummy_function' : 
            {   'class'     :   nloop_integrands.DummyNLoopIntegrand,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [100., 200., 300., 400.]
            },
        'scalar_box_full_offshell' : 
            {   'class'     :   nloop_integrands.HardCodedOffshellScalarBox,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [100., 200., 300., 400.]
            }
    }
    
    def __init__(self, *args, **opts):
        """ Define attributes of this class."""
        
        self.pyNLoop_options = {
            # At first, we will limit ourselves to using a built-in distribution of pySecDec.
            'pySecDec_path'     :   'built_in',
            'integrator'        :   'Vegas3',
            'parallelization'   :   'multicore'
        }
        super(pyNLoopInterface, self).__init__(*args, **opts)
        
    def do_set_pyNLoop_option(self, line):
        """ Logic for setting pyNLoop options."""

        args = self.split_arg(line)
        key, value = args[:2]
 
        if key == 'integrator':
            if value not in ['Vegas3']:
                raise pyNLoopInvalidCmd("pyNLoop only supports Vegas3 integrator for now.")
        
        elif key == 'parallelization':
            if value not in ['cluster', 'multicore']:
                raise pyNLoopInvalidCmd("pyNLoop only supports parallelization "+
                                                        "modes 'cluser' and 'multicore'.")

        else:
            raise pyNLoopInvalidCmd("Unrecognized pyNLoop option: %s"%key)
        
    def get_cluster(self, force_nb_cpu_cores=None):
        """setup the number of cores for multicore, and the cluster-type for cluster runs"""

        options_to_pass_to_cluster = dict(self.options)
        if force_nb_cpu_cores:
            options_to_pass_to_cluster['nb_core'] = force_nb_cpu_cores

        if self.pyNLoop_options['parallelization'] == 'cluster':
            cluster_name = options_to_pass_to_cluster['cluster_type']
            try:
                return cluster.from_name[cluster_name](**options_to_pass_to_cluster)
            except KeyError:
                # Check if a plugin define this type of cluster
                # check for PLUGIN format
                cluster_class = misc.from_plugin_import(self.plugin_path, 
                            'new_cluster', cluster_name,
                            info = 'cluster handling will be done with PLUGIN: %{plug}s' )
                if cluster_class:
                    return cluster_class(**options_to_pass_to_cluster)
        
        if self.pyNLoop_options['parallelization'] == 'multicore':
            try:
                import multiprocessing
                if not options_to_pass_to_cluster['nb_core']:
                    options_to_pass_to_cluster['nb_core'] = multiprocessing.cpu_count()
                logger.info('Using %d CPU cores' % options_to_pass_to_cluster['nb_core'])
            except ImportError:
                options_to_pass_to_cluster['nb_core'] = 1
                logger.warning('Impossible to detect the number of cores, therefore using one only.\n'+
                        'Use set nb_core X in order to set this number and be able to'+
                        'run in multicore.')

            return cluster.MultiCore(**options_to_pass_to_cluster)
        
    def parse_integrate_loop_options(self, args):
        """ Parsing arguments/options passed to the command integrate_loop."""

        options = { 
            'PS_point'              : 'random',
            'seed'                  : None,
            'sqrt_s'                : 1000.0,
            'target_accuracy'       : 1.0e-3,
            'batch_size'            : 1000,
            'verbosity'             : 1,
            'nb_CPU_cores'          : None,
            'phase_computed'        : 'Real',
            'survey_n_iterations'   : 10,
            'survey_n_points'       : 2000,
            'refine_n_iterations'   : 10,
            'refine_n_points'       : 1000
        }
        
        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args)>0:
                opt_args[-1] += ' %s'%arg
            else:
                new_args.append(arg)
        
        for arg in opt_args:
            try:
                key, value = arg.split('=')
            except:
                key, value = arg, None
            key = key[2:]
   
            if key == 'PS_point':
                if value not in ['random',]:
                    raise pyNLoopInvalidCmd('For now, pyNLoops only support the '+
                                                      'specification of a random PS point')
                options[key] = value

            elif key in ['seed','batch_size','verbosity','nb_CPU_cores', 
                         'survey_n_iterations','survey_n_points','refine_n_iterations','refine_n_points']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s'%(key,value))
                options[key] = parsed_int

            elif key in ['sqrt_s','target_accuracy']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s'%(key,value))
                options[key] = parsed_float

            elif key == 'phase_computed':
                if value.upper() in ['REAL','R','RE']:
                    options[key] == 'Real'
                elif value.upper() in ['IMAGINARY','I','IM']:
                    options[key] == 'Imaginary'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command integrated loop: %s"%key+
                      "\nAvailable commands are: %s"%(' '.join('--%s'%k for k in options)))

        return new_args, options

    def do_integrate_loop(self, line):
        """ Command for starting the numerical integration of a loop Feynman diagram."""
        
        args = self.split_arg(line)
        args, options = self.parse_integrate_loop_options(args)
        
        
        # For debugging you can easily print out the options as follows:
        #misc.sprint(options)
        
        if options['seed']:
            random.seed(options['seed'])

        chosen_topology_name = args[0]
        if chosen_topology_name not in self._hardcoded_topologies:
            raise InvalidCmd('For now, pyNLoop only support the specification of the'+
            ' following hardcoded topologies: %s'%(','.join(self._hardcoded_topologies.keys())))
        
        chosen_topology = self._hardcoded_topologies[chosen_topology_name]
        
        # Now generate the external momenta randomly, as this is the only mode we support
        # for now.
        phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            chosen_topology['masses'][:2], chosen_topology['masses'][2:],
            [options['sqrt_s']/2.,options['sqrt_s']/2.],
            beam_types=(1,1)
        )
        
        # Specifying None to get a random PS point
        random_PS_point, wgt, x1, x2 = phase_space_generator.get_PS_point(None)
        # Use the dictionary representation of the PS Point
        random_PS_point = random_PS_point.to_dict()
        
        # For loop calculations, it customary to consider all legs outgoing, so we must
        # flip the initial states directions.
        random_PS_point[1] = -random_PS_point[1]
        random_PS_point[2] = -random_PS_point[2]      
        # For debugging you can easily print out the chosen PS point as follows:
        #misc.sprint(str(random_PS_point))
   
        n_loop_integrand = chosen_topology['class'](
            n_loops            =   chosen_topology['n_loops'],
            external_momenta   =   random_PS_point,
            phase_computed     =   options['phase_computed'],
            # Data-structure for specifying a topology to be determined 
            topology           =   chosen_topology_name,
        )

        cluster = self.get_cluster(force_nb_cpu_cores=options['nb_CPU_cores'])
        integration_options = {
            'cluster'               :   cluster,
            'target_accuracy'       :   options['target_accuracy'],
            'batch_size'            :   options['batch_size'],
            'verbosity'             :   options['verbosity'],
            'survey_n_iterations'   :   options['survey_n_iterations'],
            'survey_n_points'       :   options['survey_n_points'],
            'refine_n_iterations'   :   options['refine_n_iterations'],
            'refine_n_points'       :   options['refine_n_points']
        }
        # Launch the integration
        if self.pyNLoop_options['integrator']=='Vegas3':
            integrator = Vegas3Integrator(n_loop_integrand, **integration_options)
        else:
            raise NotImplementedError
        
        # We are now finally ready to integrate :)
        logger.info("="*100)        
        logger.info('{:^100}'.format("Starting integration, lay down and enjoy..."),'$MG:color:GREEN')
        logger.info("="*100)

        xsec, error = integrator.integrate()

        run_output_path = MG5DIR
        if MPI_RANK==0:       
            logger.info("="*100)
            logger.info('{:^100}'.format("Loop integral result with integrator '%s':"%integrator.get_name()),'$MG:color:GREEN')
            logger.info('{:^100}'.format("%.5e +/- %.2e"%(xsec, error)),'$MG:color:BLUE')
            logger.info("="*100+"\n")
           
            # Write the result in 'cross_sections.dat' of the result directory
            xsec_summary = open(pjoin(run_output_path,'numerical_loop_integration.dat'),'w')
            xsec_summary_lines = []        
            xsec_summary_lines.append('%-30s%-30s%-30s'%('','Loop integral result','MC uncertainty'))
            xsec_summary_lines.append('%-30s%-30s%-30s'%('Total','%.8e'%xsec,'%.3e'%error))
            xsec_summary.write('\n'.join(xsec_summary_lines))
            xsec_summary.close()
        
    ######################################################################
    #
    # Example of the implementation of a trivial function 'hello_world'
    #
    ######################################################################

    def do_hello_world(self, line):
        """ Hello world command example."""

        logger.info('Hello World: %s'%line)

    def help_hello_world(self, line):
        """ Hello world command example."""

        logger.info('Contextual help for command hello world.')

    def complete_hello_world(self, text, line, begidx, endidx):
        """ Hello world command example."""

        return self.list_completion(text,['something', 'else'], line)

    #command to change the prompt 
    def preloop(self, *args, **opts):
        """only change the prompt after calling  the mother preloop command"""

        # The colored prompt screws up the terminal for some reason.
        #self.prompt = '\033[92mGGVV > \033[0m'
        self.prompt = 'pyNLoop > '

        # By default, load the UFO Standard Model
        logger.info("Loading default model for pyNLoop: sm")
        self.exec_cmd('import model sm', printcmd=False, precmd=True)

        # preloop mother
        madgraph_interface.CmdExtended.preloop(self)
