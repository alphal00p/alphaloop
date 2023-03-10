#####################################################
#                                                   #
#  Source file of the pyNLoop MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruijl              #
#                                                   #
#####################################################

import os
import logging
import sys
import shutil
import re
import random
import sympy
import math
import timeit
import functools
import copy

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from distutils.version import LooseVersion, StrictVersion

plugin_path = os.path.dirname(os.path.realpath( __file__ ))
import madgraph

import multiprocessing

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR, ReadWrite, MPI_RANK, MPI_SIZE, MPI_ACTIVE
import madgraph.interface.extended_cmd as cmd
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
import madgraph.various.cluster as cluster
import madgraph.integrator.vectors as vectors
import madgraph.various.progressbar as pbar

from madgraph.iolibs.files import cp, ln, mv

import madgraph.integrator.phase_space_generators as phase_space_generators
from madgraph.integrator.vegas3_integrator import Vegas3Integrator
from madgraph.integrator.pyCubaIntegrator import pyCubaIntegrator
# It will be nice to experiment with Thomas Hahn's integrators too later
# import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator

import nloop_integrands
import loop_momenta_generator

logger = logging.getLogger('pyNLoop.Interface')
try:
    import pysecdec_integrator
except ImportError:
    logger.warning("Import of pySecDec fails; you will not be able integrate loops depending on this integrator.")

pjoin = os.path.join
template_dir = pjoin(plugin_path, 'Templates')

utf8_exp_chars = {
            1 : u'\u00b9',
            2 : u'\u00b2',
            3 : u'\u00b3',
            4 : u'\u2074',
            5 : u'\u2075',
            6 : u'\u2076',
            7 : u'\u2077',
            8 : u'\u2078',
            9 : u'\u2079',
            0 : u'\u2070',
            '-' : u'\u207B',
            }
EPSILONS= { -i : (u'\u03B5\u207B%s'%utf8_exp_chars[i]).encode('utf-8') for i in range(1,9) }
EPSILONS.update({ i : (u'\u03B5%s'%utf8_exp_chars[i]).encode('utf-8') for i in range(1,9) })


class pyNLoopInterfaceError(MadGraph5Error):
    """ Error for the pyNLoop plugin """
    pass

class pyNLoopInvalidCmd(InvalidCmd):
    """ Invalid command issued to the pyNLoop interface. """
    pass

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
        'box1L_offshell_massless' :
            {   'class'     :   nloop_integrands.box1L_offshell_massless,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [100., 200., 300., 400.]
            },
        'box1L' : 
            {   'class'     :   nloop_integrands.box1L,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [100., 0., 0., 0.]
            },
        'box1L_onshell_massless' : 
            {   'class'     :   nloop_integrands.box1L_onshell_massless,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [0., 0., 0., 0.]
            },
        'box1L_subtracted' : 
            {   'class'     :   nloop_integrands.box1L_subtracted,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [0., 0., 0., 0.]
            },
        'box1L_subtracted_VH' : 
            {   'class'     :   nloop_integrands.box1L_subtracted_VH,
                'n_loops'   :   1,
                'n_legs'    :   4,
                'masses'    :   [0., 0., 0., 0.]
            },
        'box1L_direct_integration':
            {'class': nloop_integrands.box1L_direct_integration,
             'n_loops'  : 1,
             'n_legs'   : 4,
             'masses'   : [100., 200., 300., 400.]
             },
        'box1L_direct_integration_subtracted':
            {'class'    : nloop_integrands.box1L_direct_integration_subtracted,
             'n_loops'  : 1,
             'n_legs'   : 4,
             'masses'   : [0., 0., 0., 0.]
             },
        'box1L_direct_integration_one_offshell_subtracted':
            {'class'  : nloop_integrands.box1L_direct_integration_subtracted,
             'n_loops': 1,
             'n_legs' : 4,
             'masses' : [100., 0., 0., 0.]
             },
    }
    
    def validate_pySecDec_installation(self, pysecdec_path):
        
        necessary_lib_files = ['libcuba.a','libgsl.a','libgslcblas.a']
        for lib in necessary_lib_files:
            if not os.path.isfile(pjoin(pysecdec_path,'lib',lib)):
                return False
        necessary_bin_files = ['tform','form','gsl-config']
        for exe in necessary_bin_files:
            if not os.path.isfile(pjoin(pysecdec_path,'bin',exe)):
                return False
        return True
    
    def find_pySecDec_path(self):
        try:
            import pySecDec
        except ImportError:
            logger.warning(
                'pyNLoop could not successfully import the pySecDec module.\n'+
                'Make sure it is specified in your $PYTHONPATH.')
            return ''
        pySecDec_root_path = os.path.abspath(
                  pjoin(os.path.dirname(pySecDec.__file__), os.path.pardir,os.path.pardir))
        
        if LooseVersion(pySecDec.__version__) < LooseVersion("1.3.1"):
            raise pyNLoopInvalidCmd('Detected pySecDec version = %s but minimum required is 1.3.1')
        
        if not self.validate_pySecDec_installation(pySecDec_root_path):
            raise pyNLoopInvalidCmd('PySecDec installation appears to be non-standard (i.e. standalone).')

        
        return pySecDec_root_path
        
    def __init__(self, *args, **opts):
        """ Define attributes of this class."""
        
        self.pyNLoop_options = {
            # At first, we will limit ourselves to using a built-in distribution of pySecDec.
            'pySecDec_path'     :   self.find_pySecDec_path(),
            'integrator'        :   'auto',
            'parallelization'   :   'multicore'
        }
        # Generic option for all integrators
        self.integrator_options = {
            # For Vegas3
            'target_accuracy'       : 1.0e-3,
            'batch_size'            : 1000,
            'verbosity'             : 1,
            'survey_n_iterations'   : 10,
            'survey_n_points'       : 10000,
            'refine_n_iterations'   : 1,
            'refine_n_points'       : 100000,
            # For Cuba, see pyCubaIntegrator constructor for the list.
        }

        super(pyNLoopInterface, self).__init__(*args, **opts)

        # Temporary force pySecDec dependencies to be used
        os.environ['PATH']= os.environ['PATH']+':'+pjoin(self.pyNLoop_options['pySecDec_path'], 'bin')

        if self.options['nb_core'] is None:
            n_cpus = multiprocessing.cpu_count()
            logger.info("The Madgraph5_aMC@NLO option 'nb_core' is not set. PyNLoop automatically sets this option to %d."%n_cpus)
            self.do_set('nb_core %d'%n_cpus)

    def parse_set_pyNLoop_option(self, args):
        """ Parsing arguments/options passed to the command set_pyNLoop option."""

        options = { }

        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args) > 0:
                opt_args[-1] += ' %s' % arg
            else:
                new_args.append(arg)

        for arg in opt_args:
            try:
                key, value = arg.split('=')
            except:
                key, value = arg, None
            key = key[2:]

            # All options are declared valid in this contex
            options[key] = eval(str(value))

        return new_args, options

    def do_display_pyNLoop_option(self, line):
        """ Display pyNLoop options"""
        logger.info('%sGeneral pyNLoop options%s'%(misc.bcolors.GREEN, misc.bcolors.ENDC))
        logger.info('%s-----------------------%s'%(misc.bcolors.GREEN, misc.bcolors.ENDC))
        for opt in sorted(self.pyNLoop_options.keys()):
            logger.info('%-30s : %s'%(opt, str(self.pyNLoop_options[opt])))
        logger.info('%sIntegrator options%s'%(misc.bcolors.GREEN, misc.bcolors.ENDC))
        logger.info('%s------------------%s'%(misc.bcolors.GREEN, misc.bcolors.ENDC))
        for opt in sorted(self.integrator_options.keys()):
            logger.info('%-30s : %s'%(opt, str(self.integrator_options[opt])))

    def do_set_pyNLoop_option(self, line):
        """ Logic for setting pyNLoop options."""
        args = self.split_arg(line)
        args, options = self.parse_set_pyNLoop_option(args)
        key, value = args[:2]

        if key == 'integrator':
            if value not in ['Vegas3','pySecDec','Cuba','auto']:
                raise pyNLoopInvalidCmd(
"pyNLoop only supports Vegas3 and pySecDec integrator for now (or automatically set with value 'auto').")
            self.pyNLoop_options['integrator'] = value
            for opt, opt_value in options.items():
                if opt=='nb_core':
                    self.do_set('nb_core %d'%opt_value)
                else:
                    self.integrator_options[opt] = opt_value

        elif key == 'parallelization':
            if value not in ['cluster', 'multicore']:
                raise pyNLoopInvalidCmd("pyNLoop only supports parallelization "+
                                                        "modes 'cluser' and 'multicore'.")
            self.pyNLoop_options['parallelization'] = value
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
                #logger.info('Using %d CPU cores' % options_to_pass_to_cluster['nb_core'])
            except ImportError:
                options_to_pass_to_cluster['nb_core'] = 1
                logger.warning('Impossible to detect the number of cores, therefore using one only.\n'+
                        'Use set nb_core X in order to set this number and be able to'+
                        'run in multicore.')

            return cluster.MultiCore(**options_to_pass_to_cluster)

    def do_plot_deformation(self, line):
        """ Command for plotting the deformation along a given axis for a loop Feynman diagram."""

        def map_to_infinity(x):
            return ((1./(1.-x)) - 1./x)
        def map_from_infinity(x):
            return -2./(-2.+x-math.sqrt(4.+x**2))
        def find_offshell_scaling(rv, q_i, mass):
            """ Find the value of the scaling of the reference vector q_i which sends this propagator (q_i, mass) onshell."""

            m_sq = mass**2
            q_i_vec_sq = q_i[1]**2 + q_i[2]**2 + q_i[3]**2
            rv_vec_sq = rv[1]**2 + rv[2]**2 + rv[3]**2
            q_i_vec_dot_rv_vec = q_i[1]*rv[1] + q_i[2]*rv[2] + q_i[3]*rv[3]
            delta = rv[0]**2*(m_sq+q_i_vec_sq) - 2.*rv[0]*q_i[0]*q_i_vec_dot_rv_vec - \
                                                           rv_vec_sq*(m_sq-q_i[0]**2+q_i_vec_sq) + q_i_vec_dot_rv_vec**2
            if delta < 0.:
                return []

            return [
                (rv[0]*q_i[0] - q_i_vec_dot_rv_vec + math.sqrt(delta))/(rv[0]**2-rv_vec_sq),
                (rv[0]*q_i[0] - q_i_vec_dot_rv_vec - math.sqrt(delta))/(rv[0]**2-rv_vec_sq)
            ]

        ###########################################
        # Generate overall quantities for the test
        ###########################################

        args = self.split_arg(line)
        args, options = self.parse_plot_deformation_options(args)

        if options['seed']:
            random.seed(options['seed'])

        # Some sanity checks on the inputs
        if options['test']:
            if options['integration_space']:
                logger.warning('The deformation test is always performed directly in the physical loop momentum space.')
                options['integration_space'] = False

        if options['integration_space']:
            if any(item in options['items_to_plot'] for item in ['p','poles']):
                logger.warning('Poles will not be rendered when requiring to probe directly the space of integration variables.')
                options['items_to_plot'] = tuple( item for item in options['items_to_plot'] if item not in ['poles','p'])

        chosen_topology_name = args[0]
        if chosen_topology_name not in self._hardcoded_topologies:
            raise InvalidCmd('For now, pyNLoop only support the specification of the' +
                             ' following hardcoded topologies: %s' % (','.join(self._hardcoded_topologies.keys())))

        chosen_topology = self._hardcoded_topologies[chosen_topology_name]

        # Now generate the external momenta randomly, as this is the only mode we support
        # for now.
        phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            chosen_topology['masses'][:2], chosen_topology['masses'][2:],
            [options['sqrt_s'] / 2., options['sqrt_s'] / 2.],
            beam_types=(1, 1)
        )

        n_points = options['n_points']

        random_PS_point = None
        all_n_loop_integrands = []
        n_tests_done    = 0
        first           = True
        if options['test']:
            widgets = ["Performing deformation test :",
                       pbar.Percentage(), ' ', pbar.Bar(), ' ', pbar.ETA(), ' ',
                       pbar.Counter(format='%(value)d/%(max_value)d'), ' ']
            progress_bar = pbar.ProgressBar(widgets=widgets, maxval=n_points, fd=sys.stdout)
#            progress_bar = None
            if progress_bar:
                progress_bar.start()
        else:
            progress_bar = None

        # Now loop over many tests to be performed with random ref vectors and loop external kinematics
        while (first or (n_tests_done < n_points and options['test'])):
            first = False
            if progress_bar:
                progress_bar.update(n_tests_done)
            n_tests_done += 1
            if random_PS_point is None or (not options['keep_PS_point_fixed']):
                # Specifying None to get a random PS point
                random_PS_point, wgt, x1, x2 = phase_space_generator.get_PS_point(None)
                # Use the dictionary representation of the PS Point
                random_PS_point = random_PS_point.to_dict()

                # For loop calculations, it customary to consider all legs outgoing, so we must
                # flip the initial states directions.
                for i in random_PS_point:
                    if i > 2:
                        continue
                    random_PS_point[i] = -random_PS_point[i]
                if not options['test'] or options['keep_PS_point_fixed']:
                    logger.info('PS point considered:\n%s' % str(random_PS_point))

            # take a random vector
            if isinstance(options['reference_vector'], str) and options['reference_vector']=='random':
                ref_vec  = vectors.LorentzVector([random.random(),random.random(),random.random(),random.random()])
            else:
                ref_vec  = options['reference_vector']
            if isinstance(options['offset_vector'], str) and options['offset_vector']=='random':
                offset_vec  = vectors.LorentzVector([random.random(),random.random(),random.random(),random.random()])
            else:
                offset_vec  = options['offset_vector']
            # Make sure that the reference vector has each component normalized to its maximum one
            # This is to emphasize that its norm is irrelevant and only its direction matters, this is because
            # we will probe the region from -inf to +inf along this direction when working in the direct loop momentum
            # space
            ref_vec /= max(abs(k) for k in ref_vec)
            # Rescale the reference vector with sqrt_s if it does not live in the integration variable space
            if not options['integration_space']:
                ref_vec *= options['sqrt_s']
                offset_vec *= options['sqrt_s']
            else:
                # When working directly in the integration variable space, make sure that the offset vector also has
                # a maximum value for each of its components that is at most one.
                offset_vec /= max(max(abs(k) for k in offset_vec),1.0)

            if not options['test']:
                logger.info('Reference vector used: %s'%str(ref_vec))
                logger.info('Offset_vec vector used: %s'%str(offset_vec))

            # For debugging you can easily print out the options as follows:
            #misc.sprint(options)

            if len(all_n_loop_integrands)==0 or (not options['keep_PS_point_fixed']):
                all_n_loop_integrands = []
                for loop_momenta_generator_class, loop_momenta_generator_options in options['loop_momenta_generators']:
                    integrand_options = {
                        'n_loops'                       : chosen_topology['n_loops'],
                        'external_momenta'              : random_PS_point,
                        # Data-structure for specifying a topology to be determined
                        'topology'                       : chosen_topology_name,
                        'loop_momenta_generator_class'   : loop_momenta_generator_class,
                        'loop_momenta_generator_options' : loop_momenta_generator_options,
                    }
                    for opt, value in options.items():
                        if opt=='loop_momenta_generators':
                            continue
                        integrand_options[opt] = value
                    all_n_loop_integrands.append( (
                        'default' if loop_momenta_generator_class is None
                                      else loop_momenta_generator_class.__name__+'@'+str(loop_momenta_generator_options),
                        chosen_topology['class'],
                        integrand_options
                    ))

            first_integrand = all_n_loop_integrands[0][1](**all_n_loop_integrands[0][2])
            lmg = first_integrand.loop_momentum_generator

            ###############################
            # Do timing tests if asked for
            ###############################
            if options['test_timing']:
                for i_integrand, (lm_generator_name, n_loop_integrand_class, n_loop_integrand_options) in enumerate(
                        all_n_loop_integrands):
                    # Recover or build the n_loop_integrand instance
                    if i_integrand == 0:
                        n_loop_integrand = first_integrand
                    else:
                        n_loop_integrand = n_loop_integrand_class(**n_loop_integrand_options)
                    if options['test_timing']=='deformation':
                        function_to_test =  lambda: n_loop_integrand.loop_momentum_generator.generate_loop_momenta([random.random() for _ in range(4)])
                    elif options['test_timing']=='integrand':
                        function_to_test = lambda: n_loop_integrand([random.random() for _ in range(4)], [])
                    else:
                        raise pyNLoopInterfaceError("The timing test can only be performed on integrand' or 'deformation'.")
                    performance = min(timeit.repeat(functools.partial(function_to_test),
                                                        number=options['n_points'], repeat=3))/float(options['n_points'])
                    logger.info('Performance of %s for loop momenta generator %s : %.3e s / deformation.'%(
                                                                  options['test_timing'],lm_generator_name,performance))
                return

            ################
            # Analyze poles
            ################
            if (options['test'] or any(item in options['items_to_plot'] for item in ['p','poles'])):
                # compute the positions of the poles on the ref vec line
                poles = []
                imaginary_part_on_poles = []
                for i_prop, q_i in enumerate(lmg.q_is):
                    mass_prop = first_integrand.loop_propagator_masses[i_prop]
                    # Now solve for the scaling value of the ref_vector that sends this propagator onshell
                    scales_onshell = find_offshell_scaling(ref_vec, q_i-offset_vec, mass_prop)
                    for scale_onshell in scales_onshell:
#                        misc.sprint((offset_vec+scale_onshell*ref_vec-q_i).square()-mass_prop**2)
                        dp = lmg.deform_loop_momenta([offset_vec+scale_onshell*ref_vec,])[0]
                        denominator = (dp-q_i).square()-mass_prop**2
                        imaginary_part_on_poles.append(denominator.imag)
                    # Map back this scale onto the unit domain
                    for scale_onshell in scales_onshell:
                        poles.append(map_from_infinity(scale_onshell))

#                misc.sprint(len(imaginary_part_on_poles))
                if len(imaginary_part_on_poles)>0:
#                    misc.sprint(imaginary_part_on_poles)
                    # Make sure imaginary part is of the dimensionality of the momentum (not squared).
                    imaginary_part_on_poles = [math.sqrt(abs(ipp))*(-1. if ipp < 0. else 1.)
                                                                                         for ipp in imaginary_part_on_poles]
                    # Normalize to the biggest imaginary part
                    normalization_ipp = max(abs(ipp) for ipp in imaginary_part_on_poles)
                    if options['normalization'] != 'None':
                        dampening_power = 0.5
                        imaginary_part_on_poles = [((abs(ipp)/normalization_ipp)**dampening_power)*(-1. if ipp < 0. else 1.)
                                                                                     for ipp in imaginary_part_on_poles]

                    if any(ipp<0. for ipp in imaginary_part_on_poles):
                        if not options['test']:
                                logger.warning('The deformation leads to a negative imaginary part on some poles! : %s'%imaginary_part_on_poles)
                        else:
                            logger.critical('A test on the deformation %s failed!'%(all_n_loop_integrands[0][0]))
                            logger.critical('Imaginary part of the onshell propagators for that point: %s'%str(imaginary_part_on_poles))
                            logger.critical('This was the #%d test with seed %d.'%(n_tests_done, options['seed']))
                            logger.critical('External PS point:\n %s'%str(random_PS_point))
                            logger.critical('Reference vector used: %s'%str(ref_vec))
                            logger.critical('Offset vector used: %s'%str(offset_vec))
                            return

#                   misc.sprint("Poles mapped: %s"%poles)
                    if not options['test']:
                        for p, ipp in zip(poles,imaginary_part_on_poles):
                            y_location = 0. if not options['log_axis'] else 1.
                            plt.plot(p, y_location, marker='o', markersize=3, color="red")
                            plt.plot([p, p], [y_location, y_location+ipp], 'k-', lw=2, color="red")

        if progress_bar:
            progress_bar.finish()
        if options['test']:
            logger.info('A total of %d tests on the deformation %s were successfully performed.'%(
                                                                       options['n_points'], all_n_loop_integrands[0][0]))
            return

        #################################
        # Now generate the plotting data
        #################################
        scaling_range = options['range']
        if options['scale_progression']=='linear':
            x_entries = [scaling_range[0]+ (i / float(n_points))*(scaling_range[1]-scaling_range[0]) for i in range(1, n_points)]
        elif options['scale_progression']=='exponential':
            # Scale exponentially dense close towards scaling_range[0]
            x_entries = [scaling_range[0]+ 0.5*(10.**(-10.*(i/float(n_points//2))))*(scaling_range[1]-scaling_range[0]) for i in range(n_points//2,0,-1)]
            # Scale exponentially dense close towards scaling_range[1]
            x_entries.extend([scaling_range[0]+ 0.5*(2.-10.**(-10.*(i/float(n_points//2))))*(scaling_range[1]-scaling_range[0]) for i in range(0, n_points//2-1)])
        else:
            raise pyNLoopInterfaceError("Value '%s' for scaling progression not supported.")

        all_deformed_points = []
        all_jacobians       = []
        for i_integrand, (lm_generator_name, n_loop_integrand_class, n_loop_integrand_options) in enumerate(all_n_loop_integrands):
            # Recover or build the n_loop_integrand instance
            if i_integrand==0:
                n_loop_integrand = first_integrand
            else:
                n_loop_integrand = n_loop_integrand_class(**n_loop_integrand_options)

            # List the jacobians for all points for this integrand
            jacobians       = [1.0,]*len(x_entries)

            # now sample n_points points and compute the deformation
            if options['integration_space']:
                # We will probe the integration space within the unit hyperbox following this line:
                #          y_vec = offset_vec + scaling * ref_vec
                # And we want the parameter 'scaling' to be in [0,1] with 0 corresponding to one intersection of this
                # line with the unit hyperbox (call it point 'A') and 1 corresponding to the other intersection with the
                # unit hyperbox (call it point 'B'). We therefore proceed here to computing the min and max value of this
                # scaling so as to remap it to the interval [0,1].
                min_scale, max_scale = None, None
                for o, r in zip(offset_vec, ref_vec):
                    # Intersections with the planes defined by some y_vec component being 0. or 1.
                    if r==0.:
                        continue
                    intersections = sorted( [(component_value - o) / r for component_value in [0., 1.]] )
                    if min_scale is None or intersections[0] > min_scale:
                        min_scale = intersections[0]
                    if max_scale is None or intersections[1] < max_scale:
                        max_scale = intersections[1]

                if min_scale is None or max_scale is None or max_scale-min_scale <= 0.:
                    raise pyNLoopInterfaceError('Could not compute the scaling to apply in the integration space.')
                points_in_integration_space = [offset_vec + ref_vec * ( min_scale + x* (max_scale - min_scale)) for x in x_entries]

                # Then use the conformal map of the loop momentum generator of the integrand
                points = []
                for i, p in enumerate(points_in_integration_space):
                    #misc.sprint('Point before: %s'%str(p))
                    loop_momenta, conformal_jac = n_loop_integrand.loop_momentum_generator.map_to_infinite_hyperbox(p)
                    #misc.sprint('Point after: %s, %f'%(str(loop_momenta[0]),conformal_jac))
                    jacobians[i] *= conformal_jac
                    points.append(list(loop_momenta[0]))
                    if any(math.isnan(c) or math.isinf(c) for c in loop_momenta[0]):
                        raise pyNLoopInterfaceError('The conformal map failed to correctly map the following input point:\n'+
                                                                        str(p)+'\nIt returned:\n%s'%str(loop_momenta[0]))
            else:
                points = [offset_vec + ref_vec * map_to_infinity(x) for x in x_entries]

            normalizations = [1., ] * n_points
            if options['normalization'] == 'distance_real':
                normalizations = [math.sqrt(sum(k_i ** 2 for k_i in p)) for p in points]

            # Make sure the normalization is capped to be minimum 1.0
            normalizations = [max(n, 1.0e-99) for n in normalizations]

            deformed_points = []
            widgets = ["Loop deformation for generator %s :"%lm_generator_name.split('@')[0],
                pbar.Percentage(), ' ', pbar.Bar(),' ', pbar.ETA(), ' ',
                pbar.Counter(format='%(value)d/%(max_value)d'), ' ']
            progress_bar = pbar.ProgressBar(widgets=widgets, maxval=len(points), fd=sys.stdout)
#           progress_bar = None
            if progress_bar:
                progress_bar.start()
            for i_point, p in enumerate(points):
                res = n_loop_integrand.loop_momentum_generator.apply_deformation([p, ])
                deformed_points.append(res[0][0])
                #misc.sprint('Deformation jacobian:', res[1])
                jacobians[i_point] *= res[1]
                if progress_bar:
                    progress_bar.update(i_point)
            if progress_bar:
                progress_bar.finish()
            all_deformed_points.append(deformed_points)
            all_jacobians.append(jacobians)

#           misc.sprint("Last 5 points: %s"%(',\n'.join([str(p) for p in points[-5:]])))
#           misc.sprint("Last 5 deformed points: %s"%(',\n'.join([str(p) for p in deformed_points[-5:]])))

            # Note that the deformation should not have changed the real components
            # get the distance on the imaginary axis.
            for i_comp in range(4):
                if i_comp not in options['items_to_plot']:
                    continue
                plt.plot(x_entries, [d[i_comp].imag/normalizations[i_point] for i_point, d in enumerate(deformed_points)],
                         label='component #%d @%s'%(i_comp,lm_generator_name))

            if any(item in options['items_to_plot'] for item in ['d','distance']):
                dist = [ math.sqrt(sum(di.imag**2 for di in d))/normalizations[i_point] for i_point, d in enumerate(deformed_points)]
                plt.plot(x_entries, dist, linewidth=2.0, label='distance @%s'%lm_generator_name)

            # Also add the integrand (Any integrand would do for the evaluation since the deformation should be irrelevant,
            #  so we pick the first here)
            if any(item in options['items_to_plot'] for item in ['integrand_real','integrand_imag','integrand_abs']):
                integrands = [first_integrand([dp,], [], input_already_in_infinite_hyperbox=True, jacobian=jac, phase='All' )
                                                                           for dp,jac in zip(deformed_points,jacobians)]
                # Normalize to the largest value of the integrand encountered
                if 'integrand_real' in options['items_to_plot']:
                    integrand_reals = [itg.real for itg in integrands]
                    max_integrand_real = max(integrand_reals)
                    if max_integrand_real > 0.:
                        integrand_reals = [itg/max_integrand_real for itg in integrand_reals]
                    plt.plot(x_entries,integrand_reals, label='integrand_real @%s'%lm_generator_name)
                if 'integrand_imag' in options['items_to_plot']:
                    integrand_imags = [itg.imag for itg in integrands]
                    max_integrand_imag = max(integrand_imags)
                    if max_integrand_imag > 0.:
                        integrand_imags = [itg/max_integrand_imag for itg in integrand_imags]
                    plt.plot(x_entries,integrand_imags, label='integrand_imag @%s'%lm_generator_name)
                if 'integrand_abs' in options['items_to_plot']:
                    integrand_abss = [abs(itg) for itg in integrands]
                    max_integrand_abss = max(integrand_abss)
                    if max_integrand_abss > 0.:
                        integrand_abss = [itg/max_integrand_abss for itg in integrand_abss]
                    plt.plot(x_entries,integrand_abss, label='integrand_abs @%s'%lm_generator_name)
            if any(item in options['items_to_plot'] for item in [
                                               'integrand_no_jac_real','integrand_no_jac_imag','integrand_no_jac_abs']):
                integrands = [first_integrand([dp,], [], input_already_in_infinite_hyperbox=True, jacobian=1.0, phase='All' )
                                                                           for dp,jac in zip(deformed_points,jacobians)]
                # Normalize to the largest value of the integrand encountered
                if 'integrand_no_jac_real' in options['items_to_plot']:
                    integrand_reals = [itg.real for itg in integrands]
                    max_integrand_real = max(integrand_reals)
                    if integrand_reals > 0.:
                        integrand_reals = [itg/max_integrand_real for itg in integrand_reals]
                    plt.plot(x_entries,integrand_reals, label='integrand_no_jac_real @%s'%lm_generator_name)
                if 'integrand_no_jac_imag' in options['items_to_plot']:
                    integrand_imags = [itg.imag for itg in integrands]
                    max_integrand_imag = max(integrand_imags)
                    if max_integrand_imag > 0.:
                        integrand_imags = [itg/max_integrand_imag for itg in integrand_imags]
                    plt.plot(x_entries,integrand_imags, label='integrand_no_jac_imag @%s'%lm_generator_name)
                if 'integrand_no_jac_abs' in options['items_to_plot']:
                    integrand_abss = [abs(itg) for itg in integrands]
                    max_integrand_abss = max(integrand_abss)
                    if max_integrand_abss > 0.:
                        integrand_abss = [itg/max_integrand_abss for itg in integrand_abss]
                    plt.plot(x_entries,integrand_abss, label='integrand_no_jac_abs @%s'%lm_generator_name)
            if any(item in options['items_to_plot'] for item in ['jac_real','jac_imag','jac_abs']):
                # Normalize to the largest value of the integrand encountered
                if 'jac_real' in options['items_to_plot']:
                    jacobian_reals = [jac.real for jac in jacobians]
                    max_jacobian_real = max(jacobian_reals)
                    if max_jacobian_real>0.:
                        jacobian_reals = [jac/max_jacobian_real for jac in jacobian_reals]
                    plt.plot(x_entries,jacobian_reals, label='jac_real @%s'%lm_generator_name)
                if 'jac_imag' in options['items_to_plot']:
                    jacobian_imags = [jac.imag for jac in jacobians]
                    max_jacobian_imag = max(jacobian_imags)
                    if max_jacobian_imag > 0.:
                        jacobian_imags = [jac/max_jacobian_imag for jac in jacobian_imags]
                    plt.plot(x_entries,jacobian_imags, label='jac_imag @%s'%lm_generator_name)
                if 'jac_abs' in options['items_to_plot']:
                    jacobian_abss = [abs(jac) for jac in jacobians]
                    max_jacobian_abss = max(jacobian_abss)
                    if max_jacobian_abss > 0.:
                        jacobian_abss = [jac/max_jacobian_abss for jac in jacobian_abss]
                    plt.plot(x_entries,jacobian_abss, label='jac_abs @%s'%lm_generator_name)

        # Plot relative difference of deformation, normalized to its max
        if len(all_n_loop_integrands)>1 and any((isinstance(item, str) and item.startswith('difference')) for item
                                                                                           in options['items_to_plot']):
            relative_deformation_differences = []
            relative_jacobian_differences    = []
            n_lmg = len(all_n_loop_integrands)
            for i_point in range(len(all_deformed_points[0])):
                if 'difference_deformation' in options['items_to_plot']:
                    data_for_deformation_difference = []
                    for i_component in range(4):
                        data_for_deformation_difference.append(
                            [all_deformed_points[i_lmg][i_point][i_component].real for i_lmg in range(n_lmg) ] )
                        data_for_deformation_difference.append(
                            [all_deformed_points[i_lmg][i_point][i_component].imag for i_lmg in range(n_lmg)] )
                    if options['normalization']:
                        relative_deformation_differences.append(max(
                            [(max(k) - min(k))/(0.5*((abs(max(k))+abs(min(k))) if (abs(max(k))+abs(min(k)))>0. else 1.))
                                                                            for k in data_for_deformation_difference ]))
                    else:
                        relative_deformation_differences.append(max(
                                                            [(max(k) - min(k)) for k in data_for_deformation_difference ]))
                if 'difference_jacobian' in options['items_to_plot']:
                    data_for_jacobian_difference = []
                    data_for_jacobian_difference.append([all_jacobians[i_lmg][i_point].real for i_lmg in range(n_lmg)])
                    data_for_jacobian_difference.append([all_jacobians[i_lmg][i_point].imag for i_lmg in range(n_lmg)])
                    if options['normalization']:
                        relative_jacobian_differences.append(max(
                            [(max(k) - min(k))/(0.5*((abs(max(k))+abs(min(k))) if (abs(max(k))+abs(min(k)))>0. else 1.))
                                                                            for k in data_for_jacobian_difference ]))
                    else:
                        relative_jacobian_differences.append(max(
                                                            [(max(k) - min(k)) for k in data_for_jacobian_difference ]))
            # Further normalize if asked for.
            if options['normalization']:
                if 'difference_deformation' in options['items_to_plot']:
                    max_diff = max(relative_deformation_differences)
                    if max_diff > 0.:
                        relative_deformation_differences = [d/max_diff for d in relative_deformation_differences]
                    plt.plot(x_entries, relative_deformation_differences, label='Rel. deform. diff. %s' % (
                                                                        '( x%.3e )'%max_diff if max_diff > 0. else 'ZERO' ))
                if 'difference_jacobian' in options['items_to_plot']:
                    max_diff = max(relative_jacobian_differences)
                    if max_diff > 0.:
                        relative_jacobian_differences = [d/max_diff for d in relative_jacobian_differences]
                    plt.plot(x_entries, relative_jacobian_differences, label='Rel. jacobian diff. %s' % (
                                                                        '( x%.3e )'%max_diff if max_diff > 0. else 'ZERO' ))
            else:
                if 'difference_deformation' in options['items_to_plot']:
                    plt.plot(x_entries, relative_deformation_differences, label='Deform. diff.')
                if 'difference_jacobian' in options['items_to_plot']:
                    plt.plot(x_entries, relative_jacobian_differences, label='Jacobian diff.')

        plt.xlim(options['range'])
        xlabel = r'$\lambda$'
        if options['integration_space']:
            xlabel += '@IntegrationVariablesSpace'
        xlabel += ' [%s progression]'%options['scale_progression']
        plt.xlabel(xlabel)
        if options['log_axis']:
            plt.semilogy()

        fontP = FontProperties()
        fontP.set_size('small')
        legend = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, prop=fontP)
        if options['save_plot'] not in ['',None]:
            extension = os.path.basename(options['save_plot']).split('.')
            if len(extension)==1:
                extension = None
            else:
                extension = extension[-1]
            if extension is None:
                extension = 'svg'
                options['save_plot'] += '.'+extension
            logger.info("Saving plot to '%s'."%options['save_plot'])
            if extension=='svg':
                plt.savefig(options['save_plot'], bbox_extra_artists=(legend,), bbox_inches='tight')
            else:
                plt.savefig(options['save_plot'], dpi=500, bbox_extra_artists=(legend,), bbox_inches='tight')
        if options['show_plot']:
            plt.tight_layout(rect=(0,0,1,0.75))
            plt.show()
        else:
            logger.info("Display of the plot skipped according to user's request.")

    def generate_dummy_phase_space_point(self, masses, sqrt_s):
        # Generate a dummy_PS_point just for the instantiation
        phase_space_generator = phase_space_generators.FlatInvertiblePhasespace(
            masses[:2], masses[2:],
            [sqrt_s / 2., sqrt_s / 2.],
            beam_types=(1, 1)
        )
        # Specifying None to get a random PS point just for the instantiation of the integrands.
        dummy_PS_point, _, _, _ = phase_space_generator.get_PS_point(None)
        # Use the dictionary representation of the PS Point
        dummy_PS_point = dummy_PS_point.to_dict()
        # For loop calculations, it customary to consider all legs outgoing, so we must
        # flip the initial states directions.
        for i in dummy_PS_point:
            if i > 2:
                continue
            dummy_PS_point[i] = -dummy_PS_point[i]

        return dummy_PS_point


    def do_plot_box_scattering_angle(self, line):
        """ Command for plotting the integral over the scattering angle theta, for a given phi."""

        args = self.split_arg(line)
        args, options = self.parse_plot_box_scattering_angle_options(args)

        if options['seed']:
            random.seed(options['seed'])

        if options['phase_computed'] not in ['Real','Imaginary']:
            raise InvalidCmd("Command 'do_plot_box_scattering_angle' only allows option 'phase_computed' to be 'Real' or 'Imaginary'.")

        if all((not options[opt]) for opt in ['compute_numerical_result','compute_analytic_result']):
            raise InvalidCmd("The command 'do_plot_box_scattering_angle' must be instructed to compute at least the numerical or analytical result.")

        chosen_topology_name = args[0]
        if chosen_topology_name not in self._hardcoded_topologies:
            raise InvalidCmd('For now, pyNLoop only support the specification of the' +
                             ' following hardcoded topologies: %s' % (','.join(self._hardcoded_topologies.keys())))

        chosen_topology = self._hardcoded_topologies[chosen_topology_name]
        M1 = chosen_topology['masses'][0]
        M2 = chosen_topology['masses'][1]
        M3 = chosen_topology['masses'][2]
        M4 = chosen_topology['masses'][3]

        n_points = options['n_points']

        # create progress bar
        widgets = ["Performing deformation test :",
                    pbar.Percentage(), ' ', pbar.Bar(), ' ', pbar.ETA(), ' ',
                    pbar.Counter(format='%(value)d/%(max_value)d'), ' ']

        if options['verbosity'] == 0:
            progress_bar = pbar.ProgressBar(widgets=widgets, maxval=n_points, fd=sys.stdout)
        else:
            progress_bar = None

        if progress_bar:
            progress_bar.start()

        channels = [0,1,2] if options['channel'] is None else [options['channel'],]

        # Instantiate an integrator
        loop_momenta_generator_class, loop_momenta_generator_options = options['loop_momenta_generator']
        integrand_options = {
            'n_loops': chosen_topology['n_loops'],
            'external_momenta': self.generate_dummy_phase_space_point(chosen_topology['masses'], options['sqrt_s']),
            # Data-structure for specifying a topology to be determined
            'topology': chosen_topology_name,
            'loop_momenta_generator_class': loop_momenta_generator_class,
            'loop_momenta_generator_options': loop_momenta_generator_options,
            'channel': None, # will be set later
        }
        for opt, value in options.items():
            if opt == 'loop_momenta_generator':
                continue
            integrand_options[opt] = value

        loop_integrand = chosen_topology['class'](**integrand_options)

        # Choose the integrator
        if self.pyNLoop_options['integrator'] == 'auto':
            integrator_name = loop_integrand._supported_integrators[0]
        else:
            if self.pyNLoop_options['integrator'] not in loop_integrand._supported_integrators:
                integrator_name = loop_integrand._supported_integrators[0]
                logger.warning(
                    "Specified integrator '%s' is not supported by integrand '%s'. Using '%s' instead." %
                    (self.pyNLoop_options['integrator'], loop_integrand.nice_string(), integrator_name))
            else:
                integrator_name = self.pyNLoop_options['integrator']

        # Now set the integrator options
        integrator_options = dict(self.integrator_options)
        for opt, value in options.items():
            if opt in integrator_options:
                integrator_options[opt] = value
        integrator_options['cluster'] = self.get_cluster(force_nb_cpu_cores=int(self.options['nb_core']))
        integrator_options['pySecDec_path'] = self.pyNLoop_options['pySecDec_path']
        if integrator_name == 'Vegas3':
            if options['phase_computed'] == 'All':
                logger.warning('Vegas3 integrator cannot simultaneously integrate the real and imaginary part of' +
                                " the loop for now. The user's choice 'All' for 'phase_computed' is reverted to 'Real'.")
                loop_integrand.phase_computed = 'Real'
            integrator = Vegas3Integrator(loop_integrand, **integrator_options)
        elif integrator_name == 'pySecDec':
            integrator = pysecdec_integrator.pySecDecIntegrator(loop_integrand, **integrator_options)
        elif integrator_name == 'Cuba':
            integrator = pyCubaIntegrator(loop_integrand, **integrator_options)
        else:
            raise pyNLoopInterfaceError("Integrator '%s' not implemented." % integrator_name)


        # Initialize a hook to OneLOop in order to gather the analytic results
        analytic_results_available = False
        if options['compute_analytic_result']:
            loop_integrand.setup_analytic_computation(self)
            if not loop_integrand.is_analytic_result_available():
                logger.warning('AVH OneLOop library could not be properly loaded. Analytic results will not be available.')
                analytic_results_available = False
            else:
                analytic_results_available = True

        # loop over the angle theta
        results, errors = [], []
        results_analytic = []

        res_summary = open(pjoin(options['output_folder'],'scattering_angle.dat'),'w')
        res_summary.write('Integral {} with channel {}:\n'.format(chosen_topology_name, options['channel']))
        res_summary.write('Theta\t\tPhi\t\t\t{}\t\t\t\t\tError\n'.format(options['phase_computed']))

        n_tests_done = 0
        theta_range = options['range']
        thetas = [ theta_range[0] + (theta_range[1] - theta_range[0]) * i / float(n_points) for i in range(n_points) ]
        while n_tests_done < n_points:
            theta = thetas[n_tests_done]

            if progress_bar:
                progress_bar.update(n_tests_done)
            n_tests_done += 1

            # generate the phase space point
            sqrtS = options['sqrt_s']
            phi = options['phi']

            pz = math.sqrt((sqrtS - M1 - M2)*(sqrtS + M1 - M2)*(sqrtS - M1 + M2)*(sqrtS + M1 + M2))/(2.*sqrtS)
            pf = math.sqrt((sqrtS - M3 - M4)*(sqrtS + M3 - M4)*(sqrtS - M3 + M4)*(sqrtS + M3 + M4))/(2.*sqrtS)
            ps1 = vectors.LorentzVector([math.sqrt(M1**2 + pz**2), 0, 0, pz])
            ps2 = vectors.LorentzVector([math.sqrt(M2**2 + pz**2), 0, 0, -pz])
            ps3 = vectors.LorentzVector([math.sqrt(M3**2 + pf**2), pf*math.cos(phi)*math.sin(theta), pf*math.sin(phi)*math.sin(theta), pf*math.cos(theta)])
            ps4 = vectors.LorentzVector([math.sqrt(M4**2 + pf**2), -(pf*math.cos(phi)*math.sin(theta)), -(pf*math.sin(phi)*math.sin(theta)), -(pf*math.cos(theta))])
            external_momenta = vectors.LorentzVectorList([ps1, ps2, ps3, ps4]).to_dict()

            # For loop calculations, it customary to consider all legs outgoing, so we must
            # flip the initial states directions.
            for i in external_momenta:
                if i > 2:
                    continue
                external_momenta[i] = -external_momenta[i]

            if options['verbosity'] > 0:
                logger.debug('External momentum configuration probed for angle theta=%f is:\n%s'%
                                                                                         (theta, str(external_momenta)))

            if options['compute_numerical_result']:
                summed_result = 0
                summed_error = 0
                for channel in channels:

                    # Assign the new kinematic configuration
                    loop_integrand.assign_kinematic_configuration(external_momenta, channel)
                    with misc.Silence(active = options['verbosity'] <= 1):
                        amplitude, error = integrator.integrate()

                    if options['verbosity'] > 0:
                        logger.debug('Amplitude for channel {}: {} +/- {}'.format(channel, amplitude, error))
                    summed_result += amplitude
                    summed_error += error**2
                results.append(summed_result)
                errors.append(math.sqrt(summed_error))

            # Now compute the analytical result for this point
            if analytic_results_available:
                analytic_res_real, analytic_res_real_imag = loop_integrand.get_analytic_result(external_momenta)
                if options['phase_computed'] == 'Real':
                    results_analytic.append(analytic_res_real[0])
                elif options['phase_computed'] == 'Imaginary':
                    results_analytic.append(analytic_res_real_imag[0])

            if options['verbosity'] > 0:
                logger.debug('Results obtained for the PS point above:\nAnalytic: %s\npyNLoop : %s\nMC_error: %s'%
                        ('%.16e'%results_analytic[-1] if analytic_results_available else 'N/A',
                         '%.16e'%results[-1] if options['compute_numerical_result'] else 'N/A',
                         '%.2e'%errors[-1] if options['compute_numerical_result'] else 'N/A' ))

            res_summary.write('{:.3f}\t\t{:.3f}\t\t{}\t\t{}\t\t{}\n'.format(theta, phi,
                    results[-1] if options['compute_numerical_result'] else 0.,
                    errors[-1] if options['compute_numerical_result'] else 0.,
                    results_analytic[-1] if analytic_results_available else 0.)
            )

        if progress_bar:
            progress_bar.finish()

        res_summary.close()

        logger.debug('Results:\nanalytic: %s\ncentral value: %s\nMC error: %s'%(results_analytic,results, errors))
        all_handles = []
        if options['compute_numerical_result']:
            all_handles.append(plt.errorbar(x=thetas, y=results, yerr=errors, label='pyNLoop'))
        if analytic_results_available:
            all_handles.append(plt.errorbar(x=thetas, y=results_analytic, yerr=[0.,]*len(results_analytic), label='oneLOop'))
        plt.legend(handles=all_handles)
        plt.xlabel('theta')
        plt.ylabel('Amplitude (%s part)'%options['phase_computed'])


        if options['channel'] is None:
            plt.title('{}'.format(chosen_topology_name))
        else:
            plt.title('Channel {} of {}'.format(options['channel'], chosen_topology_name))

        if options['log_axis']:
            plt.yscale('log')

        plt.show()

    def parse_plot_box_scattering_angle_options(self, args):
        """ Parsing arguments/options passed to the command plot_box_scattering_angle."""

        options = { 
            'seed'                  : None,
            'sqrt_s'                : 1000.0,
            'phi'                   : math.pi / 4.,
            'cpp_integrand'         : False,
            'verbosity'             : 1,
            'output_folder'         : pjoin(MG5DIR,'MyPyNLoop_output'),
            'loop_momenta_generator' : self.parse_lmgc_specification('default'),
            'phase_computed'        : 'Real',
            'range'                 : (0., math.pi * 0.99),
            'n_points'              : 10,
            'channel'               : None,
            'log_axis'              : False,
            'show_plot'             : True,
            'save_plot'             : '',
            'compute_analytic_result': True,
            'compute_numerical_result': True,
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

            if key in ['seed','verbosity','nb_CPU_cores', 'n_points', 'channel']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s'%(key,value))
                options[key] = parsed_int

            elif key in ['reference_vector', 'offset_vector']:
                if value.lower() in ['random','r']:
                    options['reference_vector'] = 'random'
                else:
                    try:
                        ref_vec = tuple(eval(value))
                    except:
                        raise pyNLoopInvalidCmd("Cannot parse reference vector specification: %s"%str(value))
                    if len(ref_vec)!=4:
                        raise pyNLoopInvalidCmd("Reference vector must be of length 4, not %d"%len(ref_vec))
                    options[key] = vectors.LorentzVector(ref_vec)

            elif key=='scale_progression':
                if value not in ['linear','exponential']:
                    raise pyNLoopInvalidCmd("The values for the option 'scale_progression' can only be in "+
                                                                             "['linear','exponential'], not %s." %value)
                else:
                    options[key] = value

            elif key=='range':
                try:
                    range = tuple(eval(value))
                except:
                    raise pyNLoopInvalidCmd("Cannot parse reference vector specification: %s"%str(value))
                if len(range)!=2:
                    raise pyNLoopInvalidCmd("Range must be a list of length 2, not %d"%len(range))
                options['range'] = range

            elif key=='items_to_plot':
                try:
                    comps = tuple(eval(value))
                except:
                    raise pyNLoopInvalidCmd("Cannot parse items to plot: %s"%value)
                options['items_to_plot'] = comps

            elif key in ['loop_momenta_generator_class', 'lmgc']:
                try:
                    lmgc_name = eval(value)
                except:
                    raise pyNLoopInvalidCmd("Cannot parse loop momentum generator class specification: %s"%str(value))
                options['loop_momenta_generator'] = self.parse_lmgc_specification(lmgc_name)

            elif key in ['sqrt_s', 'phi']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s'%(key,value))
                options[key] = parsed_float

            elif key in ['test_timing']:
                if value is None:
                    options[key] = 'deformation'
                elif value.lower() in ['integrand','deformation']:
                    options[key] = value.lower()
                elif value.lower() in ['f','false']:
                    options[key] = False
                else:
                    raise pyNLoopInvalidCmd("Option 'test_timing' can only take values in ['integrand','deformation','False'], not %s"%value)

            elif key in ['force','log_axis','test','keep_PS_point_fixed','show_plot','integration_space',
                         'cpp_integrand','compute_analytic_result','compute_numerical_result']:
                parsed_bool = (value is None) or (value.upper() in ['T','TRUE','ON','Y'])
                options[key] = parsed_bool

            elif key in ['output_folder','save_plot']:
                if os.path.isabs(value):
                    options[key] = value
                else:
                    options[key] = pjoin(MG5DIR,value)

            elif key in ['normalization']:
                if value is None:
                    value = 'distance_real'
                _normalization_modes_supported = ['distance_real',]
                if value not in _normalization_modes_supported:
                    raise pyNLoopInvalidCmd('Normalization modes supported are %s, not %s'%(_normalization_modes_supported, value))
                options['normalization'] = value

            elif key == 'phase_computed':
                if value.upper() in ['REAL', 'R', 'RE']:
                    options[key] = 'Real'
                elif value.upper() in ['IMAGINARY', 'I', 'IM']:
                    options[key] = 'Imaginary'
                elif value.upper() in ['ALL', 'A']:
                    options[key] = 'All'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command plot deformation: %s"%key+
                      "\nAvailable commands are: %s"%(' '.join('--%s'%k for k in options)))

        return new_args, options

    def parse_plot_deformation_options(self, args):
        """ Parsing arguments/options passed to the command plot_deformation."""

        options = { 
            'seed'                  : None,
            'sqrt_s'                : 1000.0,
            'cpp_integrand'         : False,
            'verbosity'             : 1,
            'output_folder'         : pjoin(MG5DIR,'MyPyNLoop_output'),
            'force'                 : False,
            'phase_computed'        : 'Real',
            'reference_vector'      : 'random',
            'offset_vector'         : 'random',
            'loop_momenta_generators' : [self.parse_lmgc_specification('default'),],
            'n_points'              : 100,
            'items_to_plot'         : (0,1,2,3,'distance','poles'),
            'range'                 : (0.,1.),
            'normalization'         : 'None',
            'log_axis'              : False,
            'channel'               : None,
            'integration_space'     : False,
            'scale_progression'     : 'linear',
            # When 'test' is on, a battery of tests is performed instead of the plotting of the deformation
            'test'                  : False,
            'test_timing'           : None,
            'keep_PS_point_fixed'   : False,
            'show_plot'             : True,
            'save_plot'             : ''
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

            if key in ['seed','verbosity','nb_CPU_cores', 'n_points', 'channel']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s'%(key,value))
                options[key] = parsed_int

            elif key in ['reference_vector', 'offset_vector']:
                if value.lower() in ['random','r']:
                    options['reference_vector'] = 'random'
                else:
                    try:
                        ref_vec = tuple(eval(value))
                    except:
                        raise pyNLoopInvalidCmd("Cannot parse reference vector specification: %s"%str(value))
                    if len(ref_vec)!=4:
                        raise pyNLoopInvalidCmd("Reference vector must be of length 4, not %d"%len(ref_vec))
                    options[key] = vectors.LorentzVector(ref_vec)

            elif key=='scale_progression':
                if value not in ['linear','exponential']:
                    raise pyNLoopInvalidCmd("The values for the option 'scale_progression' can only be in "+
                                                                             "['linear','exponential'], not %s." %value)
                else:
                    options[key] = value

            elif key=='range':
                try:
                    range = tuple(eval(value))
                except:
                    raise pyNLoopInvalidCmd("Cannot parse reference vector specification: %s"%str(value))
                if len(range)!=2:
                    raise pyNLoopInvalidCmd("Range must be a list of length 2, not %d"%len(range))
                options['range'] = range

            elif key=='items_to_plot':
                try:
                    comps = tuple(eval(value))
                except:
                    raise pyNLoopInvalidCmd("Cannot parse items to plot: %s"%value)
                options['items_to_plot'] = comps

            elif key in ['loop_momenta_generator_classes', 'lmgc']:
                try:
                    lmgc_names = list(eval(value))
                except:
                    raise pyNLoopInvalidCmd("Cannot parse loop momentum generator class specification: %s"%str(value))
                lmgcs = []
                for lmgc_specification in lmgc_names:
                    lmgcs.append(self.parse_lmgc_specification(lmgc_specification))
                options['loop_momenta_generators'] = lmgcs

            elif key in ['sqrt_s']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s'%(key,value))
                options[key] = parsed_float

            elif key in ['test_timing']:
                if value is None:
                    options[key] = 'deformation'
                elif value.lower() in ['integrand','deformation']:
                    options[key] = value.lower()
                elif value.lower() in ['f','false']:
                    options[key] = False
                else:
                    raise pyNLoopInvalidCmd("Option 'test_timing' can only take values in ['integrand','deformation','False'], not %s"%value)

            elif key in ['force','log_axis','test','keep_PS_point_fixed','show_plot','integration_space','cpp_integrand']:
                parsed_bool = (value is None) or (value.upper() in ['T','TRUE','ON','Y'])
                options[key] = parsed_bool

            elif key in ['output_folder','save_plot']:
                if os.path.isabs(value):
                    options[key] = value
                else:
                    options[key] = pjoin(MG5DIR,value)

            elif key in ['normalization']:
                if value is None:
                    value = 'distance_real'
                _normalization_modes_supported = ['distance_real',]
                if value not in _normalization_modes_supported:
                    raise pyNLoopInvalidCmd('Normalization modes supported are %s, not %s'%(_normalization_modes_supported, value))
                options['normalization'] = value

            elif key == 'phase_computed':
                if value.upper() in ['REAL', 'R', 'RE']:
                    options[key] = 'Real'
                elif value.upper() in ['IMAGINARY', 'I', 'IM']:
                    options[key] = 'Imaginary'
                elif value.upper() in ['ALL', 'A']:
                    options[key] = 'All'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command plot deformation: %s"%key+
                      "\nAvailable commands are: %s"%(' '.join('--%s'%k for k in options)))

        return new_args, options

    def parse_integrate_loop_options(self, args):
        """ Parsing arguments/options passed to the command integrate_loop."""

        options = {
            'PS_point': 'random',
            'seed': None,
            'sqrt_s': 1000.0,
            'channel': None,
            'cpp_integrand': False,
            'verbosity': 1,
            'nb_CPU_cores': None,
            'compute_analytic_result': True,
            'phase_computed': 'Real',
            'output_folder': pjoin(MG5DIR, 'MyPyNLoop_output'),
            'force': False,
            'loop_momenta_generator': self.parse_lmgc_specification('default'),
        }

        # First combine all value of the options (starting with '--') separated by a space
        opt_args = []
        new_args = []
        for arg in args:
            if arg.startswith('--'):
                opt_args.append(arg)
            elif len(opt_args) > 0:
                opt_args[-1] += ' %s' % arg
            else:
                new_args.append(arg)

        for arg in opt_args:
            try:
                key, value = arg.split('=')
            except:
                key, value = arg, None
            key = key[2:]

            if key == 'PS_point':
                if value not in ['random', ]:
                    raise pyNLoopInvalidCmd('For now, pyNLoops only support the ' +
                                            'specification of a random PS point')
                options[key] = value

            elif key in ['seed', 'batch_size', 'verbosity', 'nb_CPU_cores', 'channel',
                         'survey_n_iterations', 'survey_n_points', 'refine_n_iterations', 'refine_n_points']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s' % (key, value))
                options[key] = parsed_int

            elif key in ['loop_momenta_generator_class','lmgc']:
                options['loop_momenta_generator'] = self.parse_lmgc_specification(eval(value))
            elif key in ['sqrt_s', 'target_accuracy']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s' % (key, value))
                options[key] = parsed_float

            elif key in ['force', 'cpp_integrand']:
                parsed_bool = (value is None) or (value.upper() in ['T', 'TRUE', 'ON', 'Y'])
                options[key] = parsed_bool

            elif key in ['output_folder']:
                if os.path.isabs(value):
                    options[key] = value
                else:
                    options[key] = pjoin(MG5DIR, value)

            elif key == 'phase_computed':
                if value.upper() in ['REAL', 'R', 'RE']:
                    options[key] = 'Real'
                elif value.upper() in ['IMAGINARY', 'I', 'IM']:
                    options[key] = 'Imaginary'
                elif value.upper() in ['ALL', 'A']:
                    options[key] = 'All'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command integrated loop: %s" % key +
                                        "\nAvailable commands are: %s" % (' '.join('--%s' % k for k in options)))

        return new_args, options

    def parse_lmgc_specification(self, specs):
        """ Given a specification of the loop momenta generator class with the format:

                <loop_momenta_generator_class_name>@{'option_a' : value_for_option_a, 'option_b' : value_for_option_b, etc... }

        returns a tuple (loop_momenta_generator_class, options_to_pass_during_instantiation)
        """

        # Hyperparameters of the deformation
        default_loop_momenta_generator_options = {
            'conformal_mapping_choice'  : 'log',
            'M1'                        : 0.035,
            'M2'                        : 0.7,
            'M3'                        : 0.035,
            'M4'                        : 0.035,
            'Gamma1'                    : 0.7,
            'Gamma2'                    : 0.008,
            'Esoft'                     : 0.003,
        }

        if specs == 'default':
            return None, default_loop_momenta_generator_options
        else:
            try:
                lmgc_name, user_lmgc_options = specs.split('@')
            except ValueError:
                lmgc_name    = specs
                user_lmgc_options = '{}'

            try:
                user_lmgc_options = eval(user_lmgc_options)
            except:
                raise pyNLoopInvalidCmd("Loop momenta generator options %s could not be parsed." % user_lmgc_options)
            try:
                lmgc_class = eval('loop_momenta_generator.%s' % lmgc_name)
            except:
                raise pyNLoopInvalidCmd("Loop momentum generator class '%s' not recognized." % lmgc_name)
            lmgc_options = dict(default_loop_momenta_generator_options)
            lmgc_options.update(user_lmgc_options)
            return (lmgc_class, lmgc_options)

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
        for i in random_PS_point:
            if i > 2:
                continue
            random_PS_point[i] = -random_PS_point[i]

        channels = [0,1,2] if options['channel'] is None else [options['channel'],]

         # For debugging you can easily print out the chosen PS point as follows:
#        misc.sprint(str(random_PS_point))
        n_loop_integrand = chosen_topology['class'](
            n_loops            =   chosen_topology['n_loops'],
            external_momenta   =   random_PS_point,
            # Data-structure for specifying a topology to be determined
            topology           =   chosen_topology_name,
            loop_momenta_generator_class = options['loop_momenta_generator'][0],
            loop_momenta_generator_options = options['loop_momenta_generator'][1],
            **options
        )
        
        all_integrands = [n_loop_integrand,]+n_loop_integrand.get_integrated_counterterms()

        all_integrands_amplitude = sympy.sympify('0')
        all_integrands_error = sympy.sympify('0')

        res_summary = None
        for i_integrand, loop_integrand in enumerate(all_integrands):
            
            logger.info("Processing integrand %s ..."%loop_integrand.nice_string())


            # Initialize a hook to OneLOop in order to gather the analytic results
            analytic_results_available = False
            if options['compute_analytic_result']:
                loop_integrand.setup_analytic_computation(self)
                if not loop_integrand.is_analytic_result_available():
                    logger.warning('AVH OneLOop library could not be properly loaded. Analytic results will not be available.')
                    analytic_results_available = False
                else:
                    analytic_results_available = True
            
            # Output low-level integrands code if necessary.
            # For integrands beyond the first, always force.
            loop_integrand.output(
                options['output_folder'], verbosity=options['verbosity'], 
                force = (i_integrand>0 or options['force']) )
    
            if i_integrand==0:
                res_summary = open(pjoin(options['output_folder'],'result.dat'),'w')
                res_summary_lines = ['PS point considered:']
                res_summary_lines.append('%-3s %s'%('#', ' '.join('%-30s'%comp for comp in ['E','p_x','p_y','p_z'])))
                for i_leg, momentum in random_PS_point.items():
                    res_summary_lines.append('%-3d %s'%(i_leg, ' '.join('%-30.16e'%comp for comp in momentum)))
                res_summary_lines.append('')
                res_summary_lines.append(('%-30s'*(n_loop_integrand.n_loops*4+3))%tuple(
                    ['','O(eps^0) Re','O(eps^0) Im']+sum([['O(eps^-%d) Re'%eps_index,'O(eps^-%d) Im'%eps_index] 
                                    for eps_index in range(1,n_loop_integrand.n_loops*2+1)],[])))
                res_summary_lines.append('')
                res_summary.write('\n'.join(res_summary_lines))

            # Choose the integrator
            if self.pyNLoop_options['integrator'] == 'auto':
                integrator_name = loop_integrand._supported_integrators[0]
            else:
                if self.pyNLoop_options['integrator'] not in loop_integrand._supported_integrators:
                    integrator_name = loop_integrand._supported_integrators[0]
                    logger.warning(
        "Specified integrator '%s' is not supported by integrand '%s'. Using '%s' instead."%
         (self.pyNLoop_options['integrator'], loop_integrand.nice_string(), integrator_name))
                else:
                    integrator_name = self.pyNLoop_options['integrator']

            # Now set the integrator options
            integrator_options = dict(self.integrator_options)
            for opt, value in options.items():
                if opt in integrator_options:
                    integrator_options[opt] = value
            integrator_options['cluster'] = self.get_cluster(force_nb_cpu_cores=int(self.options['nb_core']))
            integrator_options['pySecDec_path'] = self.pyNLoop_options['pySecDec_path']
            if integrator_name=='Vegas3':
                if options['phase_computed']=='All':
                    logger.warning('Vegas3 integrator cannot simultaneously integrate the real and imaginary part of'+
                              " the loop for now. The user's choice 'All' for 'phase_computed' is reverted to 'Real'.")
                    loop_integrand.phase_computed = 'Real'
                    options['phase_computed'] = 'Real'
                integrator = Vegas3Integrator(loop_integrand, **integrator_options)
            elif integrator_name=='pySecDec':
                integrator = pysecdec_integrator.pySecDecIntegrator(loop_integrand, **integrator_options)
            elif integrator_name=='Cuba':
                integrator = pyCubaIntegrator(loop_integrand, **integrator_options)
            else:
                raise pyNLoopInterfaceError("Integrator '%s' not implemented."%integrator_name)
            
            # We are now finally ready to integrate :)
            logger.info("="*150)        
            logger.info('{:^150}'.format("Starting integration of %s with integrator %s, lay down and enjoy..."%(
                                                            loop_integrand.nice_string(), integrator.get_name())))
            logger.info("="*150)
    
            amplitude_summed = 0
            error_summed = 0
            for channel in channels:
                # Assign the new kinematic configuration
                loop_integrand.assign_kinematic_configuration(random_PS_point, channel)

                amplitude, error = integrator.integrate()

                if options['verbosity'] > 0:
                    logger.debug('Amplitude for channel {}: {} +/- {}'.format(channel, amplitude, error))

                amplitude_summed += amplitude
                error_summed += error**2


            error_summed = math.sqrt(error_summed)

            # Make sure to cast the result to a sympy expression
            amplitude, error = self.cast_result_to_sympy_expression(amplitude_summed, error_summed)


            # Now compute the analytical result for this point
            analytical_result = None
            if analytic_results_available:
                analytic_res_real, analytic_res_imag = loop_integrand.get_analytic_result(random_PS_point)

                if options['phase_computed'] == 'Real':
                    analytical_result = analytic_res_real
                elif options['phase_computed'] == 'Imaginary':
                    analytical_result = analytic_res_imag

   
            # Aggregate this result to existing ones
            all_integrands_amplitude += amplitude
            all_integrands_error += error
            
            run_output_path = MG5DIR
            self.print_results(loop_integrand, integrator, amplitude, error, options['phase_computed'], analytic=analytical_result, channel=options['channel'])

            # Write the result in 'cross_sections.dat' of the result directory
            self.dump_result_to_file(
                res_summary, loop_integrand.n_loops, amplitude, error, loop_integrand.nice_string())

        logger.info("")
        logger.info("PS point considered in this run:\n%s" % str(random_PS_point))
        logger.info("")
        if len(all_integrands)>1:
            # Now close the result summary file after having reported the aggregated results
            self.print_results(all_integrands[0], integrator,
                    all_integrands_amplitude, all_integrands_error, options['phase_computed'], label='Aggregated results for phase %s' % options['phase_computed'])
        # Write the result in 'cross_sections.dat' of the result directory
        self.dump_result_to_file(res_summary, all_integrands[0].n_loops,
                      all_integrands_amplitude, all_integrands_error, 'Aggregated results')
        res_summary.close()

    def print_results(self, loop_integrand, integrator, amplitude, error, phase, label=None, analytic=None, channel=None,):
        """ Print to screen the results for that particular amplitude evaluation and its MC error."""
        
        if MPI_RANK==0:
            logger.info("="*150)
            if label is None:
                if channel is not None and channel > 0:
                    logger.info('{:^150}'.format("Integral of '%s' (channel %s, %s phase) with integrator '%s':"%(loop_integrand.nice_string(), channel, phase, integrator.get_name())))
                else:
                    logger.info('{:^150}'.format("Integral of '%s' (%s phase) with integrator '%s':"%(loop_integrand.nice_string(), phase, integrator.get_name())))
            else:
                logger.info('{:^150}'.format(label))
            logger.info('')
            logger.info(' '*15+'(%-56s)     +/- (%-60s)'%(amplitude.coeff('eps',0), error.coeff('eps',0)))
            for eps_index in range(1,loop_integrand.n_loops*2+1):
                logger.info(' '*13+'+ (%-56s) %s +/- (%-60s) %s'%(amplitude.coeff('eps',-eps_index), 
                                  EPSILONS[-eps_index], error.coeff('eps',-eps_index), EPSILONS[-eps_index]))
            logger.info('')

            if analytic is not None:
                logger.info('Analytical:    (%-56s)' % analytic[0])
                for eps_index in range(1,loop_integrand.n_loops*2+1):
                    logger.info(' '*15 + '(%-56s) %s' % (analytic[-eps_index], EPSILONS[-eps_index]))
            logger.info("="*150+"\n")

    def dump_result_to_file(self, stream, n_loops, amplitude, error, title):
        
        res_summary_lines = ['contribution_name = %s'%title]
        res_summary_lines.append(('%-30s'*(n_loops*4+3))%tuple(
           ['Result', amplitude.coeff('eps',0).coeff('I',0), amplitude.coeff('eps',0).coeff('I',1) ]+
           sum([[amplitude.coeff('eps',-eps_index).coeff('I',0), amplitude.coeff('eps',-eps_index).coeff('I',1)]
                                                  for eps_index in range(1,n_loops*2+1) ],[])   
        ))
        res_summary_lines.append(('%-30s'*(n_loops*4+3))%tuple(
           ['MC error', error.coeff('eps',0).coeff('I',0), error.coeff('eps',0).coeff('I',1) ]+
           sum([[error.coeff('eps',-eps_index).coeff('I',0), error.coeff('eps',-eps_index).coeff('I',1)]
                                                  for eps_index in range(1,n_loops*2+1) ],[])   
        ))
        stream.write('\n'.join(res_summary_lines))
        
    def cast_result_to_sympy_expression(self, amplitude, error):
        if isinstance(amplitude, float):
            return sympy.sympify('%.16e + 0.*I'%amplitude), sympy.sympify('%.16e + 0.*I'%error)
        elif isinstance(amplitude, tuple):
            return sympy.sympify('%.16e + %.16e*I'%amplitude), sympy.sympify('.16e + %.16e*I'%error)
        elif isinstance(amplitude, list):
            return sympy.sympify('%.16e + %.16e*I'%amplitude[0]+
                                 '(%.16e + %.16e*I)*eps**-1'%amplitude[1]+
                                 '(%.16e + %.16e*I)*eps**-2'%amplitude[2]),\
                   sympy.sympify('%.16e + %.16e*I'%error[0]+
                                 '(%.16e + %.16e*I)*eps**-1'%error[1]+
                                 '(%.16e + %.16e*I)*eps**-2'%error[2])
        else:
            return amplitude, error

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
