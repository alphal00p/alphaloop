#####################################################
#                                                   #
#  Source file of the pyNLoop GGVV MG5aMC plugin.   #
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
import copy

import matplotlib.pyplot as plt

from distutils.version import LooseVersion, StrictVersion

plugin_path = os.path.dirname(os.path.realpath( __file__ ))
import madgraph

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
# It will be nice to experiment with Thomas Hahn's integrators too later
# import madgraph.integrator.pyCubaIntegrator as pyCubaIntegrator

import nloop_integrands
import loop_momenta_generator
import pysecdec_integrator


logger = logging.getLogger('pyNLoop.Interface')

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
            raise pyNLoopInvalidCmd(
                'pyNLoop could not successfully import the pySecDec module.\n'+
                'Make sure it is specified in your $PYTHONPATH.')
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
        super(pyNLoopInterface, self).__init__(*args, **opts)

        # Temporary force pySecDec dependencies to be used
        os.environ['PATH']= os.environ['PATH']+':'+pjoin(self.pyNLoop_options['pySecDec_path'], 'bin')

    def do_set_pyNLoop_option(self, line):
        """ Logic for setting pyNLoop options."""

        args = self.split_arg(line)
        key, value = args[:2]
 
        if key == 'integrator':
            if value not in ['Vegas3','pySecDec','auto']:
                raise pyNLoopInvalidCmd(
"pyNLoop only supports Vegas3 and pySecDec integrator for now (or automatically set with value 'auto').")
        
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
            # Make sure the vector has norm of sqrt_s
            ref_vec = (ref_vec / sum(k**2 for k in ref_vec))*options['sqrt_s']
            offset_vec = offset_vec * options['sqrt_s']
            if not options['test']:
                logger.info('Reference vector used: %s'%str(ref_vec))
                logger.info('Offset_vec vector used: %s'%str(offset_vec))

            # For debugging you can easily print out the options as follows:
            #misc.sprint(options)

            if len(all_n_loop_integrands)==0 or (not options['keep_PS_point_fixed']):
                all_n_loop_integrands = []
                for loop_momenta_generator_class in options['loop_momenta_generator_classes']:
                    all_n_loop_integrands.append( (
                        loop_momenta_generator_class if loop_momenta_generator_class == 'default'
                                                      else loop_momenta_generator_class.__name__ ,
                        chosen_topology['class'](
                            n_loops=chosen_topology['n_loops'],
                            external_momenta=random_PS_point,
                            # Data-structure for specifying a topology to be determined
                            topology=chosen_topology_name,
                            loop_momenta_generator_class = loop_momenta_generator_class,
                            **options
                    ) ))

            ################
            # Analyze poles
            ################

            lmg = all_n_loop_integrands[0][1].loop_momentum_generator
            if (options['test'] or any(item in options['items_to_plot'] for item in ['p','poles'])):
                # compute the positions of the poles on the ref vec line
                poles = []
                imaginary_part_on_poles = []
                for i_prop, q_i in enumerate(lmg.q_is):
                    mass_prop = all_n_loop_integrands[0][1].loop_propagator_masses[i_prop]
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
                        imaginary_part_on_poles = [(ipp/normalization_ipp)**dampening_power for ipp in imaginary_part_on_poles]

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
        x_entries = [scaling_range[0]+ (i / float(n_points))*(scaling_range[1]-scaling_range[0]) for i in range(1, n_points)]

        # now sample n_points points and compute the deformation
        points = [ offset_vec + ref_vec * map_to_infinity(x) for x in x_entries]
        normalizations = [1.,]*n_points
        if options['normalization']=='distance_real':
            normalizations = [ math.sqrt(sum(k_i**2 for k_i in p)) for p in points ]

        # Make sure the normalization is capped to be minimum 1.0
        normalizations = [max(n,1.0e-99) for n in normalizations]

        for lm_generator_name, n_loop_integrand in all_n_loop_integrands:

            deformed_points = []
            jacobians       = []
            widgets = ["Loop deformation for generator %s :"%lm_generator_name,
                pbar.Percentage(), ' ', pbar.Bar(),' ', pbar.ETA(), ' ',
                pbar.Counter(format='%(value)d/%(max_value)d'), ' ']
            progress_bar = pbar.ProgressBar(widgets=widgets, maxval=len(points), fd=sys.stdout)
#           progress_bar = None
            if progress_bar:
                progress_bar.start()
            for i_point, p in enumerate(points):
                res = n_loop_integrand.loop_momentum_generator.apply_deformation([p, ])
                deformed_points.append(res[0][0])
                jacobians.append(res[1])
                if progress_bar:
                    progress_bar.update(i_point)
            if progress_bar:
                progress_bar.finish()

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

            # Also add the integrand
            if any(item in options['items_to_plot'] for item in ['integrand_real','integrand_imag','integrand_abs']):
                integrands = (all_n_loop_integrands[0][1]([dp,], [], input_already_in_infinite_hyperbox=True, jacobian=jac )
                                                                           for dp,jac in zip(deformed_points,jacobians))
                # Normalize to the largest value of the integrand encountered
                if 'integrand_real' in options['items_to_plot']:
                    integrand_reals = [itg.real for itg in integrands]
                    max_integrand_real = max(integrand_reals)
                    integrand_reals = [itg/max_integrand_real for itg in integrand_reals]
                    plt.plot(x_entries,integrand_reals, label='integrand_real @%s'%lm_generator_name)
                if 'integrand_imag' in options['items_to_plot']:
                    integrand_imags = [itg.imag for itg in integrands]
                    max_integrand_imag = max(integrand_imags)
                    integrand_imags = [itg/max_integrand_imag for itg in integrand_imags]
                    plt.plot(x_entries,integrand_imags, label='integrand_imag @%s'%lm_generator_name)
                if 'integrand_abs' in options['items_to_plot']:
                    integrand_abss = [abs(itg) for itg in integrands]
                    max_integrand_abss = max(integrand_abss)
                    integrand_abss = [itg/max_integrand_abss for itg in integrand_abss]
                    plt.plot(x_entries,integrand_abss, label='integrand_abs @%s'%lm_generator_name)
            if any(item in options['items_to_plot'] for item in [
                                               'integrand_no_jac_real','integrand_no_jac_imag','integrand_no_jac_abs']):
                integrands = (all_n_loop_integrands[0][1]([dp,], [], input_already_in_infinite_hyperbox=True, jacobian=1.0 )
                                                                           for dp,jac in zip(deformed_points,jacobians))
                # Normalize to the largest value of the integrand encountered
                if 'integrand_no_jac_real' in options['items_to_plot']:
                    integrand_reals = [itg.real for itg in integrands]
                    max_integrand_real = max(integrand_reals)
                    integrand_reals = [itg/max_integrand_real for itg in integrand_reals]
                    plt.plot(x_entries,integrand_reals, label='integrand_no_jac_real @%s'%lm_generator_name)
                if 'integrand_no_jac_imag' in options['items_to_plot']:
                    integrand_imags = [itg.imag for itg in integrands]
                    max_integrand_imag = max(integrand_imags)
                    integrand_imags = [itg/max_integrand_imag for itg in integrand_imags]
                    plt.plot(x_entries,integrand_imags, label='integrand_no_jac_imag @%s'%lm_generator_name)
                if 'integrand_no_jac_abs' in options['items_to_plot']:
                    integrand_abss = [abs(itg) for itg in integrands]
                    max_integrand_abss = max(integrand_abss)
                    integrand_abss = [itg/max_integrand_abss for itg in integrand_abss]
                    plt.plot(x_entries,integrand_abss, label='integrand_no_jac_abs @%s'%lm_generator_name)
            if any(item in options['items_to_plot'] for item in ['jac_real','jac_imag','jac_abs']):
                # Normalize to the largest value of the integrand encountered
                if 'jac_real' in options['items_to_plot']:
                    jacobian_reals = [jac.real for jac in jacobians]
                    max_jacobian_real = max(jacobian_reals)
                    jacobian_reals = [jac/max_jacobian_real for jac in jacobian_reals]
                    plt.plot(x_entries,jacobian_reals, label='jac_real @%s'%lm_generator_name)
                if 'jac_imag' in options['items_to_plot']:
                    jacobian_imags = [jac.imag for jac in jacobians]
                    max_jacobian_imag = max(jacobian_imags)
                    jacobian_imags = [jac/max_jacobian_imag for jac in jacobian_imags]
                    plt.plot(x_entries,jacobian_imags, label='jac_imag @%s'%lm_generator_name)
                if 'jac_abs' in options['items_to_plot']:
                    jacobian_abss = [abs(jac) for jac in jacobians]
                    max_jacobian_abss = max(jacobian_abss)
                    jacobian_abss = [jac/max_jacobian_abss for jac in jacobian_abss]
                    plt.plot(x_entries,jacobian_abss, label='jac_abs @%s'%lm_generator_name)



        if options['log_axis']:
            plt.semilogy()
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3)
        plt.show()

    def parse_integrate_loop_options(self, args):
        """ Parsing arguments/options passed to the command integrate_loop."""

        options = {
            'PS_point': 'random',
            'seed': None,
            'sqrt_s': 1000.0,
            'target_accuracy': 1.0e-3,
            'batch_size': 1000,
            'verbosity': 1,
            'nb_CPU_cores': None,
            'phase_computed': 'All',
            'survey_n_iterations': 10,
            'survey_n_points': 2000,
            'refine_n_iterations': 10,
            'refine_n_points': 1000,
            'output_folder': pjoin(MG5DIR, 'MyPyNLoop_output'),
            'force': False,
            'loop_momenta_generator_class': None,
            # Hyperparameters of the deformation
            'conformal_mapping_choice' : 'log',
            'M1_factor'     : 0.035,
            'M2_factor'     : 0.7,
            'M3_factor'     : 0.035,
            'gamma1'        : 0.7,
            'gamma2'        : 0.008,
            'Esoft_factor'  : 0.003,
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

            elif key in ['M1_factor', 'M2_factor', 'M3_factor', 'gamma1', 'gamma2', 'Esoft_factor']:
                try:
                    options[key] = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot set deformation parameter %s to %s.'%(key, value))

            elif key=='conformal_mapping_choice':
                _valid_conformal_mapping_choices=['lin','log']
                if value.lower() not in _valid_conformal_mapping_choices:
                    raise pyNLoopInvalidCmd('Conformal mapping choice can only be one of: %s'%_valid_conformal_mapping_choices)
                options[key] = value.lower()

            elif key in ['seed', 'batch_size', 'verbosity', 'nb_CPU_cores',
                         'survey_n_iterations', 'survey_n_points', 'refine_n_iterations', 'refine_n_points']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s' % (key, value))
                options[key] = parsed_int

            elif key in ['loop_momenta_generator_class','lmgc']:
                try:
                    lmgc = eval('loop_momenta_generator.%s'%value)
                except:
                    raise pyNLoopInvalidCmd("Loop momentum generator class '%s' not reckognized."%value)
                options['loop_momenta_generator_class'] = lmgc

            elif key in ['sqrt_s', 'target_accuracy']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s' % (key, value))
                options[key] = parsed_float

            elif key in ['force']:
                parsed_bool = (value is None) or (value.upper() in ['T', 'TRUE', 'ON', 'Y'])
                options[key] = parsed_bool

            elif key in ['output_folder']:
                if os.path.isabs(value):
                    options[key] = value
                else:
                    options[key] = pjoin(MG5DIR, value)

            elif key == 'phase_computed':
                if value.upper() in ['REAL', 'R', 'RE']:
                    options[key] == 'Real'
                elif value.upper() in ['IMAGINARY', 'I', 'IM']:
                    options[key] == 'Imaginary'
                elif value.upper() in ['ALL', 'A']:
                    options[key] == 'All'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command integrated loop: %s" % key +
                                        "\nAvailable commands are: %s" % (' '.join('--%s' % k for k in options)))

        return new_args, options

    def parse_plot_deformation_options(self, args):
        """ Parsing arguments/options passed to the command integrate_loop."""

        options = { 
            'seed'                  : None,
            'sqrt_s'                : 1000.0,
            'verbosity'             : 1,
            'output_folder'         : pjoin(MG5DIR,'MyPyNLoop_output'),
            'force'                 : False,
            'phase_computed'        : 'All',
            'reference_vector'      : 'random',
            'offset_vector'         : 'random',
            'loop_momenta_generator_classes' : [None,],
            'n_points'              : 100,
            'items_to_plot'         : (0,1,2,3,'distance','poles'),
            'range'                 : (0.,1.),
            'normalization'         : 'None',
            'log_axis'              : False,
            # When 'test' is on, a battery of tests is performed instead of the plotting of the deformation
            'test'                  : False,
            'keep_PS_point_fixed'   : False,
            # Hyperparameters of the deformation
            'conformal_mapping_choice': 'log',
            'M1_factor'             : 0.035,
            'M2_factor'             : 0.7,
            'M3_factor'             : 0.035,
            'gamma1'                : 0.7,
            'gamma2'                : 0.008,
            'Esoft_factor'          : 0.003,
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

            if key in ['seed','verbosity','nb_CPU_cores', 'n_points']:
                try:
                    parsed_int = int(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s integer: %s'%(key,value))
                options[key] = parsed_int

            elif key in ['M1_factor', 'M2_factor', 'M3_factor', 'gamma1', 'gamma2', 'Esoft_factor']:
                try:
                    options[key] = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot set deformation parameter %s to %s.'%(key, value))

            elif key=='conformal_mapping_choice':
                _valid_conformal_mapping_choices=['lin','log']
                if value.lower() not in _valid_conformal_mapping_choices:
                    raise pyNLoopInvalidCmd('Conformal mapping choice can only be one of: %s'%_valid_conformal_mapping_choices)
                options[key] = value.lower()

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
                    raise pyNLoopInvalidCmd("Cannot parse reference vector specification: %s"%str(value))
                lmgc = []
                for class_name in lmgc_names:
                    if class_name=='default':
                        lmgc.append(None)
                    else:
                        try:
                            lmgc.append(eval('loop_momenta_generator.%s'%class_name))
                        except:
                            raise pyNLoopInvalidCmd("Loop momentum generator class '%s' not reckognized."%class_name)
                options['loop_momenta_generator_classes'] = lmgc

            elif key in ['sqrt_s']:
                try:
                    parsed_float = float(value)
                except ValueError:
                    raise pyNLoopInvalidCmd('Cannot parse specified %s float: %s'%(key,value))
                options[key] = parsed_float

            elif key in ['force','log_axis','test','keep_PS_point_fixed']:
                parsed_bool = (value is None) or (value.upper() in ['T','TRUE','ON','Y'])
                options[key] = parsed_bool

            elif key in ['output_folder']:
                if os.path.isabs(value):
                    options[key] = value
                else:
                    options[key] = pjoin(MG5DIR,value)

            elif key in ['normalization']:
                _normalization_modes_supported = ['distance_real',]
                if value not in _normalization_modes_supported:
                    raise pyNLoopInvalidCmd('Normalization modes supported are %s, not %s'%(_normalization_modes_supported, value))
                options['normalization'] = value

            elif key == 'phase_computed':
                if value.upper() in ['REAL', 'R', 'RE']:
                    options[key] == 'Real'
                elif value.upper() in ['IMAGINARY', 'I', 'IM']:
                    options[key] == 'Imaginary'
                elif value.upper() in ['ALL', 'A']:
                    options[key] == 'All'
                else:
                    raise pyNLoopInvalidCmd("The phase computed can only be 'Real' or 'Imaginary'.")

            else:
                raise pyNLoopInvalidCmd("Unrecognized option for command plot deformation: %s"%key+
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
        for i in random_PS_point:
            if i > 2:
                continue
            random_PS_point[i] = -random_PS_point[i]

        # For debugging you can easily print out the chosen PS point as follows:
#        misc.sprint(str(random_PS_point))
        n_loop_integrand = chosen_topology['class'](
            n_loops            =   chosen_topology['n_loops'],
            external_momenta   =   random_PS_point,
            # Data-structure for specifying a topology to be determined
            topology           =   chosen_topology_name,
            **options
        )
        
        all_integrands = [n_loop_integrand,]+n_loop_integrand.get_integrated_counterterms()

        all_integrands_amplitude = sympy.sympify('0')
        all_integrands_error = sympy.sympify('0')

        res_summary = None
        for i_integrand, loop_integrand in enumerate(all_integrands):
            
            logger.info("Processing integrand %s ..."%loop_integrand.nice_string())
            
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
                 
            # Generic option for all integrators
            cluster = self.get_cluster(force_nb_cpu_cores=options['nb_CPU_cores'])
            integration_options = {
                'cluster'               :   cluster,
                'target_accuracy'       :   options['target_accuracy'],
                'batch_size'            :   options['batch_size'],
                'verbosity'             :   options['verbosity'],
                'survey_n_iterations'   :   options['survey_n_iterations'],
                'survey_n_points'       :   options['survey_n_points'],
                'refine_n_iterations'   :   options['refine_n_iterations'],
                'refine_n_points'       :   options['refine_n_points'],
                'pySecDec_path'         :   self.pyNLoop_options['pySecDec_path']
            }
                    
            if integrator_name=='Vegas3':
                if options['phase_computed']=='All':
                    logger.warning('Vegas3 integrator cannot simultaneously integrate the real and imaginary part of'+
                              " the loop for now. The user's choice 'All' for 'phase_computed' is reverted to 'Real'.")
                    loop_integrand.phase_computed = 'Real'
                integrator = Vegas3Integrator(loop_integrand, **integration_options)
            elif integrator_name=='pySecDec':
                integrator = pysecdec_integrator.pySecDecIntegrator(loop_integrand, **integration_options)
            else:
                raise pyNLoopInterfaceError("Integrator '%s' not implemented."%integrator_name)
            
            # We are now finally ready to integrate :)
            logger.info("="*150)        
            logger.info('{:^150}'.format("Starting integration of %s with integrator %s, lay down and enjoy..."%(
                                                            loop_integrand.nice_string(), integrator.get_name())))
            logger.info("="*150)
    
            amplitude, error = integrator.integrate()
            # Make sure to cast the result to a sympy expression
            amplitude, error = self.cast_result_to_sympy_expression(amplitude, error)
    
            # Aggregate this result to existing ones
            all_integrands_amplitude += amplitude
            all_integrands_error += error
            
            run_output_path = MG5DIR
            self.print_results(loop_integrand, integrator, amplitude, error)
            
            # Write the result in 'cross_sections.dat' of the result directory
            self.dump_result_to_file(
                res_summary, loop_integrand.n_loops, amplitude, error, loop_integrand.nice_string())

        logger.info("")
        logger.info("PS point considered in this run:\n%s" % str(random_PS_point))
        logger.info("")
        if len(all_integrands)>1:
            # Now close the result summary file after having reported the aggregated results
            self.print_results(all_integrands[0], integrator,
                    all_integrands_amplitude, all_integrands_error, label='Aggregated results')
        # Write the result in 'cross_sections.dat' of the result directory
        self.dump_result_to_file(res_summary, all_integrands[0].n_loops,
                      all_integrands_amplitude, all_integrands_error, 'Aggregated results')
        res_summary.close()

    def print_results(self, loop_integrand, integrator, amplitude, error, label=None):
        """ Print to screen the results for that particular amplitude evaluation and its MC error."""
        
        if MPI_RANK==0:
            logger.info("="*150)
            if label is None:
                logger.info('{:^150}'.format("Integral of '%s' with integrator '%s':"%(loop_integrand.nice_string(), integrator.get_name())))
            else:
                logger.info('{:^150}'.format(label))
            logger.info('')
            logger.info(' '*15+'(%-56s)     +/- (%-60s)'%(amplitude.coeff('eps',0), error.coeff('eps',0)))
            for eps_index in range(1,loop_integrand.n_loops*2+1):
                logger.info(' '*13+'+ (%-56s) %s +/- (%-60s) %s'%(amplitude.coeff('eps',-eps_index), 
                                  EPSILONS[-eps_index], error.coeff('eps',-eps_index), EPSILONS[-eps_index]))
            logger.info('')
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
