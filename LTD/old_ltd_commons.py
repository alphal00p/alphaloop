#!/usr/bin/env python3
import os
import sys

from ltd_utils import HyperParameters

#############################################################################################################
# Define and store hyperparameters
#############################################################################################################

hyperparameters = HyperParameters({

    'General'       :   {
        # Consider a multi-channeling treatment of the integrand with shifted parametrisations regularising
        # integrable singularities for which one has no deformation.
        'multi_channeling'      :   False,
        # Instead of None, one can specify here a list of indices, like 0,2,7 which corresponds to
        # the channel IDs to consider. A channel ID corresponds to its index in the list produced by the
        # cartesian product of the cut_structure with the propagators in each of the loop lines.
        'multi_channeling_channel': None,
        # can be additive, fixed, constant, duals, intersections or none
        'deformation_strategy'  :   'fixed',
        'topology'              :   'Box',
        'amplitude'             :  '',
        # specify the name of a python module that contains a function numerator(complex_loop_momenta)
        # that will be called for every cut
        'python_numerator'      :   None,
        # counterterms parameters
        'use_ct'                :   False,
        'use_collinear_ct'      :   False,
        'mu_uv_sq_re_im'        :   [1e4,0],
        # only evaluate the cuts in this list. empty means all
        'cut_filter'            :   [],
        'numerical_threshold'   :   0.,
        # number of digits that should be the same between integrand and rotated version
        'relative_precision'    :   3.,
        # absolute precision, heavily dependent on integral value
        'absolute_precision'    :   1e+99,
        'unstable_point_warning_percentage'  :   1.,
        'numerical_instability_check': True,
        'num_digits_different_for_inconsistency': 10.,
        # return the unstable point only if it has more stable digits than specified below
        'minimal_precision_for_returning_result': 2.,
        # which core to log to screen, None logs all cores
        'screen_log_core'       :   1,
        # log max and unstable points to screen
        'log_points_to_screen'  :   False,
        # log statistics to screen
        'log_stats_to_screen'   :   False,
        'log_file_prefix'       :   'stats/statistics',
        'res_file_prefix'       :   '',        
        'integration_statistics':   False,
        'statistics_interval'   :   100000,
        'debug'                 :   0
    },

    'Integrator'    :   {
        # The integrator can be vegas, divonne, cuhre or suave
        'integrator'        :   'cuhre',
        'n_start'           :   int(1.0e5),
        'n_max'             :   int(1.0e10),
        'n_increase'        :   int(1.0e5),
        # can be set to high values for use with MPI, otherwise leave it at 1
        'n_vec'             :   1,
        'seed'              :   1,
        'integrated_phase'  :  'both',
        'state_filename_prefix' :   None,
        'survey_n_points'   :   0,
        'survey_n_iterations':  0,
        'refine_n_runs'      :  0,
        'refine_n_points'    :  0,
        # Vegas related parameters
        'keep_state_file'   : False,
        'reset_vegas_integrator' : True,
        'use_only_last_sample' : False,
        # Non-vegas related integrator parameters
        'eps_rel'           :   1e-8,
        'eps_abs'           :   0.,
        # A border set different to zero allows to not probe particular problematic kinematics
        'border'            :   1.0e-3,
        'n_new'             :   100000,
        'n_min'             :   2,
        'flatness'          :   50.,
        'max_pass'          :   5,
        'maxchisq'          :   0.,
        'mindeviation'      :   0.025,
    },

    'Deformation'   :   {
        # can be constant, linear or sigmoid
        'overall_scaling' : 'constant',
        # fraction of e_cm used for scaling
        'overall_scaling_constant'  : 1.0,
        # A negative number indicates this normalisation is disabled
        # A positive number indicate the value to use in the T function for this normalisation strategy
        'normalize_on_E_surfaces_m' : -1.0,

        # Maximum number of iterations for refinement of the deformation direction by taking shifts
        # of previous deformations into account
        'max_iterations': 1,
        # The stability threshold for when the difference between deformations is considered negligible.
        'stability_threshold': -1,

        # optionally set a different lambda per surface
        'lambdas'   : [],

        'scaling'   :   {
            # positive value: maximum lambda in auto scaling
            # negative value: no auto scaling, lambda is set to abs(lambda)
            'lambda'                    : 1.0,
            # sigma=0 means normal min. sigma large decreases steepness
            'softmin_sigma'             : 0.0,
            'expansion_check'           : True,
            # The expansion check strategy can either be
            # first_lambda_order : c * (q_i^2^cut+m_i^2)/|kappa_i^cut * q_i^cut| 
            # full_lambda_dependence : lambda^2 < (-2*kappa_i.q_i^2 + 
            #           sqrt(4*kappa_i.q_i^4 + kappa^4 c^2 (q_i^2+m_i^2)^2))/kappa^4
            # magic_fudge : c * c * (d / (b + d))
            # magic_fudge_with_min : min( c * c * (d / (b + d)), c * c * (d / (a + d)))
            # ratio: lambda^2 < c * a * a / (a * d - b * b), if c < 1/2 the branch cut check is always satisfied
            'expansion_check_strategy'  : 'ratio',
            # The expansion_threshold considered is min(e_th, 0.5*e_th_max), with e_th_max computed 
            # during the pre-processing. Note instead that a negative value is allowed and the absolute value will
            # then directly be used as an input.
            'expansion_threshold'       : 0.15,
            'positive_cut_check'        : True,
            # The two branchcut M parameters below allow the argument of the square roots
            # to visit all four complex quadrants while still never crossing a branchcut
            'branch_cut_m'              : -1.0,
            'source_branch_cut_m'       : -1.0,
            'cut_propagator_check'      : False,
            'non_cut_propagator_check'  : False,
            'skip_hyperboloids'         : True,
            'source_branch_cut_threshold'  : 0.8,
            'source_branch_cut_multiplier' : 0.8,
            # the strategy can be real_solution (the old method), exact (for 1 loop), and tangent_check
            'pole_check_strategy'       : 'real_solution',
            'theta_r_out'               : 2.,
            'theta_r_in'                : 2.,
            'theta_c'                   : 2.,
        },

        'additive'              :   {
            # can be exponential, hyperbolic, or unity
            'mode'  :   'exponential',
            'a_ij'  :   0.01,
            # set aijs per surface. if the entry isn't there, a_ij is used instead
            'a_ijs' :   [],
        },

        'fixed' : {
            # Maximum allowed value for the anti-selector on an E-surface.
            'delta' : 0.03,
            # Argument of the anti-selector of the form  E_i^2 / ( E_i^2 + M_ij^2 * M^\star ^2 * (p_i^0)^2 )
            # where: E_i is the equation of the E-surface #i
            #        M^\star is the maximum value that M can take for the anti-selector to be < \delta on E-surfaces.
            #        M_ij is the multiplier of M^\star (should typically be >=1.0 for the detla constraint to be satisfied)
            #        p_i^0 is the "surface shift" of E_surface #i.
            'M_ij'  :   1.0,
            # The parameter below is used only for the `softmin` anti-selection.
            'sigma' :   0.0,
            # can be hyperbolic, softmin, or unity
            'mode'  :   'hyperbolic',
            # if not empty, use a different m_ij per numerical surface id (the index in the surface list)
            'm_ijs' : [],
            # localizes the deformation around the ellipsoids using an exponential function with variance a_ij
            'local' : False,
            'a_ijs' : [],
            # Normalize to each source
            'normalize_per_source': False,
            # Normalise the deformation of each subspace with the number of sources in that subspace
            'normalisation_of_subspace_components' : True,
            # Normalise the overall fixed deformation with respect to the number of its building elements (i.e. subspaces)
            'normalisation_per_number_of_sources' : True,
            # add the normal vector for every ellipsoid which is excluded on all other ellipsoids
            'include_normal_source': False,
            # dampen the fixed sources away from intersections. Use this parameter only when include_normal_source is true
            'source_dampening_factor': -0.0005,
        }
    },

    'Parameterization'   :   {
        # can be cartesian or spherical
        'mode'      :   'spherical',
        # can be log or linear,
        # Warning: log is badly behaved in the UV as it does not allow to probe that region enough. Keep linear for safety.
        'mapping'   :   'linear',
        # controls the UV behaviour of the spherical log map
        'b'         :   1.0,
        # rescale the input from [0,1] to [lo,hi]
        'input_rescaling' : [
            [[0. ,1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]]
        ],
        # shift the loop momenta. the first component is a rescaling of the radius
        'shifts' : [
            [1.0, 0., 0., 0.],            
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.]
        ]
    },

})


def synchronize(root_dir = '', sync_topologies=False):
    # Synchronise the database of hard-coded topologies to the yaml data base that Rust can easily import
    # The synchronisation of topologies can take some time, so let's not do that just systematically
    if sync_topologies:
        from topologies import hard_coded_topology_collection
        hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
        print("Synchronised topologies.yaml")
    
    hyperparameters.export_to(os.path.join(root_dir, 'hyperparameters.yaml'))
    print("Synchronised hyperparameters.yaml")

# Main synchronises yaml file to python records
if __name__ == '__main__':
    sync_topologies = any(('sync_topo' in arg.lower() or 'full' in arg.lower()) for arg in sys.argv)
    synchronize(sync_topologies=sync_topologies)
