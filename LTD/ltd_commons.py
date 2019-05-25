#!/usr/bin/env python2
import os
import sys

from ltd_utils import HyperParameters

#############################################################################################################
# Define and store hyperparameters
#############################################################################################################

hyperparameters = HyperParameters({

    'General'       :   {
        # can be additive, cutgroups, constant, duals or none
        'deformation_strategy'  :   'cutgroups',
        'topology'              :   'Box',
        # specify the name of a python module that contains a function numerator(complex_loop_momenta)
        # that will be called for every cut
        'python_numerator'      :   None,
        # only evaluate the cuts in this list. empty means all
        'cut_filter'            :   [],
        'numerical_threshold'   :   0.,
        # number of digits that should be the same between integrand and rotated version
        'relative_precision'    :   5.,
        # absolute precision, heavily dependent on integral value
        'absolute_precision'    :   1e-10,
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
        'integration_statistics':   True,
        'statistics_interval'   :   100000,
        'debug'                 :   0
    },

    'Integrator'    :   {
        # The integrator can be vegas or cuhre or suave
        'integrator'        :   'vegas',
        'n_start'           :   int(1.0e6),
        'n_max'             :   int(1.0e10),
        'n_increase'        :   int(1.0e6),
        # can be set to high values for use with MPI, otherwise leave it at 1
        'n_vec'             :   1,
        'n_new'             :   100000,
        'n_min'             :   2,
        'flatness'          :   50.,
        'seed'              :   1,
        'integrated_phase'  :  'imag',
        'state_filename'    :   None,
        'survey_n_points'   :   0,
        'survey_n_iterations':  0,
        'refine_n_runs'      :  0,
        'refine_n_points'    :  0,
    },

    'Deformation'   :   {
        # can be constant, linear or sigmoid
        'overall_scaling' : 'linear',
        # fraction of e_cm used for scaling
        'overall_scaling_constant'  : 1.,

        # optionally set a lambda per surface/cutgroup
        'lambdas'   : [],

        'scaling'   :   {
            # positive value: maximum lambda in auto scaling
            # negative value: no auto scaling, lambda is set to abs(lambda)
            'lambda'                    : 1000.0,
            # sigma=0 means normal min. sigma large decreases steepness
            'softmin_sigma'             : 0.0,
            'expansion_check'           : False,
            'expansion_threshold'       : 0.1,
            'positive_cut_check'        : True ,
            'cut_propagator_check'      : True,
            'non_cut_propagator_check'  : True,
            'skip_hyperboloids'         : True,
        },

        'additive'              :   {
            # can be exponential, hyperbolic, or unity
            'mode'  :   'unity',
            'a_ij'  :   0.0000001,
            # set aijs per surface. if the entry isn't there, a_ij is used instead
            'a_ijs' :   [],
        },

        'cutgroups' : {
            'M_ij'  :   0.00001,
            'sigma' :   0.0,
            # can be hyperbolic, softmin, or unity
            'mode'  :   'softmin',
        }
    },

    'Parameterization'   :   {
        # can be cartesian or spherical
        'mode'      :   'spherical',
        # can be log or linear
        'mapping'   :   'linear',
        # controls the UV behaviour of the spherical log map
        'b'         :   1.0,
        # rescale the input from [0,1] to [lo,hi]
        'input_rescaling' : [
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
            [[0., 1.], [0., 1.], [0., 1.]],
        ],
        # shift the loop momenta. the first component is a rescaling of the radius
        'shifts' : [
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.], 
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
