#!/usr/bin/env python2
import os

from ltd_utils import HyperParameters 
from topologies import hard_coded_topology_collection

#############################################################################################################
# Define and store hyperparameters
#############################################################################################################

hyperparameters = HyperParameters({

    'General'       :   {
        # can be multiplicative, additive, cutgroups, duals or none
        'deformation_strategy'  :   'cutgroups',
        'topology'              :   'Box',
        # only evaluate the cuts in this list. empty means all
        'cut_filter'            :   [],
        'numerical_threshold'   :   0.,
        # number of digits that should be the same between integrand and rotated version
        'relative_precision'    :   5.,
        # absolute precision, heavily dependent on integral value
        'absolute_precision'    :   1e-5,
        'unstable_point_warning_percentage'  :   1.,
        'numerical_instability_check': True,
        # return the unstable point if true, else return 0
        'return_unstable_point':    False,
        # which core to log to screen, None logs all cores
        'screen_log_core'       :   1,
        # log max and unstable points to screen
        'log_points_to_screen'  :   False,
        # log statistics to screen
        'log_stats_to_screen'   :   True,
        'log_file_prefix'       :   'stats/statistics',
        'integration_statistics':   True,
        'statistics_interval'   :   100000,
        'debug'                 :   0
    },

    'Integrator'    :   {
        # The integrator can be vegas or cuhre or suave
        'integrator'        :   'vegas',
        'n_start'           :   int(1.0e5),
        'n_max'             :   int(1.0e9),
        'n_increase'        :   0,
        'n_new'             :   100000,
        'n_min'             :   2,
        'flatness'          :   50.,
        'seed'              :   1,
        'integrated_phase'  :  'both'
    },

    'Deformation'   :   {
        # can be constant, linear or sigmoid
        'overall_scaling' : 'linear',
        # fraction of e_cm used for scaling
        'overall_scaling_constant': 1.,

        'scaling'   :   {
            # positive value: maximum lambda in auto scaling
            # negative value: no auto scaling, lambda is set to abs(lambda)
            'lambda'                    : 1.0,
            # sigma=0 means normal min. sigma large decreases steepness
            'softmin_sigma'             : 0.0,
            'expansion_check'           : False,
            'expansion_threshold'       : 0.1,
            'positive_cut_check'        : True ,
            'cut_propagator_check'      : True,
            'non_cut_propagator_check'  : True,
            'skip_hyperboloids'         : False,
        },

        'additive'              :   {
            # can be exponential, hyperbolic, or unity
            'mode'  :   'unity',
            'a_ij'  :   0.0000001,
        },

        'multiplicative'        :   {
            'M_ij'  :   0.1
        },

        'cutgroups' : {
            'M_ij'  :   0.00001,
            'sigma' :   0.0,
            # can be hyperbolic, softmin, or unity
            'mode'  :   'unity',
        }
    },

    'Parameterization'   :   {
        # can be cartesian or spherical
        'mode'      :   'spherical',
        # can be log or linear
        'mapping'   :   'log',
        # shift the loop momenta. the first component is a rescaling of the radius
        'shifts' : [
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
            [1.0, 0., 0., 0.],
        ]
    },

})


def synchronize(root_dir = ''):
    # Synchronise the database of hard-coded topologies to the yaml data base that Rust can easily import
    hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
    hyperparameters.export_to(os.path.join(root_dir, 'hyperparameters.yaml'))

# Main synchronises yaml file to python records
if __name__ == '__main__':
    synchronize()
