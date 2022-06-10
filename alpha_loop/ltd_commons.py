#!/usr/bin/env python3
import os
import sys
import logging

logger = logging.getLogger('alphaLoop.ltd_commons')

#############################################################################################################
# HyperParameters
#############################################################################################################
class HyperParameters(dict):

    def __init__(self, *args, **opts):
        super(HyperParameters, self).__init__(*args, **opts)

    def to_flat_format(self):
        """ Turn this instance into a flat dictionary made out of simple lists or dictionaries only."""
        
        # For now the hyperparameters dict is already supposed to be directly exportable to yaml
        return dict(self)

    @staticmethod
    def from_flat_format(flat_dict):
        """ Creates an instance of this class from a flat dictionary record."""
        
        # Directly get HyperParameters from the flat_dict
        return HyperParameters(flat_dict)

    def export_to(self, output_path, format='yaml'):
        """ Exports these hyperparameters to a given format."""
        
        export_format = format.lower()
        allowed_export_format = ['yaml']
        if export_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be exported in the following formats: %s"%(', '.join(allowed_export_format)))

        if export_format=='yaml':
            try:
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to export hyperparameters to yaml.")

        if output_path is not None:
            with open(output_path,'w') as f:
                f.write(yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False))
        else:
            return yaml.dump(self.to_flat_format(), Dumper=Dumper, default_flow_style=False)


    def update_folder(self, folder, dict_update, allow_undefined=False):

        for k, v in dict_update.items():
            if k not in folder:
                if not allow_undefined:
                    logger.exception("Could not find specified directory '%s'."%k)
                else:
                    folder[k] = v
            else:
                if not isinstance(v, dict):
                    folder[k] = v
                else:
                    self.update_folder(folder[k],v,allow_undefined=allow_undefined)

    def update(self, dict_update, allow_undefined=False):
        """Give the possibility of overwriting parameters with a dictionary."""

        self.update_folder(self, dict_update, allow_undefined=allow_undefined)

    def set_parameter(self, param, value, allow_undefined=False):
        """Give the possibility of overwriting options by specifyin param path with format A.B.C..."""
        
        directories = param.split('.')
        folder = self
        for curr_dir in directories[:-1]:
            try:
                folder = folder[curr_dir]
            except KeyError:
                logger.exception("Could not find specified option directory '%s'."%curr_dir)
                return False

        # Set value specified
        if (directories[-1] not in folder) and not allow_undefined:
            logger.exception("Could not find specified option '%s' in folder '%s'."%(
                directories[-1], '.'.join(directories[:-1])
            ))
            return False
        else:
            folder[directories[-1]] = value

        return True

    @staticmethod
    def import_from(input_path, format='yaml'):
        """ Imports this topology from a given format."""
        
        import_format = format.lower()
        allowed_import_format = ['yaml']
        if import_format not in ['yaml']:
            raise BaseException("Hyperparameters can only be imported from the following formats: %s"%(', '.join(allowed_import_format)))

        if import_format=='yaml':
            try: 
                import yaml
                from yaml import Loader, Dumper
            except ImportError:
                raise BaseException("Install yaml python module in order to import hyperparameters from yaml.")
 
        if '\n' in input_path:
            return HyperParameters.from_flat_format(yaml.load(input_path, Loader=Loader))
        else:
            return HyperParameters.from_flat_format(yaml.load(open(input_path,'r'), Loader=Loader))

#############################################################################################################
# Define and store hyperparameters
#############################################################################################################

hyperparameters = HyperParameters({

    'General'       :   {
        # Consider a multi-channeling treatment of the integrand with shifted parametrisations regularising
        # integrable singularities for which one has no deformation.
        'multi_channeling'      :   True,
        # Instead of None, one can specify here a list of indices, like 0,2,7 which corresponds to
        # the channel IDs to consider. A channel ID corresponds to its index in the list produced by the
        # cartesian product of the cut_structure with the propagators in each of the loop lines.
        'multi_channeling_channel': None,
        'multi_channeling_alpha': -1.0,
        'use_optimal_channels': True,
        'use_lmb_channels' : False,
        # Derive the overlap structure required for the fixed deformation in Rust.
        'derive_overlap_structure':  False,
        # can be additive, fixed, constant or none
        'deformation_strategy'  :   'none',
        'topology'              :   'Box',
        # scale at which to use partial fractioning for the integrand
        # and choose whether or not to use the latest partial fractioning
        # implementation
        'partial_fractioning_threshold' :  -1,
        'partial_fractioning_multiloop' :  True,
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
        # absolute precision, heavily dependent on integral value
        'absolute_precision'    :   1e+99,
        # randomly shift each component in x-space by plus or minus stability_nudge_size
        'stability_nudge_size'  :   0.,
        # return the unstable point only if it has more stable digits than specified below
        'minimal_precision_for_returning_result': 999.,
        # the stability pipeline
        'stability_checks'      : [
            {
                # number of samples to take for the numerical stability check
                'n_samples': 3,
                'prec': 16,
                'use_pf': False,
                # number of digits that should be the same between rotated versions
                'relative_precision': 10.0,
                # force an upgrade when a new weight is this threshold times the current maximum weight
                'escalate_for_large_weight_threshold': 0.8,
                'minimal_precision_to_skip_further_checks': 999.0,
                'accepted_radius_range_in_x_space': [0, 0.9]
            },
            {            
                'n_samples': 3,
                'prec': 16,
                'use_pf': True,
                'relative_precision': 5.0,
                'escalate_for_large_weight_threshold': 0.8,
                'minimal_precision_to_skip_further_checks': 999.0,
                'accepted_radius_range_in_x_space': [0, 1]
            },
            {
                'n_samples': 3,
                'prec': 32,
                'use_pf': True,
                'relative_precision': 8.0,
                'escalate_for_large_weight_threshold': -1.,
                'minimal_precision_to_skip_further_checks': 999.0,
                'accepted_radius_range_in_x_space': [0, 1]
            }
        ],
        'res_file_prefix'       :   '',
        'debug'                 :   0
    },

    'Integrator'    :   {
        # Rust will do the parallelization instead of Cuba. Make sure to set n_vec to about 50 to 100 times
        # the numer of cores for best performance.
        'internal_parallelization'  : True,
        # Use the dashboard. Make sure to run Cuba with 0 cores or enable internal_parallelization.
        'dashboard'         :   True,
        # Do not output integrand statistics when the dashboard is disabled.
        'quiet_mode'        :   False,
        'show_plot'         :   False,
        # The integrator can be havana, vegas, divonne, cuhre or suave
        'integrator'        :   'havana',
        'n_start'           :   int(1.0e5),
        'n_max'             :   int(1.0e10),
        'n_increase'        :   int(1.0e5),
        'max_discrete_bin_probability_ratio'    : 100.,
        'min_samples_for_update'    : 1000,
        'n_bins'            : 16,
        'train_on_avg'      : False,
        'learning_rate'     : 1.5,
        # can be set to high values for use with MPI or internal_parallelization, otherwise leave it at 1
        'n_vec'             :   100000,
        'seed'              :   1,
        'integrated_phase'  :  'real',
        'state_filename_prefix' :   None,
        'survey_n_points'   :   0,
        'survey_n_iterations':  0,
        'refine_n_runs'      :  0,
        'refine_n_points'    :  0,
        # Vegas related parameters
        'load_from_state_file': False,
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
        'maxpass'          :   5,
        'maxchisq'          :   0.,
        'mindeviation'      :   0.025,
    },

    'Deformation'   :   {
        # can be constant, linear, sigmoid, or exp_dampening
        'overall_scaling' : 'constant',
        # fraction of e_cm used for scaling
        'overall_scaling_constant'  : 10.0,
        # A negative number indicates this normalisation is disabled
        # A positive number indicate the value to use in the T function for this normalisation strategy
        'normalize_on_E_surfaces_m' : -1.0,

        # optionally set a different lambda per surface
        'lambdas'   : [],

        'scaling'   :   {
            # positive value: maximum lambda in auto scaling
            # negative value: no auto scaling, lambda is set to abs(lambda)
            'lambda'                    : 10.0,
            # sigma=0 means normal min. sigma large decreases steepness
            'softmin_sigma'             : 0.0,
            # The expansion check strategy can either be
            # first_lambda_order : c * (q_i^2^cut+m_i^2)/|kappa_i^cut * q_i^cut| 
            # full_lambda_dependence : lambda^2 < (-2*kappa_i.q_i^2 + 
            #           sqrt(4*kappa_i.q_i^4 + kappa^4 c^2 (q_i^2+m_i^2)^2))/kappa^4
            # magic_fudge : c * c * (d / (b + d))
            # magic_fudge_with_min : min( c * c * (d / (b + d)), c * c * (d / (a + d)))
            # ratio: lambda^2 < c * a * a / (a * d - b * b), if c < 1/2 the branch cut check is always satisfied
            # none
            'expansion_check_strategy'  : 'ratio',
            # The expansion_threshold considered is min(e_th, 0.5*e_th_max), with e_th_max computed 
            # during the pre-processing. Note instead that a negative value is allowed and the absolute value will
            # then directly be used as an input.
            'expansion_threshold'       : -0.3,
            'branch_cut_check'          : True,
            # take the branch cut lambda to a higher power to dampen faster around focal points
            'branch_cut_alpha'          : 1.0,
            # The two branchcut M parameters below allow the argument of the square roots
            # to visit all four complex quadrants while still never crossing a branchcut
            'branch_cut_m'              : -1.0,
            'source_branch_cut_m'       : -1.0,
            'source_branch_cut_threshold'  : 0.8,
            'source_branch_cut_multiplier' : 0.8,
            # the strategy can be none, real_solution (the old method), exact (for 1 loop), and tangent_check
            'pole_check_strategy'       : 'real_solution',
            'soft_dampening_power'      : 1.5,
            'theta_r_out'               : 10.,
            'theta_r_in'                : 10.,
            'theta_c'                   : 10.,
        },

        'fixed' : {
            # Maximum allowed value for the anti-selector on an E-surface.
            'delta' : 0.3,
            # Argument of the anti-selector of the form  E_i^2 / ( E_i^2 + M_ij^2 * M^\star ^2 * (p_i^0)^2 )
            # where: E_i is the equation of the E-surface #i
            #        M^\star is the maximum value that M can take for the anti-selector to be < \delta on E-surfaces.
            #        M_ij is the multiplier of M^\star (should typically be >=1.0 for the detla constraint to be satisfied)
            #        p_i^0 is the "surface shift" of E_surface #i.
            # If M_ij or m_ijs are set negative then the corresponding absolute value is directly used as an input,
            # and not as a rescaling of the m parameter deduced from the value of delta.
            'M_ij'  :   -0.069,
            # The parameter below is used only for the `softmin` anti-selection.
            'sigma' :   0.0,
            # can be hyperbolic, softmin, or unity
            'mode'  :   'hyperbolic',
            # dampen the deformation on pinches
            'dampen_on_pinch': True,
            # dampen the deformation on pinches (if dampen_on_pinch is set) after the lambda scaling
            'dampen_on_pinch_after_lambda': True,
            'pinch_dampening_alpha': 2.0,
            'pinch_dampening_k_com': 1.,
            'pinch_dampening_k_shift': 0.,
            # Specify a strategy to veto small problematic regions close to IR singularities
            # The strategy can either be "dismiss_point", "dismiss_deformation" or "none"
            'ir_handling_strategy': 'none',
            'ir_alpha': 1.5,
            'ir_k_com': 0.,
            'ir_k_shift': 1.0,
            'ir_beta_energy': 1.0,
            'ir_beta_ellipse': 1.0,
            'ir_beta_pinch': 1.0,
            'ir_interpolation_length': 0.0,
            'ir_threshold': 1.0e-8,
            # if not empty, use a different m_ij per numerical surface id (the index in the surface list)
            'm_ijs' : [],
            # localizes the deformation around the ellipsoids using an exponential function with variance a_ij
            'local' : False,
            'a_ijs' : [],
            # use heuristics for the center finding instead of ECOS. Heuristics are built for one ellipsoid
            # and it will test when 0 is on the inside.
            'use_heuristic_centers': True,
            # Introduce a radius variable into the SOCP problem and maximize it
            'maximize_radius': False,
            # Normalize to each source
            'normalize_per_source': False,
            # Normalise the deformation of each subspace with the number of sources in that subspace
            'normalisation_of_subspace_components' : True,
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
            [[0., 1.], [0., 1.], [0., 1.]],
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

    'CrossSection'   :   {
        'incoming_momenta'                      :   [[500., 0., 0., 500.],[500., 0., 0., -500.]],
        # used to compute an amplitude instead of a cross section
        'fixed_cut_momenta'                     :   [],
        'm_uv_sq'                               :   155.**2,
        'mu_r_sq'                               :   155.**2,
        'gs'                                    :   1.2177157847767195,
        # give massless propagators a small mass (only for FORM_integrand)
        'small_mass_sq'                         :   0.,
        'uv_cutoff_scale_sq'                    :   500.*500.,
        'picobarns'                             :   False,
        'inherit_deformation_for_uv_counterterm':   False,
        'do_rescaling'                          :   True,        
        'NormalisingFunction' : {
            # Two possible normalising functions for now: 'left_right_exponential', 'left_right_polynomial' or 'right_exponential'
            # The former dampens both t=0 and t=infty, while the second only dampens t=infty
            'name'                              :   'left_right_exponential',
            'center'                            :   1.0,
            'spread'                            :   1.0,
        },
        # Can be LTD, PF
        'integrand_type'                        :   'PF',
        # evaluate the C expression for the sum of diagrams
        'sum_diagram_sets'                      :   False,
        # compare locally against the same topology written in another loop momentum basis
        'compare_with_additional_topologies'    :   False,
    },

    'Observables'   :  [
        # options: jet1pt, AFB, cross_section
        {
            'type': 'cross_section'
        },
#        {
#            'type'                  :  'jet1pt',
#            'x_min'                 :   0.,
#            'x_max'                 :   0.8,
#            'n_bins'                :   50,
#            'dR'                    :   0.4,
#            'min_jpt'               :   0.,
#            'use_fastjet'           :   True,
#            'write_to_file'         :   True,
#            'filename'              :   'Jet1PT.HwU'
#        }
    ],

    'Selectors'   :   [
        # options: jet, ranged
#        {
#            'type'                  :   'ranged',
#            'pdgs'                  :   [-6,6,25],
#            # options: E, CosThetaP, pT
#            'filter'                :   'E',
#            'min_value'             :   0.,
#            'max_value'             :   1000.,
#        },
#        {
#            'type'                  :   'jet',
#            'min_jets'              :   2,
#            'max_jets'              :   100,
#            'min_j1pt'              :   0.,
#            # A negative maximum means no cut is applied
#            'max_j1pt'              :   -1.,
#            'dR'                    :   0.4,
#            'min_jpt'               :   10.0,
#            'use_fastjet'           :   True,
#        }
    ],

})


def synchronize(root_dir = '', sync_topologies=False):
    # Synchronise the database of hard-coded topologies to the yaml data base that Rust can easily import
    # The synchronisation of topologies can take some time, so let's not do that just systematically
    if sync_topologies:
        from LTD.topologies import hard_coded_topology_collection
        hard_coded_topology_collection.export_to(os.path.join(root_dir, 'topologies.yaml'))
        print("Synchronised topologies.yaml")
    
    hyperparameters.export_to(os.path.join(root_dir, 'hyperparameters.yaml'))
    print("Synchronised hyperparameters.yaml")

# Main synchronises yaml file to python records
if __name__ == '__main__':
    sync_topologies = any(('sync_topo' in arg.lower() or 'full' in arg.lower()) for arg in sys.argv)
    synchronize(sync_topologies=sync_topologies)
