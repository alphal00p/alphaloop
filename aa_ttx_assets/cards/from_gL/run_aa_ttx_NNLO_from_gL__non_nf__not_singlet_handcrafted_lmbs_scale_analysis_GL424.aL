import model aL_sm-no_widths

launch /shared/vhirschi/aa_ttx/runs/NNLO_from_gL__non_nf__not_singlet_handcrafted_lmbs_GL424
set_hyperparameter Selectors []
set_hyperparameter CrossSection.incoming_momenta [[250.0,0.,0.,-250.0],[250.0,0.,0.,250.0]]
set_hyperparameter CrossSection.picobarns True
set_hyperparameter General.deformation_strategy "none"
#set_hyperparameter CrossSection.inherit_deformation_for_uv_counterterm True
#set_hyperparameter Deformation.fixed.dampen_on_pinch True
#set_hyperparameter Deformation.fixed.dampen_on_pinch_after_lambda True
#set_hyperparameter Deformation.fixed.pinch_dampening_alpha 1.0
#set_hyperparameter Deformation.scaling.branch_cut_alpha 1.0
#set_hyperparameter Deformation.scaling.lambda 10.0
#set_hyperparameter Deformation.overall_scaling "constant"
#set_hyperparameter Deformation.overall_scaling_constant 10.0
set_hyperparameter CrossSection.m_uv_sq 160000.0
set_hyperparameter CrossSection.mu_r_sq 8315.251344
set_hyperparameter General.multi_channeling True
set_hyperparameter General.use_optimal_channels False
set_hyperparameter General.use_lmb_channels False

#refresh_derived_data

#uv_profile SG_QG0 --min_scaling 1.0e2 --max_scaling 1.0e4

#ir_profile SG_QG0 --ir_limits S(pq7) --f128 --min 1.0e-02 --max 1.0e-04 --n_points 20 --seed 2 --full_report --show_momenta --no_warnings

#integrate SG_QG0 -i inspect -xs 0.6 0.7 0.8 0.9 0.1 0.2

# all \
#--fresh

#--integrand_descriptions m_uv_sq_x_1
#--integrands CrossSection.m_uv_sq=160000.0
integrate \
all \
-i havana -bs 1000 --n_start 10000000 --n_max 1000000000000 --n_increase 0 \
-nw 500 -c 1 --cluster_type local \
--run_id 1447 \
--run_description aa_ttx_NNLO__non_nf__not_singlet_scale_analysis \
--no_mc_over_channels \
--show_selected_phase_only \
--integrand_descriptions mu_r_x_1 mu_r_x_10 mu_r_x_0_1 \
--integrands \
CrossSection.mu_r_sq=8315.251344 \
CrossSection.mu_r_sq=831525.1344 \
CrossSection.mu_r_sq=83.15251344
