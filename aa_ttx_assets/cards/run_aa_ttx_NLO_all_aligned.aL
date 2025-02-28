import model aL_sm-no_widths

launch LU_aa_ttx_NLO_all_aligned
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
set_hyperparameter CrossSection.mu_r_sq 160000.0
set_hyperparameter General.multi_channeling True
set_hyperparameter General.use_optimal_channels False
set_hyperparameter General.use_lmb_channels False

refresh_derived_data

#uv_profile SG_QG0 --min_scaling 1.0e2 --max_scaling 1.0e4

#ir_profile SG_QG0 --ir_limits S(pq7) --f128 --min 1.0e-02 --max 1.0e-04 --n_points 20 --seed 2 --full_report --show_momenta --no_warnings

#integrate SG_QG0 -i inspect -xs 0.1 0.2 0.3 0.4 0.5 0.6  
#integrate SG_QG0 -i inspect -xs 0.2 0.3 0.4 0.5 0.6 0.7
#integrate SG_QG0 -i inspect -xs 0.3 0.4 0.5 0.6 0.7 0.8
#integrate SG_QG0 -i inspect -xs 0.4 0.5 0.6 0.7 0.8 0.9
#integrate SG_QG0 -i inspect -xs 0.5 0.6 0.7 0.8 0.9 0.1
#integrate SG_QG0 -i inspect -xs 0.6 0.7 0.8 0.9 0.1 0.2

# all \
#--fresh
integrate \
SG_QG0 SG_QG2 SG_QG4 SG_QG8 SG_QG24 SG_QG26 SG_QG30 SG_QG32 SG_QG38 SG_QG40 SG_QG44 SG_QG46 SG_QG48 SG_QG50 SG_QG52 SG_QG56 SG_QG60 SG_QG64 SG_QG68 SG_QG70 \
--fresh \
-i havana -bs 1000000 --n_start 1000000000 --n_max 1000000000000 --n_increase 0 \
-nw 45 -c 1 --cluster_type local \
--run_id 1337 \
--run_description aa_ttx \
--no_mc_over_channels \
--show_selected_phase_only \
--integrand_descriptions m_uv_sq_x_1 \
--integrands \
CrossSection.m_uv_sq=160000.0
#--integrand_descriptions m_uv_sq_x_1 m_uv_sq_x_0.25 m_uv_sq_x_4 \
#--integrands \
#CrossSection.m_uv_sq=160000.0 \
#CrossSection.m_uv_sq=40000.0 \
#CrossSection.m_uv_sq=640000.0
