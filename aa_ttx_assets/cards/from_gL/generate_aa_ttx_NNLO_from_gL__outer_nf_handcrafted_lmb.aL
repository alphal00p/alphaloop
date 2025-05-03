import model aL_sm-no_widths
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model epem
set_FORM_option number_of_lmbs None

set_FORM_option allow_external_momenta_isomorphic_mapping False
set_alphaLoop_option FORM_compile_cores 40
set_FORM_option cores 20
#set_FORM_option max_n_lines_in_C_source 300

set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option checkpoint_lvl 1
set_FORM_option reference_lmb {}
set_qgraf_spoof ./PLUGIN/new_aL_deformation/aa_ttx_assets/graphs/NNLO_from_gL__outer_nf.py
set_FORM_option reference_lmb {\
"GL60":["pq1","pq2","pq9"],\
"GL62":["pq1","pq2","pq9"],\
"GL64":["pq1","pq2","pq6"],\
"GL66":["pq1","pq2","pq6"],\
"GL72":["pq1","pq2","pq9"],\
"GL74":["pq1","pq2","pq9"],\
"GL76":["pq1","pq2","pq10"],\
"GL78":["pq1","pq2","pq10"],\
}
qgraf_generate  a a > t t~ / u d s c b QCD^2==1 QED^2==2 [QCD]
!rm -rf /shared/vhirschi/aa_ttx/runs/NNLO_from_gL__outer_nf_handcrafted_lmb
output qgraf /shared/vhirschi/aa_ttx/runs/NNLO_from_gL__outer_nf_handcrafted_lmb
