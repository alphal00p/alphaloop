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
set_alphaLoop_option SG_name_list ["SG_QG0","SG_QG6","SG_QG50","SG_QG110","SG_QG122","SG_QG218","SG_QG220","SG_QG240","SG_QG260","SG_QG336","SG_QG350","SG_QG356","SG_QG373","SG_QG790","SG_QG1180","SG_QG1196","SG_QG1216"]
set_FORM_option reference_lmb {\
"SG_QG0":['pq4','pq6','pq7'],\
"SG_QG6":['pq2','pq5','pq7'],\
"SG_QG50":['pq3','pq10','pq7'],\
"SG_QG110":['pq2','pq4','pq10'],\
"SG_QG122":['pq5','pq6','pq9'],\
"SG_QG218":['pq3','pq4','pq6'],\
"SG_QG220":['pq1','pq5','pq7'],\
"SG_QG240":['pq2','pq4','pq9'],\
"SG_QG260":['pq3','pq6','pq9'],\
"SG_QG336":['pq2','pq6','pq8'],\
"SG_QG350":['pq2','pq8','pq10'],\
"SG_QG356":['pq1','pq8','pq10'],\
"SG_QG373":['pq1','pq7','pq10'],\
"SG_QG790":['pq1','pq8','pq10'],\
"SG_QG1180":['pq1','pq7','pq9'],\
"SG_QG1196":['pq2','pq7','pq9'],\
"SG_QG1216":['pq4','pq6','pq9']\
}
set_alphaLoop_option checkpoint_lvl 2
qgraf_generate  a a > t t~ / u d s c b QCD^2==2 QED^2==2 [QCD QCD]
!rm -rf LU_aa_ttx_NNLO_selection
output qgraf LU_aa_ttx_NNLO_selection
