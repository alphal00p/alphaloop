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
set_alphaLoop_option FORM_compile_cores 1
set_FORM_option cores 1
#set_FORM_option max_n_lines_in_C_source 300

set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option SG_name_list ["SG_QG0","SG_QG2","SG_QG4","SG_QG8","SG_QG24","SG_QG26","SG_QG30","SG_QG32","SG_QG38","SG_QG40","SG_QG44","SG_QG46","SG_QG48","SG_QG50","SG_QG52","SG_QG56","SG_QG60","SG_QG64","SG_QG68","SG_QG70"]
set_FORM_option reference_lmb {\
"SG_QG0":['pq2','pq5'],\
"SG_QG2":['pq2','pq5'],\
"SG_QG4":['pq1','pq4'],\
"SG_QG8":['pq1','pq4'],\
"SG_QG24":['pq3','pq5'],\
"SG_QG26":['pq3','pq5'],\
"SG_QG30":['pq3','pq5'],\
"SG_QG32":['pq3','pq5'],\
"SG_QG38":['pq4','pq6'],\
"SG_QG40":['pq4','pq6'],\
"SG_QG44":['pq4','pq6'],\
"SG_QG46":['pq4','pq6'],\
"SG_QG48":['pq2','pq7'],\
"SG_QG50":['pq2','pq7'],\
"SG_QG52":['pq1','pq7'],\
"SG_QG56":['pq1','pq7'],\
"SG_QG60":['pq1','pq7'],\
"SG_QG64":['pq1','pq7'],\
"SG_QG68":['pq2','pq7'],\
"SG_QG70":['pq2','pq7']\
}
set_alphaLoop_option checkpoint_lvl 2
qgraf_generate  a a > t t~ / u d s c b QCD^2==1 QED^2==2 [QCD]
!rm -rf LU_aa_ttx_NLO_all_aligned
output qgraf LU_aa_ttx_NLO_all_aligned
