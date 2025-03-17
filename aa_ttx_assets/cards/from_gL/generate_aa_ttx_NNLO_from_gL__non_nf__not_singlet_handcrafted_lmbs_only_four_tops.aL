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
set_qgraf_spoof ./PLUGIN/new_aL_deformation/aa_ttx_assets/graphs/NNLO_from_gL__non_nf__not_singlet.py
set_FORM_option reference_lmb {\
"GL584":['pq2','pq4','pq8'],\
"GL588":['pq2','pq4','pq10'],\
"GL594":['pq2','pq4','pq8'],\
"GL620":['pq2','pq4','pq10'],\
"GL860":['pq1','pq4','pq10'],\
"GL862":['pq1','pq4','pq10'],\
"GL514":['pq5','pq4','pq2'],\
"GL511":['pq2','pq4','pq10'],\
"GL493":['pq2','pq4','pq5'],\
"GL491":['pq2','pq4','pq10'],\
"GL470":['pq2','pq4','pq10'],\
"GL453":['pq2','pq4','pq8'],\
"GL448":['pq2','pq4','pq5'],\
"GL442":['pq2','pq4','pq10'],\
"GL440":['pq2','pq4','pq10'],\
"GL434":['pq2','pq4','pq10'],\
"GL428":['pq2','pq4','pq10'],\
"GL410":['pq2','pq4','pq10'],\
"GL409":['pq2','pq4','pq10'],\
"GL374":['pq2','pq4','pq7'],\
"GL372":['pq2','pq4','pq7'],\
"GL370":['pq2','pq4','pq7'],\
"GL368":['pq2','pq4','pq7'],\
"GL316":['pq3','pq4','pq10'],\
"GL313":['pq3','pq4','pq10'],\
"GL290":['pq3','pq4','pq8'],\
"GL284":['pq3','pq4','pq10'],\
"GL280":['pq3','pq4','pq8'],\
"GL244":['pq3','pq4','pq10'],\
"GL236":['pq3','pq4','pq10'],\
"GL235":['pq3','pq4','pq10'],\
"GL226":['pq3','pq4','pq10'],\
"GL221":['pq3','pq4','pq10'],\
"GL86":['pq2','pq5','pq10'],\
"GL69":['pq2','pq5','pq10'],\
"GL336":['pq5','pq4','pq6'],\
"GL340":['pq5','pq4','pq6'],\
"GL344":['pq4','pq6','pq10'],\
"GL348":['pq4','pq6','pq10'],\
"GL614":['pq2','pq6','pq10'],\
"GL627":['pq2','pq6','pq10'],\
}
#define j = j t t~
qgraf_generate  a a > t t~ t t~ / u d s c b QCD^2==1 QED^2==2 [QCD]
!rm -rf ./PLUGIN/new_aL_deformation/aa_ttx_assets/runs/NNLO_from_gL__non_nf__not_singlet_handcrafted_lmbs_only_four_tops
output qgraf ./PLUGIN/new_aL_deformation/aa_ttx_assets/runs/NNLO_from_gL__non_nf__not_singlet_handcrafted_lmbs_only_four_tops
