import model ./PLUGIN/alphaloop/models/aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':4}
# Set the option below to -1 in order to generate all rust inputs as well
set_alphaLoop_option n_rust_inputs_to_generate -1
# Set the option below to False as for self-energies it is not correctly supported yet
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
prepare_model_for_self_energies 2 ((22,6,-6),)
# 
# Self-energy corresponding to two-loop diagrams off the top quark and anti-top quark:
#
add process e+ e- > a > t t~ h t SE_1_t~ SE_1_g SE_1_g se_1_a_1 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @10
add process e+ e- > a > t t~ h t~ SE_1_t SE_1_g SE_1_g se_1_a_1 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @11
#
# Finally output that glorious mess
#
!rm -rf EXAMPLE_httx_2LSE
output alphaLoop EXAMPLE_httx_2LSE
