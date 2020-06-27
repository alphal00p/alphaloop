import model sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':2}
# Set the option below to -1 in order to generate all rust inputs as well
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
set_alphaLoop_option consider_vertex_id_in_graph_isomorphism False 
set_alphaLoop_option consider_edge_orientation_in_graph_isomorphism False
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
prepare_model_for_self_energies 2 ((22,6,-6),)
set_FORM_option number_of_lmbs None
set_FORM_option reference_lmb None
define SE_1 = SE_1_g SE_1_t SE_1_t~
define SE_2 = SE_2_g SE_2_t SE_2_t~
define SE_start = t t~ g 
# Generate the non-self energy contributions
# ------------------------------------------
#
# We can ignore the qq~ final state since it is only from self-energies.
# 
add process e+ e- > a > t t~ h g QCD=1 @0
#
# Now generate the self-energies
# ------------------------------
#
# Self-energies in real-virtual contributions, off the top quark, anti-top quark and gluon
#
add process e+ e- > a > t t~ h SE_start SE_1 SE_1 se_1_a_1 QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @40
#
# Finally output that glorious mess
#
!rm -rf TEST_httx_NLO_auto_SEs
output alphaLoop TEST_httx_NLO_auto_SEs
