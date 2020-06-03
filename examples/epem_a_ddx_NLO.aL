import model aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':2}
set_alphaLoop_option n_rust_inputs_to_generate -1
# Set the option below to False as for self-energies it is not correctly supported yet
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False
prepare_model_for_self_energies 2 ((22,1,-1),)
# The double-triangle topology
add process e+ e- > a > d d~ g QCD=1 SE_1_QCD=0 SE_2_QCD=0 @1
# The self-energy topology attached to the anti-down quark
add process e+ e- > a > d d~ d~ SE_1_g SE_1_d se_1_a_1 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @2
# The self-energy topology attached to the down quark
add process e+ e- > a > d d~ d SE_1_g SE_1_d~ se_1_a_1 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @3
rm -rf TEST_epem_a_ddx_NLO
output alphaLoop TEST_epem_a_ddx_NLO
