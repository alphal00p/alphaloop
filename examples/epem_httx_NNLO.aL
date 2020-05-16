import model sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':4}
# Set the option below to -1 in order to generate all rust inputs as well
set_alphaLoop_option n_rust_inputs_to_generate 0
prepare_model_for_self_energies 2 ((22,6,-6),)
# Generate the non-self energy contributions
# ------------------------------------------
#
# We can ignore the qq~ final state since it is only from self-energies.
# 
add process e+ e- > a > t t~ h g g QCD=2 @0
#
# Now generate the self-energies
# ------------------------------
# 
# Self-energy corresponding to two-loop diagrams:
#
add process e+ e- > a > t t~ h se_1_a_1 t SE_1_g SE_1_g SE_1_t~ QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @1
#
# Next, we want the case of iterated 1-loop self-energies (in this case 2, at NNLO):
#
add process e+ e- > a > t t~ h se_1_a_1 t SE_1_g SE_1_t~ se_2_a_1 t SE_2_g SE_2_t~ QCD=0 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 SE_1_QCD=1 SE_2_CONTACT_1=1 SE_2_BRIDGE_1=1 SE_2_ANCHOR=1 @2
#
# The case of one 2-loop SE, with a nested 1-loop SE within:
# 
add process e+ e- > a > t t~ h se_1_a_2 t SE_1_g SE_1_t~ SE_1_g SE_1_t SE_1_t~ QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_CONTACT_2=1 SE_1_BRIDGE_2=1 SE_1_ANCHOR=1 @3
#
# Self-energies in real-virtual contributions, firs of the top quark
#
add process e+ e- > a > t t~ h g se_1_a_1 t SE_1_g SE_1_t~ QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @4
#
# And finally of the gluon
#
add process e+ e- > a > t t~ h g se_1_a_1 g SE_1_g SE_1_g QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @5
#
# Finally output that glorious mess
#
!rm -rf EXAMPLE_httx_NNLO
output alphaLoop EXAMPLE_httx_NNLO
