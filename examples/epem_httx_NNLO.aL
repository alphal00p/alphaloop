import model ./PLUGIN/alphaloop/models/aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':4}
# Set the option below to -1 in order to generate all rust inputs as well
set_alphaLoop_option n_rust_inputs_to_generate -1
# Set the option below to False as for self-energies it is not correctly supported yet
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism True 
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
# Self-energy corresponding to two-loop diagrams off the top quark and anti-top quark:
#
add process e+ e- > a > t t~ h t SE_1_g SE_1_g SE_1_t~ se_1_a_1 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @10
add process e+ e- > a > t t~ h t~ SE_1_g SE_1_g SE_1_t se_1_a_1 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @11
#
# Next, we want the case of iterated 1-loop self-energies (in this case 2, at NNLO), off the top quark twice, or anti-top quark twice or once top, once anti-top
#
add process e+ e- > a > t t~ h t SE_1_g SE_1_t~ t SE_2_g SE_2_t~ se_1_a_1 se_2_a_1 QCD=0 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 SE_1_QCD=1 SE_2_CONTACT_1=1 SE_2_BRIDGE_1=1 SE_2_ANCHOR=1 @20
add process e+ e- > a > t t~ h t~ SE_1_g SE_1_t t~ SE_2_g SE_2_t se_1_a_1 se_2_a_1 QCD=0 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 SE_1_QCD=1 SE_2_CONTACT_1=1 SE_2_BRIDGE_1=1 SE_2_ANCHOR=1 @21
add process e+ e- > a > t t~ h t SE_1_g SE_1_t~ t~ SE_2_g SE_2_t se_1_a_1 se_2_a_1 QCD=0 SE_1_QCD=1 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 SE_1_QCD=1 SE_2_CONTACT_1=1 SE_2_BRIDGE_1=1 SE_2_ANCHOR=1 @22
#
# The case of one 2-loop SE, with a nested 1-loop SE within, off the top and anti-top quark.
# Note: this should in principle also generate the self-energy correction of the gluon from closed bubble top quark loop.
# 
add process e+ e- > a > t t~ h t SE_1_g SE_1_t~ SE_1_g SE_1_t SE_1_t~ se_1_a_2 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_CONTACT_2=1 SE_1_BRIDGE_2=1 SE_1_ANCHOR=1 @30
add process e+ e- > a > t t~ h t SE_1_g SE_1_t~ SE_1_g SE_1_g SE_1_g se_1_a_2 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_CONTACT_2=1 SE_1_BRIDGE_2=1 SE_1_ANCHOR=1 @31
add process e+ e- > a > t t~ h t~ SE_1_g SE_1_t SE_1_g SE_1_t SE_1_t~ se_1_a_2 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_CONTACT_2=1 SE_1_BRIDGE_2=1 SE_1_ANCHOR=1 @32
add process e+ e- > a > t t~ h t~ SE_1_g SE_1_t SE_1_g SE_1_g SE_1_g se_1_a_2 QCD=0 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_CONTACT_2=1 SE_1_BRIDGE_2=1 SE_1_ANCHOR=1 @33
#
# Self-energies in real-virtual contributions, off the top quark, anti-top quark and gluon
#
add process e+ e- > a > t t~ h g t SE_1_g SE_1_t~ se_1_a_1 QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @40
add process e+ e- > a > t t~ h g t~ SE_1_g SE_1_t se_1_a_1 QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @41
add process e+ e- > a > t t~ h g g SE_1_g SE_1_g se_1_a_1 QCD=1 SE_1_QCD=2 SE_1_CONTACT_1=1 SE_1_BRIDGE_1=1 SE_1_ANCHOR=1 @42
#
# Finally output that glorious mess
#
!rm -rf EXAMPLE_httx_NNLO
output alphaLoop EXAMPLE_httx_NNLO
