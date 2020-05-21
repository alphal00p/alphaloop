import model ./PLUGIN/alphaloop/models/aL_sm-no_widths
# Set LO cuts only
set_alphaLoop_option perturbative_orders {'QCD':0}
# Also enable the generation of Rust yaml inputs
set_alphaLoop_option n_rust_inputs_to_generate -1
# For LO comparison we can enable the isomorphism accounting for particles and anti-particles
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism True
# For LO we can also safely include self-energies
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
# And finally let us also enable the FORM processing of numerators
set_alphaLoop_option FORM_processing_output_format c 
generate e+ e- > a > d d~
!rm -rf TEST_epem_a_ddx
output alphaLoop TEST_epem_a_ddx
