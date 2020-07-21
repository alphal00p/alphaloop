import model sm-no_widths
# Set LO cuts only
set_alphaLoop_option perturbative_orders {'QCD':2}
#set_alphaLoop_option n_jets 0
#set_alphaLoop_option final_state_pdgs []
# Also enable the generation of Rust yaml inputs
set_alphaLoop_option n_rust_inputs_to_generate -1
# For LO comparison we can enable the isomorphism accounting for particles and anti-particles
set_alphaLoop_option differentiate_particle_from_antiparticle_in_graph_isomorphism False 
# For LO we can also safely include self-energies
set_alphaLoop_option include_self_energies_from_squared_amplitudes True
# And finally let us also enable the FORM processing of numerators
set_alphaLoop_option FORM_processing_output_format c
#set_alphaLoop_option supergraphs_selection [(1,2),(2,2)]
set_FORM_option number_of_lmbs -1
set_FORM_option reference_lmb 2 

set_alphaLoop_option FORM_integrand_type PF
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_compile_optimization 3

generate e+ e- > a > d d~ g
!rm -rf TEST_epem_a_ddxg_NLO
output alphaLoop TEST_epem_a_ddxg_NLO
