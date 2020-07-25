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
set_alphaLoop_option apply_graph_isomorphisms True
# And finally let us also enable the FORM processing of numerators
set_alphaLoop_option FORM_processing_output_format c
#set_alphaLoop_option supergraphs_selection [(1,2),(2,2)]
set_FORM_option number_of_lmbs 0
set_FORM_option reference_lmb 6
generate e+ e- > a > d d~ g
!rm -rf epem_a_ddx_paper_rendering 
output alphaLoop epem_a_ddx_paper_rendering
!echo "\n\n\n\033[31mDO NOT FORGET TO ADJUST THE MULTIPLICITY OF THE SELF-ENERGY SUPERGRAPHS!\033[0m\n\n\n"
