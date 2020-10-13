import model aL_sm-no_widths
set_alphaLoop_option n_rust_inputs_to_generate -1
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s
set_alphaLoop_option qgraf_model SM_HEFT
set_alphaLoop_option loop_induced False 
set_alphaLoop_option qgraf_cut_filter False 
set_FORM_option number_of_lmbs None
set_FORM_option generate_integrated_UV_CTs True 
#qgraf_define j = d d~ g gh gh~
qgraf_generate h > d d~ g / u s c b QCD^2==3 QED^2==1 []
!rm -rf h_jj_LO_HEFT_ddxg
output qgraf h_jj_LO_HEFT_ddxg
