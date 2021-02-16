import model aL_sm-no_widths
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option FORM_compile_arg all
set_alphaLoop_option FORM_integrand_type both
set_alphaLoop_option FORM_compile_optimization 3
set_FORM_option extra-options {'OPTIMITERATIONS':1000}
set_FORM_option optimisation_strategy 'CSEgreedy'
!rm -rf SCALAR_PENTAGON_BABIS 
output_scalar_integral SCALAR_PENTAGON_BABIS \
 --topology=\
 [('q1', 1, 101), ('q2', 2, 102),\
  ('p1', 102, 101), ('p2', 101, 103), ('p3', 103, 104), ('p4', 104, 105), ('p5', 105, 102),\  
('pout3', 103, 3), ('pout4', 104, 4), ('pout5', 105, 5)] \
 --name=scalar_pentagon_babis \
 --analytical_result=0.0159245\
 --externals=[('q1','q2'),('pout3','pout4','pout5')] \
 --default_kinematics={'q1':[1.,0.,0.,1.],\
                       'q2':[-1.,0.,0.,1.],\
                       'pout3':[0.,2.3,0.2,1.1],\
                       'pout4':[0.,0.2,1.3,3.123]} \
 --lmb=['p1',] \
 --numerator='1-xZenoConst1*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)+xZenoConst2*(k1.k1+2*k1.p1+p1.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst3*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)+xZenoConst4*(k1.k1+k2.k2+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1)*(k1.k1-2*k1.p2+p2.p2)+xZenoConst5*(k1.k1+k2.k2+k3.k3+p1.p1+2*k1.p1-2*k1.k2-2*k2.p1-2*k3.k1+2*k3.k2-2*k3.p1)*(k1.k1-2*k1.p2+p2.p2)'
