import inspect
import unittest
from . import ProcessTester, _aL_dir, _mg5_dir, ProcessTesterError
import multiprocessing

class a_ddx_LO(ProcessTester,unittest.TestCase):

    _generation_card = """
import model aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':0}
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s

set_FORM_option cores 1
set_FORM_option FORM_parallel_cores 1
set_FORM_option FORM_use_max_mem_fraction 0.2

set_alphaLoop_option qgraf_model SM
set_FORM_option generate_integrated_UV_CTs True
set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF

qgraf_define j = d d~ g gh gh~
qgraf_generate a > j j / u s c b t QCD^2==0 QED^2==1
output qgraf %(output_path)s
        """
    
    _global_hyperparameters = {
        'CrossSection.picobarns' : False,
        'CrossSection.incoming_momenta' : [[1.,0.,0.,0.],],
        'CrossSection.m_uv_sq' : 1.,
        'CrossSection.mu_r_sq' : 1.,
    }

    def test_local_evaluations_without_deformation(self):
        super(a_ddx_LO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "none",
            },
            default_xs = (0.1,0.2,0.3),
            targets = {
                'SG_QG0' : [ (0.1,0.2,0.3), 1.0, -0.0019998761241987392+0.j ],
            }
        )

    def test_local_evaluations_with_deformation(self):
        super(a_ddx_LO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "fixed",
            },
            default_xs = (0.1,0.2,0.3),
            targets = {
                'SG_QG0' : [ (0.1,0.2,0.3), 1.0, -0.0019998761241987392+0.j ],
            }
        )

class a_ddx_NLO(ProcessTester,unittest.TestCase):

    _generation_card = """
import model aL_sm-no_widths
set_alphaLoop_option perturbative_orders {'QCD':2}
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s

set_FORM_option cores 2
set_FORM_option FORM_parallel_cores 1
set_FORM_option FORM_use_max_mem_fraction 0.2

set_alphaLoop_option qgraf_model SM
set_FORM_option generate_integrated_UV_CTs True
set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF

qgraf_define j = d d~ g gh gh~
qgraf_generate a > j j / u s c b t [QCD] QCD^2==1 QED^2==1
output qgraf %(output_path)s
        """
    
    _global_hyperparameters = {
        'CrossSection.picobarns' : False,
        'CrossSection.incoming_momenta' : [[1.,0.,0.,0.],],
        'CrossSection.m_uv_sq' : 1.,
        'CrossSection.mu_r_sq' : 1.,
    }

    def test_local_evaluations_without_deformation(self):
        super(a_ddx_NLO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "none",
            },
            default_xs = (0.1,0.2,0.3,0.4,0.5,0.6),
            targets = 
{'SG_QG0': [(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 1, (6.853163409073177e-05+0j)],
 'SG_QG1': [(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 2.0, (-8.666839338579953e-07+0j)]}
        )

    def test_local_evaluations_with_deformation(self):
        super(a_ddx_NLO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "fixed",
            },
            default_xs = (0.1,0.2,0.3,0.4,0.5,0.6),
            targets = 
{'SG_QG0': [(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),1,(1.726694894421919e-05-5.3341651518304635e-05j)],
 'SG_QG1': [(0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 2.0, (-8.666839338585902e-07+0j)]}
        )

    def test_UV(self):
        super(a_ddx_NLO,self).run_UV_tests(warmup_cmds=
"""
set_hyperparameter General.deformation_strategy "none"
set_hyperparameter Selectors []
""", 
            uv_test_cmd = "uv_profile ALL -nsof -min 1.0e3 -max 1.0e5 --n_points 100 --f128 --seed 1"
        )

    def test_IR(self):
        super(a_ddx_NLO,self).run_IR_tests(warmup_cmds=
"""
set_hyperparameter General.deformation_strategy "none"
set_hyperparameter Selectors []
""", 
            ir_test_cmd = "ir_profile ALL -nsof --perturbative_order 1 --ir_limits cutkosky --f128 --min 1.0e-03 --max 1.0e-05 --n_points 30 --seed 1"
        )

    def test_deformation(self):
        super(a_ddx_NLO,self).run_deformation_tests(warmup_cmds=
"""
set_hyperparameter General.deformation_strategy "fixed"
set_hyperparameter Selectors []
""", 
            deformation_test_cmd = "deformation_profile ALL -nsof -nsf -max 1.0e-05 --seed 1"
        )


class a_ddx_NNLO(ProcessTester,unittest.TestCase):

    _generation_card = """
import model aL_sm-no_widths
set_alphaLoop_option perturbative_orders {{'QCD':4}}
set_alphaLoop_option FORM_processing_output_format c
set_alphaLoop_option apply_graph_isomorphisms True
set_alphaLoop_option qgraf_template_model no_s

set_FORM_option cores {n_cores:d}
set_FORM_option FORM_parallel_cores 2
set_FORM_option FORM_use_max_mem_fraction 0.2

set_alphaLoop_option qgraf_model SM
set_FORM_option generate_integrated_UV_CTs True
set_alphaLoop_option FORM_compile_optimization 1
set_alphaLoop_option FORM_integrand_type PF

qgraf_define j = d d~ g gh gh~
qgraf_generate a > j j / u s c b t [QCD QCD] QCD^2==2 QED^2==1
output qgraf %(output_path)s
        """.format(n_cores=multiprocessing.cpu_count())
    
    _global_hyperparameters = {
        'CrossSection.picobarns' : False,
        'CrossSection.incoming_momenta' : [[1.,0.,0.,0.],],
        'CrossSection.m_uv_sq' : 1.,
        'CrossSection.mu_r_sq' : 1.,
    }

    def test_local_evaluations_without_deformation(self):
        super(a_ddx_NNLO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "none",
            },
            default_xs = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
            targets = {
                'SG_QG0': [None, 1, (1.2878730433862703e-07+0j)],
                'SG_QG1': [None, 1, (1.2956260667898326e-07+0j)],
                'SG_QG13': [None, 1, (9.85427413596499e-06+0j)],
                'SG_QG14': [None, 1, (1.6066129308132972e-06+0j)],
                'SG_QG15': [None, 1, (-7.066347023176101e-07+0j)],
                'SG_QG16': [None, 2.0, (-1.7775582075698923e-06+0j)],
                'SG_QG17': [None, 2.0, (4.050862983405686e-07+0j)],
                'SG_QG20': [None, 2.0, (-3.147244169458765e-08+0j)],
                'SG_QG22': [None, 2.0, (1.5244279669489691e-06+0j)],
                'SG_QG23': [None, 2.0, (-2.7354919487919653e-07+0j)],
                'SG_QG24': [None, 2.0, (-4.4059230673593037e-08+0j)],
                'SG_QG25': [None, 2.0, (4.356407266547342e-07+0j)],
                'SG_QG3': [None, 1, (-1.5187822308664218e-06+0j)],
                'SG_QG5': [None, 2.0, (7.958636305485362e-08+0j)],
                'SG_QG6': [None, 2.0, (6.277272025452023e-07+0j)],
                'SG_QG9': [None, 4.0, (-2.9815452531633664e-07+0j)]
            }
        )

    @unittest.skip("Deformation currently bugged at NNLO")
    def test_local_evaluations_with_deformation(self):
        super(a_ddx_NNLO,self).run_local_evaluations(
            test_name = inspect.stack()[0][3],
            extra_hyperparameters = {
                'General.deformation_strategy': "fixed",
            },
            default_xs = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
            targets = {} # Will need to be updated when deformation will be fixed
        )