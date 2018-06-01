#                                                                    
#                       PyNLoop MG5aMC plugin                        
#                                                                    
# Numerical integration of multi-loop integrals steered with Python  
#                                                                    


1. Installation and requirements
--------------------------------

 The PyNLoop directory must be placed int PLUGIN folder of your 
 distribution of the MadNkLO branch of MG5aMC.
 ( bzr branch lp:~madnklo/mg5amcnlo/MadNkLO )
 
 The PLUGIN's python dependencies are:

`python2.7+`
`Vegas3`
`numpy`
`bidict`
`sympy`
`pySecDec v1.3.1+`


2. Usage
--------

 The pyNLoop plugin can be started by running:
   
   `./bin/mg5_aMC --mode=pyNLoop`

 from within your MadNkLO distribution. You can then test the following
 new command added by the pyNLoop plugin:

```
pyNLoop > hello_world something
pyNLoop.Interface: Hello World: something
```

 The integration can be launched with the `integrate_loop` command:

`pyNLoop > integrate_loop <integrand_specification> <options>`

 where each element of `<options>` if of the form `--<option_name>=<options_value>' 
 with the following default key-value pairs:

```
        options = { 
            'PS_point'              : 'random',
            'seed'                  : None,
            'sqrt_s'                : 1000.0,
            'target_accuracy'       : 1.0e-3,
            'batch_size'            : 1000,
            'verbosity'             : 1,
            'nb_CPU_cores'          : None,
            'phase_computed'        : 'All',
            'survey_n_iterations'   : 10,
            'survey_n_points'       : 2000,
            'refine_n_iterations'   : 10,
            'refine_n_points'       : 1000,
            'output_folder'         : pjoin(MG5DIR,'MyPyNLoop_output'),
			# This is to bypass questions, for example related to overwriting existing output
            'force'                 : False
        }
```

The value of `<integrand_specification>` are currently only the two harcoded example: `dummy_function` and `box1L`.

The first one is integrated using Vegas3:

```
pyNLoop > integrate_loop dummy_function --nb_CPU_cores=3 --survey_n_points=1000 --survey_n_iterations=5 --refine_n_points=1000 --refine_n_iterations=5 --verbosity=0 --output_folder=DummyExample
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface:                                                      Starting integration, lay down and enjoy...
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface:                                              Integral of 'DummyNLoopIntegrand' with integrator 'Vegas3':
pyNLoop.Interface:
pyNLoop.Interface:                (1.7003451099941695                                )     +/- (0.0063324417793653581                             )
pyNLoop.Interface:              + (0                                                 ) ε⁻¹ +/- (0                                                 ) ε⁻¹
pyNLoop.Interface:              + (0                                                 ) ε⁻² +/- (0                                                 ) ε⁻²
pyNLoop.Interface:
pyNLoop.Interface: ======================================================================================================================================================
```

while the second one reproduces the canonical example of pySecDec:

```
pyNLoop > integrate_loop box1L --verbosity=0 --output_folder=pySecDecExample --force
pyNLoop.Integrand: Generating pySecDec output for integrand 'box1L' ...
pyNLoop.pySecDecIntegrator: Compiling pySecDec output of integrand 'box1L' ...
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface:                                                      Starting integration, lay down and enjoy...
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface: ======================================================================================================================================================
pyNLoop.Interface:                                               Integral of 'box1L' with integrator 'pySecDecIntegrator':
pyNLoop.Interface:
pyNLoop.Interface:                (-0.425518029321860614 + 1.86892573491243086*I     )     +/- (0.00706836320272220993 + 0.0186494646389539255*I  )
pyNLoop.Interface:              + (0.639402100255728412 + 1.95699136112339149e-6*I   ) ε⁻¹ +/- (0.00650732337473908327 + 0.000971486097553553372*I) ε⁻¹
pyNLoop.Interface:              + (-0.14287722321378088 - 1.8115437022083544e-7*I    ) ε⁻² +/- (0.00118488128199477761 + 0.000210769586888484833*I) ε⁻²
pyNLoop.Interface:
pyNLoop.Interface: ======================================================================================================================================================
```
