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

 A simple integration of a dummy function is already implemented:

`pyNLoop > integrate_loop dummy_function <options>`

 where each element of <options> if of the form `--<option_name>=<options_value>' 
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
            'phase_computed'        : 'Real',
            'survey_n_iterations'   : 10,
            'survey_n_points'       : 2000,
            'refine_n_iterations'   : 10,
            'refine_n_points'       : 1000
        }
```

Example of a quick run:

```
pyNLoop > integrate_loop dummy_function --nb_CPU_cores=3 --survey_n_points=1000 --survey_n_iterations=5 --refine_n_points=1000 --refine_n_iterations=5 --verbosity=0
pyNLoop.Interface: Using 3 CPU cores
pyNLoop.Interface: ====================================================================================================
pyNLoop.Interface:                             Starting integration, lay down and enjoy...
pyNLoop.Interface: ====================================================================================================
pyNLoop.Interface: ====================================================================================================
pyNLoop.Interface:                            Loop integral result with integrator 'Vegas3':
pyNLoop.Interface:                                       1.70723e+00 +/- 6.42e-03
pyNLoop.Interface: ====================================================================================================
```
