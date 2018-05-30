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

 Vegas3
 numpy


2. Usage
--------

 The pyNLoop plugin can be started by running:
   
   `./bin/mg5_aMC --mode=pyNLoop`

 from within your MadNkLO distribution. You can then test the following
 new command added by the pyNLoop plugin:

   `pyNLoop > hello_world something`
   `pyNLoop.Interface: Hello World: something`
