## import the required files
import os
import sys

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, root_path)

import alpha_loop.interface as interface
import alpha_loop.exporters as exporters

# Apply patches to some functions of madgraph to suit them to the LTD^2 programme.
# Patch #1: Turn off the optimization of wavefunction recycling when writing matrix.f
# Patch #2: Remove the denominator of all propagators in aloha. 
import alpha_loop.madgraph_patches as madgraph_patches
madgraph_patches.remove_wavefunction_recycling()
madgraph_patches.remove_propagator_denominators()
madgraph_patches.force_aloha_to_loop_mode()
madgraph_patches.enable_manual_complex_conjugation_in_aloha()

# Three types of functionality are allowed in a plugin
#   1. new output mode
#   2. new cluster support
#   3. new interface

# 1. Define new output mode
#    example: new_output = {'myformat': MYCLASS}
#    madgraph will then allow the command "output myformat PATH"
#    MYCLASS should inherated of the class madgraph.iolibs.export_v4.VirtualExporter 
new_output = {
    'alphaLoop' : exporters.alphaLoopExporter
}

# 2. Define new way to handle the cluster.
#    example new_cluster = {'mycluster': MYCLUSTERCLASS}
#    allow "set cluster_type mycluster" in madgraph
#    MYCLUSTERCLASS should inherated from madgraph.various.cluster.Cluster
new_cluster = {}

# 3. Define a new interface (allow to add/modify MG5 command)
#    This can be activated via ./bin/mg5_aMC --mode=PLUGINNAME
## Put None if no dedicated command are required
new_interface = interface.alphaLoopInterface
 
 
########################## CONTROL VARIABLE ####################################
__author__ = 'Valentin Hirschi'
__email__ = 'valentin.hirschi@gmail.com'
__version__ = (1,0,0)

#
# To be updated later
#

minimal_mg5amcnlo_version = (2,0,0) 
maximal_mg5amcnlo_version = (1000,1000,1000)
latest_validated_version = (2,0,0)
