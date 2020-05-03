#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import madgraph.iolibs.helas_call_writers as helas_call_writers

class alphaLoopHelasCallWriter(helas_call_writers.FortranUFOHelasCallWriter):
    """ HelasCallWriter for alphaLoop with functions that can be specialised as needed."""

    def __init__(self, *args, **opts):

        super(alphaLoopHelasCallWriter, self).__init__(*args, **opts)