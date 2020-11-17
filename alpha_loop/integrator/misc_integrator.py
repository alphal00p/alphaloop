################################################################################
#
# Copyright (c) 2017 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

import numpy as np

def check_input_types(continuous_inputs, discrete_inputs):
    """ Wrapper to easily check input types."""
    assert(isinstance(discrete_inputs, np.ndarray))
    assert(not discrete_inputs or isinstance(i, np.int64) for i in discrete_inputs)
    assert(isinstance(continuous_inputs, np.ndarray))
    assert(not continuous_inputs or isinstance(f, numpy.float64) for f in continuous_inputs)        
    return True

