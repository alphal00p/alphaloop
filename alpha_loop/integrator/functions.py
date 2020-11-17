################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
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


import os
import logging
import math
import cmath
import random
import shutil
import numpy as np

import madgraph

MADEVENT= False
import madgraph.various.misc as misc
import madgraph.iolibs.files as files
import madgraph.various.cluster as cluster
import madgraph.various.lhe_parser as lhe_parser
import alpha_loop.integrator.misc_integrator as misc_integrator

logger = logging.getLogger('madgraph.integrator')
pjoin = os.path.join

class VirtualFunction(object):
    """A mother base class that specified mandatory feature that any function block should implement."""
        
    def __call__(self, continuous_inputs, discrete_inputs, **opts):
        """ Integrand function call, with list of continuous and discrete input values for all dimensions."""
        self.check_input_types(continuous_inputs, discrete_inputs)
        
        # A unique float must be returned, but additional return values can 
        # be specified in the body of the function or via options.
        return 0.0

    def check_input_types(self, *args, **opts):
        return misc_integrator.check_input_types(*args, **opts)

class FunctionList(list):
    """A mother base class that specified mandatory feature that any observable list should implement."""

    def __init__(self, *args, **opts):
        super(FunctionList, self).__init__(*args, **opts)

    def append(self, arg, **opts):
        """ Type-checking. """
        assert(isinstance(arg, VirtualFunction))
        super(FunctionList, self).append(arg, **opts)

class FunctionFromPythonExpression(VirtualFunction):
    """A simple function from python expression using the variables ci[k] and di[k] for the continuous and
    discrete inputs respectively"""

    def __init__(self, expr, dimensions=None, **opts):
        self.expr = expr
        if dimensions:
            self.ci_labels = [d.name for d in dimensions.get_continuous_dimensions()]
            self.di_labels = [d.name for d in dimensions.get_discrete_dimensions()]
        else:
            self.ci_labels = None
            self.di_labels = None            
        super(FunctionFromPythonExpression,self).__init__(**opts)

    def __call__(self, ci, di, **opts):
        """ Evaluate Python expression from its expression provided when instantiated
        and the discrete and continous dimensions provided."""
        locals = {'ci':ci,'di':di, 'math': math, 'cmath': cmath, 'random': random,'numpy':np}

        if self.ci_labels:
            for i, input in enumerate(ci):
                locals[self.ci_labels[i]] = input
        if self.di_labels:
            for i, input in enumerate(di):
                locals[self.di_labels[i]] = input

        return {'weight': eval(self.expr, locals)}
    
