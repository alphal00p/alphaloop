#!/usr/bin/env python2
import os
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath( __file__ )),os.path.pardir))
import vectors 

class Deformation(object):

    def __init__(self, topology, **opts):
        self.topology = topology

    def __call__(self, loop_momenta, *args, **opts):
        raise BaseException("This function should be implemented in a daughter class.") 

class Contribution(object):

    def __init__(self, identifier, *args, **opts):
        self.identifier = identifier
    
    def __call__(self, PS_point, loop_momenta, *args, **opts):
        raise BaseException("This function should be implemented in a daughter class.")

    @staticmethod
    def prop(p, m_squared):
        return 1./(p.square()-m_squared)

    @staticmethod
    def inv_prop(p, m_squared):
        return p.square()-m_squared

class Diagram(Contribution):

    def __init__(self, identifier, parameters, **opts):
        super(Diagram, self).__init__(identifier, **opts)
        self.parameters = parameters

    def __call__(self, PS_point, loop_momenta, *args, **opts):
        raise BaseException("This function should be implemented in a daughter class.")

class Counterterm(Contribution):

    def __init__(self, identifier, parameters, **opts):
        super(Counterterm, self).__init__(identifier, **opts)
        self.parameters = parameters

    def __call__(self, PS_point, loop_momenta, *args, **opts):
        raise BaseException("This function should be implemented in a daughter class.")

class Amplitude(Contribution):

    def __init__(self, identifier, n_loops=1, **opts):
        super(Amplitude, self).__init__(identifier, **opts)
        self.n_loops = n_loops

        self.parameters = {
            'g_f'       :   1.166390e-05,
            'alpha_s'   :   1.180000e-01,
            'alpha_ew'  :   1./1.325070e+02,
            'C_F'       :   4./3.,
            'm_IR2'     :   100.
        }

        for opt in opts:
            if opt in self.parameters:
                self.parameters[opt] = opts[opt]

    def __call__(self, PS_point, loop_momenta):
        raise BaseException("This function should be implemented in a daughter class.")

if __name__ == "__main__":

    import QQAA
    
    _PERFORM_INTEGRATION = False 
    _PERFORM_TEST_POINT  = True
    _N_LOOPS             = 1
    _PROCESS             = 'qq_aa'
    amp = QQAA.Amplitude(_PROCESS, n_loops=_N_LOOPS)

    PS_point = [
        vectors.LorentzVector([5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00,  5.0000000000000000e+02]),
        vectors.LorentzVector([5.0000000000000000e+02,  0.0000000000000000e+00,  0.0000000000000000e+00, -5.0000000000000000e+02]),
        vectors.LorentzVector([5.0000000000000000e+02,  1.1092428444383280e+02,  4.4483078948812141e+02, -1.9955292993087880e+02]),
        vectors.LorentzVector([5.0000000000000000e+02, -1.1092428444383280e+02, -4.4483078948812141e+02,  1.9955292993087869e+02])
    ]


    # Example of a test_point evaluation
    if _PERFORM_TEST_POINT:

        loop_momenta = [
            vectors.LorentzVector([1.0000000000000000e+02,  2.0000000000000000e+02,  3.0000000000000000e+02,  4.0000000000000000e+02]),
        ] 

        print(amp(PS_point, loop_momenta))

    # Example of integration
    if _PERFORM_INTEGRATION:
        
        try:
            import vegas 
        except ImportError:
            print("WARNING :: Vegas 3 could not be found, install it if needed.")
            sys.exit(1)
        try:
            import numpy as np
        except ImportError:
            print("ERROR :: Vegas 3 requires numpy, install it first.")
            sys.exit(1)
        try:
            import gvar as gv
        except ImportError:
            print("WARNING :: Vegas 3 requires gvar (gaussian variable handling package), install it if needed.")
            sys.exit(1)
       
        
        n_dimensions = 3*_N_LOOPS
        vegas3_integrator = vegas.Integrator(n_dimensions * [[0., 1.]],
                analyzer        = vegas.reporter(), 
                nhcube_batch    = 1000,
                # We must disable the hypercube optimisation for now when plotting observables
                max_nhcube      = 1e9,
                beta            = 0.75,
                alpha           = 0.5,
                sync_ran        = True
        )

        import parametrisation
        parametriser = parametrisation.CartesianParametrisation(n_loops=1)
        import LTD
        my_topology = LTD.LoopTopology(
            n_loops     = _N_LOOPS,
            loop_lines  =   [
                LTD.LoopLine(
                    loop_momenta = (LTD.POSITIVE_FLOW,),
                    q            = PS_point[1]+PS_point[2],
                    m_squared    = 0.
                ),
                LTD.LoopLine(
                    loop_momenta = (LTD.POSITIVE_FLOW,),
                    q            = PS_point[1]+PS_point[3],
                    m_squared    = 0.
                ),
                # TODO: etc...
            ]
        )
        deformer = LTD.LTDDeformation(my_topology)

        LTD_integrand_function = LTD.LTDIntegrand(
            my_topology, parametriser, deformer, amp, PS_point)
    
        result = vegas3_integrator(
            LTD_integrand_function,     # Function to integrate
            nitn    =   10,             # Number of iterations
            neval   =   10000            # Number of evaluations per iteration
        )

        print("LTD integration result for process '%s':\n%s"%(_PROCESS,str(result)))
