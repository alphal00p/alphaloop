#####################################################
#                                                   #
#  Source file of the alphaLoop MG5aMC plugin.      #
#                                                   #
#####################################################

import madgraph.iolibs.helas_call_writers as helas_call_writers
import madgraph.core.helas_objects as helas_objects
import madgraph.loop.loop_helas_objects as loop_helas_objects


class alphaLoopHelasCallWriter(helas_call_writers.FortranUFOHelasCallWriter):
    """ HelasCallWriter for alphaLoop with functions that can be specialised as needed."""

    def __init__(self, *args, **opts):
        super(alphaLoopHelasCallWriter, self).__init__(*args, **opts)

    def get_amplitude_calls(self, matrix_element,*args, **opts):
        """Return a list of strings, corresponding to the Helas calls
        for the matrix element"""
        
        # This function does not seem to need to be re-implemented.
        return super(alphaLoopHelasCallWriter, self).get_amplitude_calls(
            matrix_element,*args, **opts)

    def get_wavefunction_call(self, wavefunction):
        """Customize Helas wavefunction call for LTD^2"""

        call = super(alphaLoopHelasCallWriter, self).get_wavefunction_call(wavefunction)

        # Special treatment for external final_state legs
        if len(wavefunction.get('mothers'))==0 and wavefunction.get('leg_state')==True:
            # This is only needed for non-scalar
            if wavefunction.get('spin') in [2,3]:
                # Prefix the polarisation vector routine differently when on 
                # the left or the right of the Cutkosky cut so as to be able to emulate
                # any propagator
                call = call.split(',')
                call = ','.join(call[:-2]+['K']+call[-2:])
                call = call.replace('CALL ','CALL PROP_')

        return call

    def get_matrix_element_calls(self, matrix_element):
        """Return a list of strings, corresponding to the Helas calls
        for the matrix element"""

        assert isinstance(matrix_element, helas_objects.HelasMatrixElement), \
                  "%s not valid argument for get_matrix_element_calls" % \
                  type(matrix_element)

        # Do not reuse the wavefunctions for loop matrix elements
        if isinstance(matrix_element, loop_helas_objects.LoopHelasMatrixElement):
            return self.get_loop_matrix_element_calls(matrix_element)
        
        me = matrix_element.get('diagrams')
        matrix_element.reuse_outdated_wavefunctions(me)

        res = []
        for diagram in matrix_element.get('diagrams'):
            res.append("%d CONTINUE"%diagram.get('number'))
            res.extend([ self.get_wavefunction_call(wf) for \
                         wf in diagram.get('wavefunctions') ])
            res.append("# Amplitude(s) for diagram number %d" % \
                       diagram.get('number'))
            for amplitude in diagram.get('amplitudes'):
                res.append(self.get_amplitude_call(amplitude))
            res.append("GOTO 9999")

        return res