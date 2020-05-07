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

    def __init__(self, *args, use_physical_gluon_helicity_sum=False, **opts):
        super(alphaLoopHelasCallWriter, self).__init__(*args, **opts)
        self.use_physical_gluon_helicity_sum = use_physical_gluon_helicity_sum

    def add_leftright_arg(self, call):
        """ Add the leftright 'K' argument to a HELAS call."""

        # TODO maybe do something less hacky

        call_elems= call.split('(')
        call_elems = ['('.join(call_elems[:-1]),call_elems[-1]]
        part_one_call_split = call_elems[0].split(',')
        part_one_call_split = ','.join(part_one_call_split[:-1]+['K']+part_one_call_split[-1:])
        return '%s(%s'%(part_one_call_split,call_elems[1])

    def get_amplitude_call(self, amplitude,**opts):
        """ Specialize funxtion as to add leftright suffix"""        

        # Add to the HELAS call the information of whether it sits to the left or 
        # or to the right of the Cutkosky cut so as to be able to apply the
        # complex conjugation "manually" without conjugating complex momenta.
        call = super(alphaLoopHelasCallWriter, self).get_amplitude_call(amplitude,**opts)

        # Add to the HELAS call the information of whether it sits to the left or 
        # or to the right of the Cutkosky cut so as to be able to apply the
        # complex conjugation "manually" without conjugating complex momenta.
        return self.add_leftright_arg(call)

    def get_wavefunction_call(self, wavefunction):
        """Customize Helas wavefunction call for LTD^2"""

        call = super(alphaLoopHelasCallWriter, self).get_wavefunction_call(wavefunction)

        # Add to the HELAS call the information of whether it sits to the left or 
        # or to the right of the Cutkosky cut so as to be able to apply the
        # complex conjugation "manually" without conjugating complex momenta.
        call = self.add_leftright_arg(call)

        # Special treatment for external legs
        if len(wavefunction.get('mothers'))==0:
            # Different treatment for the final and initial-state legs for our
            # current treatment
            if wavefunction.get('leg_state')==True:
                if wavefunction.get('spin')==3 and wavefunction.get('mass').upper()=='ZERO' \
                                    and self.use_physical_gluon_helicity_sum:
                    call = call.replace('CALL ','CALL PROPPHYS_')
                else:
                    if ( wavefunction.get('spin')==1 and (not wavefunction.get('self_antipart')) and wavefunction.get('color')==8 ):
                        # This is a ghost, use external "polarisation" emulating the -1 factor of closed loops
                        if wavefunction.get('is_part'):
                            call = call.replace('CALL SXXXXX','CALL PROP_GHXXXX')
                        else:
                            call = call.replace('CALL SXXXXX','CALL PROP_GHBARX')
                    else:
                        call = call.replace('CALL ','CALL PROP_')
            else:
                call = call.replace('CALL ','CALL EXTERNAL_')

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