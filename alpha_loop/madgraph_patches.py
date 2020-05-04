import alpha_loop.utils as utils
import alpha_loop.interface as interface

# Patch the class method `generate_matrix_elements` of helas_objects
# so as to force the wavefunction recycling option to be disabled.
# Ideally this should be done more cleanly via an actual option of the PLUGIN.
# This is sadly not (yet?) exposed to PLUGINS in MG5aMC.
interface.logger.warning("%salphaLoop temporarily dynamically patched 'helas_objects.HelasMultiProcess.generate_matrix_elements' so as to turn off wavefunction recycling.%s"%(
utils.bcolors.RED, utils.bcolors.ENDC))
import madgraph.core.helas_objects as helas_objects
orginal_generate_matrix_elements=helas_objects.HelasMultiProcess.generate_matrix_elements
def patched_generate_matrix_elements(cls, *args, matrix_element_opts={},**opts):
    matrix_element_opts['optimization']=0
    return orginal_generate_matrix_elements(cls,*args, matrix_element_opts=matrix_element_opts,**opts)
helas_objects.HelasMultiProcess.generate_matrix_elements = patched_generate_matrix_elements

interface.logger.warning("%salphaLoop temporarily dynamically patched 'aloha_writers.ALOHAWriterForFortran.define_expression' so as to turn off all propagator denominators.%s"%(
utils.bcolors.RED, utils.bcolors.ENDC))
# Patch the class method `aloha_writers.ALOHAWriterForFortran.define_expression` so as
# to disable all propagator denominators.
# TODO: Ideally find a less brutal way of achieving this, probably by finding 
# how to assign a custom denominator to all routines.
import aloha
import aloha.aloha_object as aloha_object
import aloha.aloha_lib as aloha_lib
import aloha.aloha_writers as aloha_writers
from six import StringIO
import madgraph.various.misc as misc
def patched_define_expression(self):
    """Define the functions in a 100% way """
    out = StringIO()
    if self.routine.contracted:
        all_keys = list(self.routine.contracted.keys())
        all_keys.sort()
        for name in all_keys:
            obj = self.routine.contracted[name]
            out.write(' %s = %s\n' % (name, self.write_obj(obj)))
            self.declaration.add(('complex', name))
            
    
    def sort_fct(a, b):
        if len(a) < len(b):
            return -1
        elif len(a) > len(b):
            return 1
        elif a < b:
            return -1
        else:
            return +1
        
    keys = list(self.routine.fct.keys())        
    keys.sort(key=misc.cmp_to_key(sort_fct))
    for name in keys:
        fct, objs = self.routine.fct[name]
        format = ' %s = %s\n' % (name, self.get_fct_format(fct))
        try:
            text = format % ','.join([self.write_obj(obj) for obj in objs])
        except TypeError:
            text = format % tuple([self.write_obj(obj) for obj in objs])
        finally:
            out.write(text)
    
    numerator = self.routine.expr
    if not 'Coup(1)' in self.routine.infostr:
        coup_name = 'COUP'
    else:
        coup_name = '%s' % self.change_number_format(1)
    
    
    if not self.offshell:
        if coup_name == 'COUP':
            formatted = self.write_obj(numerator.get_rep([0]))
            if formatted.startswith(('+','-')):
                out.write(' vertex = COUP*(%s)\n' % formatted)
            else:
                out.write(' vertex = COUP*%s\n' % formatted)
        else:
            out.write(' vertex = %s\n' % self.write_obj(numerator.get_rep([0])))
    else:
        OffShellParticle = '%s%d' % (self.particles[self.offshell-1],\
                                                              self.offshell)
        is_loop = False
        if 'L' in self.tag:
            if self.tag.count('L') == 1 and 'PL' in self.tag:
                is_loop = False
            else:
                is_loop = True
                
        if not is_loop:
            coeff = 'denom*'    
            if not aloha.complex_mass:
                if self.routine.denominator:
                    out.write('    denom = %(COUP)s/(%(denom)s)\n' % {'COUP': coup_name,\
                            'denom':self.write_obj(self.routine.denominator)}) 
                else:
                    # Below is the hack to de-activate the denominators.
                    # TODO: find a more elegant way to do this, probably by finding how to assign a custom denominator to all routines.
                    #out.write('    denom = %(COUP)s/VVVV(P%(i)s(0)**2-P%(i)s(1)**2-P%(i)s(2)**2-P%(i)s(3)**2 - M%(i)s * (M%(i)s -CI* W%(i)s))\n' % \
                    #          {'i': self.outgoing, 'COUP': coup_name})
                    out.write('    denom = %(COUP)s\n' % {'COUP': coup_name})
            else:
                if self.routine.denominator:
                    raise Exception('modify denominator are not compatible with complex mass scheme')                
                out.write('    denom = %(COUP)s/(P%(i)s(0)**2-P%(i)s(1)**2-P%(i)s(2)**2-P%(i)s(3)**2 - M%(i)s**2)\n' % \
                  {'i': self.outgoing, 'COUP': coup_name})
            self.declaration.add(('complex','denom'))
            if aloha.loop_mode:
                ptype = 'list_complex'
            else:
                ptype = 'list_double'
            self.declaration.add((ptype,'P%s' % self.outgoing))
        else:
            if coup_name == 'COUP':
                coeff = 'COUP*'
            else:
                coeff = ''
        to_order = {}  
        for ind in numerator.listindices():
            formatted = self.write_obj(numerator.get_rep(ind))
            if formatted.startswith(('+','-')):
                if '*' in formatted:
                    formatted = '(%s)*%s' % tuple(formatted.split('*',1))
                else:
                    if formatted.startswith('+'):
                        formatted = formatted[1:]
                    else:
                        formatted = '(-1)*%s' % formatted[1:]
            to_order[self.pass_to_HELAS(ind)] = \
                    '    %s(%d)= %s%s\n' % (self.outname, self.pass_to_HELAS(ind)+1, 
                    coeff, formatted)
        key = list(to_order.keys())
        key.sort()
        for i in key:
            out.write(to_order[i])
    return out.getvalue()

aloha_writers.ALOHAWriterForFortran.define_expression = patched_define_expression
