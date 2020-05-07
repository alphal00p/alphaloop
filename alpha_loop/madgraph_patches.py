import alpha_loop.utils as utils
import alpha_loop.interface as interface
import madgraph.various.misc as misc

import madgraph.core.helas_objects as helas_objects
original_generate_matrix_elements = helas_objects.HelasMultiProcess.generate_matrix_elements
def remove_wavefunction_recycling(restore=False):
    if restore:
        interface.logger.info("%sPatch of 'helas_objects.HelasMultiProcess.generate_matrix_elements' reverted.%s"%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        helas_objects.HelasMultiProcess.generate_matrix_elements = original_generate_matrix_elements
        return

    # Patch the class method `generate_matrix_elements` of helas_objects
    # so as to force the wavefunction recycling option to be disabled.
    # Ideally this should be done more cleanly via an actual option of the PLUGIN.
    # This is sadly not (yet?) exposed to PLUGINS in MG5aMC.
    interface.logger.warning("%salphaLoop temporarily dynamically patched 'helas_objects.HelasMultiProcess.generate_matrix_elements' so as to turn off wavefunction recycling.%s"%(
    utils.bcolors.RED, utils.bcolors.ENDC))
    orginal_generate_matrix_elements=helas_objects.HelasMultiProcess.generate_matrix_elements
    def patched_generate_matrix_elements(cls, *args, matrix_element_opts={},**opts):
        # MODIF IS THE LINE BELOW
        matrix_element_opts['optimization']=0
        return orginal_generate_matrix_elements(cls,*args, matrix_element_opts=matrix_element_opts,**opts)
    helas_objects.HelasMultiProcess.generate_matrix_elements = patched_generate_matrix_elements

import aloha.aloha_writers as aloha_writers 
original_define_expression = aloha_writers.ALOHAWriterForFortran.define_expression
def remove_propagator_denominators(restore=False):
    if restore:
        interface.logger.info("%sPatch of 'aloha_writers.ALOHAWriterForFortran.define_expression' reverted.%s"%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        aloha_writers.ALOHAWriterForFortran.define_expression = original_define_expression
        return

    interface.logger.warning("%salphaLoop temporarily dynamically patched 'aloha_writers.ALOHAWriterForFortran.define_expression' so as to turn off all propagator denominators.%s"%(
    utils.bcolors.RED, utils.bcolors.ENDC))
    # Patch the class method `aloha_writers.ALOHAWriterForFortran.define_expression` so as
    # to disable all propagator denominators.
    # TODO: Ideally find a less brutal way of achieving this, probably by finding 
    # how to assign a custom denominator to all routines.
    import aloha
    import aloha.aloha_object as aloha_object
    import aloha.aloha_lib as aloha_lib
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
                        # MODIF IS THE LINE BELOW
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

import aloha.create_aloha as create_aloha
original_compute_subset = create_aloha.AbstractALOHAModel.compute_subset
def force_aloha_to_loop_mode(restore=False):
    if restore:
        interface.logger.info("%sPatch of 'create_aloha.AbstractALOHAModel.compute_subset' reverted.%s"%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        create_aloha.AbstractALOHAModel.compute_subset = original_compute_subset
        return
    interface.logger.warning("%salphaLoop temporarily dynamically patched 'create_aloha.AbstractALOHAModel.compute_subset' so as to force it to use the loop mode in all cases.%s"%(
    utils.bcolors.RED, utils.bcolors.ENDC))
    # Force aloha to store momenta in the first four and not two components of HELAS wavefunctions
    # For this we must force aloha.loop_mode to be set to True even in the absence of loops
    import aloha
    import aloha.create_aloha as create_aloha
    from aloha.create_aloha import AbstractRoutineBuilder, CombineRoutineBuilder
    import logging
    import time
    import operator
    logger = logging.getLogger('PATCHED_ALOHA')
    def patched_compute_subset(self, data):
        """ create the requested ALOHA routine. 
        data should be a list of tuple (lorentz, tag, outgoing)
        tag should be the list of special tag (like conjugation on pair)
        to apply on the object """
        logger.info('aloha starts to compute helicity amplitudes')
        start = time.time()
        # Search identical particles in the vertices in order to avoid
        #to compute identical contribution
        self.look_for_symmetries()
        # reorganize the data (in order to use optimization for a given lorentz
        #structure
        # MODIF BELOW:
        #aloha.loop_mode = False 
        aloha.loop_mode = True
        # self.explicit_combine = False
        request = {}
        for list_l_name, tag, outgoing in data:
            #allow tag to have integer for retro-compatibility
            all_tag = tag[:]
            conjugate = [i for i in tag if isinstance(i, int)]
            
            tag =  [i for i in tag if isinstance(i, str) and not i.startswith('P')]
            tag = tag + ['C%s'%i for i in conjugate]             
            tag = tag + [i for i in all_tag if isinstance(i, str) and  i.startswith('P')] 
            
            conjugate = tuple([int(float(c[1:])) for c in tag if c.startswith('C')])
            loop = any((t.startswith('L') for t in tag))
            if loop:
                aloha.loop_mode = True
                self.explicit_combine = True
                    
            for l_name in list_l_name:
                try:
                    request[l_name][conjugate].append((outgoing,tag))
                except Exception:
                    try:
                        request[l_name][conjugate] = [(outgoing,tag)]
                    except Exception:
                        request[l_name] = {conjugate: [(outgoing,tag)]}
                        
        # Loop on the structure to build exactly what is request
        for l_name in request:
            lorentz = eval('self.model.lorentz.%s' % l_name)
            if lorentz.structure == 'external':
                for tmp in request[l_name]:
                    for outgoing, tag in request[l_name][tmp]:
                        name = aloha_writers.get_routine_name(lorentz.name,outgoing=outgoing,tag=tag)
                        if name not in self.external_routines:
                            self.external_routines.append(name)
                continue
            
            builder = AbstractRoutineBuilder(lorentz, self.model)
            
            for conjg in request[l_name]:
                #ensure that routines are in rising order (for symetries)
                def sorting(a,b):
                    if a[0] < b[0]: return -1
                    else: return 1
                routines = request[l_name][conjg]
                routines.sort(key=misc.cmp_to_key(sorting))
                if not conjg:
                    # No need to conjugate -> compute directly
                    self.compute_aloha(builder, routines=routines)
                else:
                    # Define the high level conjugate routine
                    conjg_builder = builder.define_conjugate_builder(conjg)
                    # Compute routines
                    self.compute_aloha(conjg_builder, symmetry=lorentz.name,
                                        routines=routines)
            
        
        # Build mutiple lorentz call
        for list_l_name, tag, outgoing in data:
            if len(list_l_name) ==1:
                continue
            #allow tag to have integer for retrocompatibility
            conjugate = [i for i in tag if isinstance(i, int)]
            all_tag = tag[:]
            tag =  [i for i in tag if isinstance(i, str) and not i.startswith('P')]
            tag = tag + ['C%s'%i for i in conjugate] 
            tag = tag + [i for i in all_tag if isinstance(i, str) and  i.startswith('P')] 
            
            if not self.explicit_combine:
                lorentzname = list_l_name[0]
                lorentzname += ''.join(tag)
                if (lorentzname, outgoing) in self:
                    self[(lorentzname, outgoing)].add_combine(list_l_name[1:])
                else:
                    lorentz = eval('self.model.lorentz.%s' % list_l_name[0])
                    assert lorentz.structure == 'external'
            else:
                l_lorentz = []
                for l_name in list_l_name: 
                    l_lorentz.append(eval('self.model.lorentz.%s' % l_name))
                builder = CombineRoutineBuilder(l_lorentz)
                            
                for conjg in request[list_l_name[0]]:
                    #ensure that routines are in rising order (for symetries)
                    def sorting(a,b):
                        if a[0] < b[0]: return -1
                        else: return 1
                    routines = request[list_l_name[0]][conjg]
                    routines.sort(key=operator.itemgetter(0))
                    if not conjg:
                        # No need to conjugate -> compute directly
                        self.compute_aloha(builder, routines=routines)
                    else:
                        # Define the high level conjugate routine
                        conjg_builder = builder.define_conjugate_builder(conjg)
                        # Compute routines
                        self.compute_aloha(conjg_builder, symmetry=lorentz.name,
                                        routines=routines)
        
        logger.info("aloha creates %s routines in  %0.3f s", AbstractRoutineBuilder.counter, time.time()-start)

    create_aloha.AbstractALOHAModel.compute_subset = patched_compute_subset


import aloha.aloha_writers as aloha_writers
original_get_header_txt = aloha_writers.ALOHAWriterForFortran.get_header_txt
original_get_declaration_txt = aloha_writers.ALOHAWriterForFortran.get_declaration_txt

def enable_manual_complex_conjugation_in_aloha(restore=False):
    if restore:
        interface.logger.info("%sPatch of 'aloha_writers.ALOHAWriterForFortran.get_header_txt' and 'aloha_writers.ALOHAWriterForFortran.get_declaration_txt' reverted.%s"%(utils.bcolors.GREEN, utils.bcolors.ENDC))
        aloha_writers.ALOHAWriterForFortran.get_header_txt = original_get_header_txt
        aloha_writers.ALOHAWriterForFortran.get_declaration_txt = original_get_declaration_txt
        return
    text = "alphaLoop temporarily dynamically patched 'aloha_writers.ALOHAWriterForFortran.get_header_txt' and "+\
        "'aloha_writers.ALOHAWriterForFortran.get_declaration_txt' so as to allow for 'manual' complex conjugation in ALOHA."
    interface.logger.warning("%s%s%s"%(utils.bcolors.RED, text, utils.bcolors.ENDC))

    # Force aloha to store momenta in the first four and not two components of HELAS wavefunctions
    # For this we must force aloha.loop_mode to be set to True even in the absence of loops
    import aloha
    import aloha.create_aloha as create_aloha
    import aloha.aloha_lib as aloha_lib
    KERNEL = aloha_lib.KERNEL
    from aloha.create_aloha import AbstractRoutineBuilder, CombineRoutineBuilder
    import logging
    import operator
    import cmath
    from six import StringIO

    def patched_get_header_txt(self, name=None, couplings=None, **opt):
        """Define the Header of the fortran file. 
        """
        if name is None:
            name = self.name
           
        out = StringIO()
        # define the type of function and argument
        
        arguments = [arg for format, arg in self.define_argument_list(couplings)]
        if not self.offshell:
            output = 'vertex'
            self.declaration.add(('complex','vertex'))
        else:
            output = '%(spin)s%(id)d' % {
                     'spin': self.particles[self.outgoing -1],
                     'id': self.outgoing}
            self.declaration.add(('list_complex', output))
       
#       ======================
#       START MODIFS for LTD^2
#       ======================
        new_arguments = []
        for arg in arguments:
            if arg.upper().startswith('COUP'):
                new_arguments.append(arg+'_IN')
            else:
                new_arguments.append(arg)
        arguments=new_arguments

        out.write('subroutine %(name)s(%(args)s,leftright,%(output)s)\n' % \
                  {'output':output, 'name': name, 'args': ', '.join(arguments)})

#        out.write('subroutine %(name)s(%(args)s,%(output)s)\n' % \
#                  {'output':output, 'name': name, 'args': ', '.join(arguments)})
        
#       ======================
#       END MODIFS for LTD^2
#       ======================

        return out.getvalue() 

    def patched_get_declaration_txt(self):
        """ Prototype for how to write the declaration of variable
            Include the symmetry line (entry FFV_2)
        """
        
        out = StringIO()
        out.write('implicit none\n')
        # Check if we are in formfactor mode
        if self.has_model_parameter:
            out.write(' include "../MODEL/input.inc"\n')
            out.write(' include "../MODEL/coupl.inc"\n')
        argument_var = [name for type,name in self.call_arg]

#       ======================
#       START MODIFS for LTD^2
#       ======================

        # define the complex number CI = 0+1j
        if 'MP' in self.tag:
            out.write(' complex*32 CI, CI_IN\n')
            if KERNEL.has_pi:
                out.write(' REAL ( KIND = 16 ) PI\n')
        else:
            out.write(' complex*16 CI, CI_IN\n')
            if KERNEL.has_pi:
                out.write(' double precision PI\n')

        out.write(' integer leftright\n')
        out.write(' parameter (CI_IN=(%s,%s))\n' % 
                    (self.change_number_format(0),self.change_number_format(1)))

#       ======================
#       END MODIFS for LTD^2
#       ======================

        if KERNEL.has_pi:
            out.write(' parameter (PI=%s)\n' % self.change_number_format(cmath.pi))
        
        
        for type, name in self.declaration.tolist():
            if type.startswith('list'):
                type = type[5:]
                #determine the size of the list
                if name in argument_var:
                    size ='*'
                elif name.startswith('P'):
                    size='0:3'
                elif name[0] in ['F','V']:
                    if aloha.loop_mode:
                        size = 8
                    else:
                        size = 6
                elif name[0] == 'S':
                    if aloha.loop_mode:
                        size = 5
                    else:
                        size = 3
                elif name[0] in ['R','T']: 
                    if aloha.loop_mode:
                        size = 20
                    else:
                        size = 18
                else:
                    size = '*'
    
                out.write(' %s %s(%s)\n' % (self.type2def[type], name, size))
            elif type == 'fct':
                if name.upper() in ['EXP','LOG','SIN','COS','ASIN','ACOS']:
                    continue
                out.write(' %s %s\n' % (self.type2def['complex'], name))
                out.write(' external %s\n' % (name))
            else:
                out.write(' %s %s\n' % (self.type2def[type], name))
#       ======================
#       START MODIFS for LTD^2
#       ======================
                if name.startswith('COUP'):
                    out.write(' %s %s_IN\n' % (self.type2def[type], name))
#       ======================
#       END MODIFS for LTD^2
#       ======================

        # Add the lines corresponding to the symmetry
        
        #number = self.offshell
        #arguments = [name for format, name in self.define_argument_list()]
        #new_name = self.name.rsplit('_')[0] + '_%s' % new_nb
        #return '%s\n    call %s(%s)' % \
        #    (self.get_header_txt(new_name, couplings), self.name, ','.join(arguments))
        couplings = [name for type, name in self.declaration if name.startswith('COUP') ]
        couplings.sort()
        for elem in self.routine.symmetries:
            new_name = self.name.rsplit('_',1)[0] + '_%s' % elem
            out.write('%s\n' % self.get_header_txt(new_name, couplings).replace('subroutine','entry'))


#       ======================
#       START MODIFS for LTD^2
#       ======================
        
        conjg_code = ["IF (leftright.eq.1) THEN"]
        conjg_code.append("CI = CI_IN")
        for vtype, name in self.declaration.tolist():
            if name.startswith('COUP'):
                conjg_code.append("%s = %s_IN"%(name,name))
        conjg_code.append("ELSE")
        conjg_code.append("CI = -CI_IN")
        for vtype, name in self.declaration.tolist():
            if name.startswith('COUP'):
                if 'MP' in self.tag:
                    conjg_code.append("%s = CONJG(%s_IN)"%(name,name))
                else:
                    conjg_code.append("%s = DCONJG(%s_IN)"%(name,name))
        conjg_code.append("ENDIF\n")
        out.write('\n'.join(conjg_code))

#       ======================
#       END MODIFS for LTD^2
#       ======================

        return out.getvalue()

    aloha_writers.ALOHAWriterForFortran.get_header_txt = patched_get_header_txt
    aloha_writers.ALOHAWriterForFortran.get_declaration_txt = patched_get_declaration_txt