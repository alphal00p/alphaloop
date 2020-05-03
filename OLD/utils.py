#####################################################
#                                                   #
#  Source file of the pyNLoop MG5aMC plugin.   #
#  Use only with consent of its authors.            #
#                                                   #
#  author: Valentin Hirschi, Ben Ruijl              #
#                                                   #
#####################################################

import os
import logging
import sys
import importlib

import madgraph.various.misc as misc
import madgraph.core.base_objects as base_objects

from madgraph import InvalidCmd, MadGraph5Error, MG5DIR

plugin_path = os.path.dirname(os.path.realpath( __file__ ))

pjoin = os.path.join

logger = logging.getLogger('pyNLoop.utils')

class utilsError(MadGraph5Error):
    """ Error for the pyNLoop plugin """
    pass

def find_one_loop_path(heptools_install_dir=None):
    """ Tries to automatically find the path of OneLOop."""

    # First look in HEPTools
    if not heptools_install_dir:
        heptools_install_dir = pjoin(MG5DIR, 'HEPTools')
    if os.path.exists(pjoin(heptools_install_dir, 'oneloop','libavh_olo.a')) and \
        os.path.exists(pjoin(heptools_install_dir, 'oneloop', 'avh_olo_dp.mod')):
        return pjoin(heptools_install_dir, 'oneloop'), pjoin(heptools_install_dir, 'oneloop')

    # If not there, then look in environment variables
    lib_path = misc.which_lib('libavh_olo.a')
    if lib_path is not None:
        lib_dir = os.path.dirname(lib_path)
        inc_dir = None
        if os.path.exists(pjoin(lib_dir, 'avh_olo_dp.mod')):
            inc_dir = lib_dir
        elif os.path.exists(pjoin(os.path.dirname(lib_dir),'inc','avh_olo_dp.mod')):
            inc_dir = pjoin(os.path.dirname(lib_dir),'inc')
        return lib_dir, inc_dir

    # Not found
    return None, None

class AVHOneLOopHook(object):

    def __init__(self, hook_path=None, OneLOop_lib_path=None, OneLOop_include_path=None, heptools_install_dir=None, f2py_compiler='f2py'):
        """ Given the path of both the hook interface and of the OneLOop library, this instantiates an object
        ready to compute arbitrary scalar loops."""

        if hook_path is None:
            hook_path = pjoin(plugin_path,'AVHOneLOopHook')

        if not os.path.exists(pjoin(hook_path,'AVHOneLOopHook.f90')):
            raise utilsError("The AVHOneLOopHook path specified '%s' does not appear valid."%hook_path)

        self.avh_oneloop_hook = None

        if any(path is None for path in [OneLOop_lib_path,OneLOop_include_path]):
            lib_path_found, inc_path_found = find_one_loop_path(heptools_install_dir=heptools_install_dir)
            if OneLOop_lib_path is None:
                OneLOop_lib_path = lib_path_found
            if OneLOop_include_path is None:
                OneLOop_include_path = inc_path_found
            if any(path is None for path in [OneLOop_lib_path, OneLOop_include_path]):
                logger.warning('The library OneLOop does not seem available on your system. Consider installing it'+
                               " using the command 'install oneloop'.")
                logger.warning(
                    'Access to OneLOop computations is therefore disabled (only necessary for some features of pyNLoop.')
                self.avh_oneloop_hook = None
                return

        if not os.path.exists(pjoin(OneLOop_lib_path, 'libavh_olo.a')):
            raise utilsError("The OneLOop library path specified '%s' does not appear valid." % OneLOop_lib_path)
        if not os.path.exists(pjoin(OneLOop_include_path, 'avh_olo_dp.mod')):
            raise utilsError("The OneLOop include path specified '%s' does not appear valid." % OneLOop_lib_path)

        # First make the shared object library with f2PY:
        try:
            logger.info('Now compiling shared library pyAVH_OneLOop_hook.so in path %s'%hook_path)
            misc.compile(arg=['pyAVH_OneLOop_hook.so',
                              'ONELOOP_LIB_PATH=%s'%OneLOop_lib_path,
                              'ONELOOP_INC_PATH=%s'%OneLOop_include_path,
                              'F2PY=%s'%f2py_compiler], cwd=hook_path)


            if os.path.dirname(hook_path) not in sys.path:
                sys.path.append(os.path.dirname(hook_path))
            self.avh_oneloop_hook = importlib.import_module('AVHOneLOopHook.pyAVH_OneLOop_hook')
        except Exception as e:
            logger.warning('The following exception occurred when trying to load the AVHOneLOopHook library: %s'%str(e))
            logger.warning('Access to OneLOop computations is therefore disabled (only necessary for some features of pyNLoop.')
            self.avh_oneloop_hook = None

    def is_available(self):
        """ Simply returns whether the library was correctly loaded."""

        return (self.avh_oneloop_hook is not None)

    def example_call(self):
        """ Just a dummy function to test the interface call to a one-loop box with hard-coded values."""

        if self.avh_oneloop_hook is None:
            logger.warning('The AVHOneLOopHook could not be properly loaded. Returning dummy zero result then.')
            return base_objects.EpsilonExpansion({  0 : 0.+0.j,
                                                   -1 : 0. + 0.j,
                                                   -2 : 0. + 0.j   })
        p_sq1 = 100.0**2
        p_sq2 = 200.0**2
        p_sq3 = 300.0**2
        p_sq4 = 400.0**2
        p_12  = 812438.8277451771e0
        p_23  = -406042.8596833937e0
        m1, m2 , m3 , m4 = 0., 0., 0., 0.

        result = self.avh_oneloop_hook.compute_one_loop_box(p_sq1,p_sq2,p_sq3,p_sq4,p_12,p_23,m1,m2,m3,m4)
        logger.info('The real part of the result should be: %s'%base_objects.EpsilonExpansion(
            {   0 :  -5.0195667135415558e-11,
               -1 : 0.,
               -2 : 0.   }))
        logger.info('The imaginary part of the result should be: %s'%base_objects.EpsilonExpansion(
            {   0 : -9.2157738446831588e-11,
               -1 : 0.,
               -2 : 0.   }))
        return base_objects.EpsilonExpansion({0: result[0].real, -1: result[1].real, -2: result[2].real}), \
               base_objects.EpsilonExpansion({0: result[0].imag, -1: result[1].imag, -2: result[2].imag})

    def compute_one_loop_box(self, PS_point, loop_propagator_masses=(0.,0.,0.,0.)):
        """ Computes an arbitrary box from the PS_point and loop_propagator_masses specified.
        Returns two epsilon expansion, one for the real part and one for the imaginary part."""

        if self.avh_oneloop_hook is None:
            logger.warning('The AVHOneLOopHook could not be properly loaded. Returning dummy zero result then.')
            return base_objects.EpsilonExpansion({0: 0., -1: 0., -2: 0.}), \
                   base_objects.EpsilonExpansion({0: 0., -1: 0., -2: 0.})

        p_sq1 = PS_point[1].square()
        p_sq2 = PS_point[2].square()
        p_sq3 = PS_point[3].square()
        p_sq4 = PS_point[4].square()
        p_12  = (PS_point[1] + PS_point[2]).square()
        p_23  = (PS_point[2] + PS_point[3]).square()
        m1, m2 , m3 , m4 = loop_propagator_masses

        result = self.avh_oneloop_hook.compute_one_loop_box(p_sq1,p_sq2,p_sq3,p_sq4,p_12,p_23,m1,m2,m3,m4)

        return base_objects.EpsilonExpansion({0: float(result[0].real), -1: float(result[1].real), -2: float(result[2].real)}), \
               base_objects.EpsilonExpansion({0: float(result[0].imag), -1: float(result[1].imag), -2: float(result[2].imag)})

    def compute_one_loop_triangle(self, PS_point, loop_propagator_masses=(0.,0.,0.)):
        """ Computes an arbitrary triangle from the PS_point and loop_propagator_masses specified.
        Returns two epsilon expansion, one for the real part and one for the imaginary part."""

        if self.avh_oneloop_hook is None:
            logger.warning('The AVHOneLOopHook could not be properly loaded. Returning dummy zero result then.')
            return base_objects.EpsilonExpansion({0: 0., -1: 0., -2: 0.}), \
                   base_objects.EpsilonExpansion({0: 0., -1: 0., -2: 0.})

        p_sq1 = PS_point[1].square()
        p_sq2 = PS_point[2].square()
        p_sq3 = PS_point[3].square()
        m1, m2 , m3 = loop_propagator_masses

        result = self.avh_oneloop_hook.compute_one_loop_triangle(p_sq1,p_sq2,p_sq3,m1,m2,m3)

        return base_objects.EpsilonExpansion({0: float(result[0].real), -1: float(result[1].real), -2: float(result[2].real)}), \
               base_objects.EpsilonExpansion({0: float(result[0].imag), -1: float(result[1].imag), -2: float(result[2].imag)})

    def compute_one_loop_bubble(self, PS_point, loop_propagator_masses=(0.,0.)):
        raise NotImplementedError('Function compute_one_loop_bubble not implemented yet.')

    def compute_one_loop_tadpole(self, PS_point, loop_propagator_masses=(0.,)):
        raise NotImplementedError('Function compute_one_loop_tadpole not implemented yet.')

