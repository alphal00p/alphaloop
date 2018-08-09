import ctypes
import numpy as np

_deformer = ctypes.CDLL("./DCD_interface.so")
#append Q
_deformer.append_Q.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int);
_deformer.append_Q.restype  = (ctypes.c_int);
#feed P+ and P-
_deformer.set_Pp.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int);
_deformer.set_Pm.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int);
_deformer.set_Pp.restype  = (ctypes.c_int);
_deformer.set_Pm.restype  = (ctypes.c_int);
#init class
_deformer.init.argtypes = ();
_deformer.init.restype  = (ctypes.c_int);
#delete previous defined class
_deformer.clear.argtypes = ();
_deformer.clear.restype  = (ctypes.c_void_p);
#get the array 
_deformer.get_deformed_loop_momentum.argtypes = ()
_deformer.get_deformed_loop_momentum.restype = (ctypes.POINTER(ctypes.c_double))
#get jacobian
_deformer.get_jacobian.argtypes = ()
_deformer.get_jacobian.restype  = (ctypes.c_double)
#get the array 
_deformer.deform_loop_momentum.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.c_int)
_deformer.deform_loop_momentum.restype  = (ctypes.c_int)


def append_Q(q):
    dim = len(q)
    array_type = ctypes.c_double * dim
    return _deformer.append_Q(array_type(*q),dim)

def set_Ppm(P_plus,P_minus):
    array_type = ctypes.c_double * 4
    return_val = []
    #set P+
    dim = len(P_plus)
    return_val += [_deformer.set_Pp(array_type(*P_plus),dim)]
    #set P-
    dim = len(P_minus)
    return_val += [_deformer.set_Pm(array_type(*P_minus),dim)]
    
    return max(return_val)
    
    
def deform_loop_momentum(k):
    dim = len(k)
    array_type = ctypes.c_double * dim
    return _deformer.deform_loop_momentum(array_type(*k),dim)

def get_deformed_loop_momentum():
    array_pointer = ctypes.cast(_deformer.get_deformed_loop_momentum(),ctypes.POINTER(ctypes.c_double * 8))
    array = np.frombuffer(array_pointer.contents)
    return [x + y*1j  for x, y in zip(array[:4],array[4:])]
def get_jacobian():
    return _deformer.get_jacobian()        

def init():
    return _deformer.init()
def clear():
    return _deformer.clear()


if __name__ == "__main__":
    #INPUTS 
    qs=[]
    qs+=[[-0.5, -0.5, 0.0, 0.0]]
    qs+=[[-1.0, 0.0, 0.0, 0.0]]
    qs+=[[-0.5, 0.0,-0.5, 0.0]]
    qs+=[[ 0.0, 0.0, 0.0, 0.0]]
    
    P_plus = [-1.0, 0.0, 0.0, 0.0]
    P_minus = [0.0, 0.0, 0.0, 0.0]
    
    loop_mom=[ 0.1, 0.2, 10.0, 0.4]
    
    #SET UP DEFORMER
    for q in qs:
        append_Q(q)
    set_Ppm(P_plus,P_minus)
    init()
    
    #USE IT FOR AS MANY LOOP MOMENTA AS YOU WISH
    print("loop_mom: %s"%str(loop_mom))
    deform_loop_momentum(loop_mom)
    print("deformed: %s"%str(get_deformed_loop_momentum()))
    print("jacobian: %s"%str(get_jacobian()))
    
    #RESET VARIABLES FOR NEW Qs
    clear()    
    