import logging
logger = logging.getLogger('alphaLoop.EsurfaceIntersectionFinder')
import math
import random

try:
    import cvxpy
    import scipy.optimize as optimize
except:
    logger.critical("Error: could not import package cvxpy and scipy.optimize necessary for processing the E-surfaces for the IR profile. Make sure it is installed.")

from pprint import pprint, pformat

class EsurfaceIntersectionFinder(object):

    def __init__(self, E_surfaces, cvxpy_coordinates, E_cm, seed_point_shifts=None, debug=False, **opts):
        
        self.E_surfaces = E_surfaces
        self.cvxpy_coordinates = cvxpy_coordinates
        self.debug = debug
        self.E_cm = E_cm
        if seed_point_shifts is None or isinstance(seed_point_shifts, int):
            n_shifts = seed_point_shifts if isinstance(seed_point_shifts, int) else 10
            self.seed_point_shifts = [
                [
                    [ random.random()*(self.E_cm/100.) for i_comp in range(3) ]
                    for i_loop in range(len(self.cvxpy_coordinates))
                ] for i_shift in range(0,n_shifts)
            ]
        else:
            self.seed_point_shifts = seed_point_shifts
        self.scipy_method = 'hybr'#'hybr'
    
    def find_intersection(self):

        if self.debug:
            logger.info("EsurfaceIntersectionFinder will now attempt to find the intersection of the following E_surfaces:\n%s"%(
                '\n'.join(pformat(E_surf) for E_surf in self.E_surfaces)
            ))

        # First test with cvxpy if we can find a point in the interior of all specified E-surfaces.
        # This is a prerequisite for the intersection to exist (when including the boudary as part of the interior) and 
        # it is useful to server as a powerful see for the scipy optimize function
        p = cvxpy.Problem(cvxpy.Minimize(1), [ E_surf['cxpy_expression'] <= 0 for E_surf in self.E_surfaces ] )

        seed_point = None
        try:
            p.solve()
        except cvxpy.SolverError:
            p.solve(solver=cvxpy.SCS, eps=1e-9)
        except Exception as e:
            pass

        if self.debug:
            logger.info("The attempt at finding a seed in the interior of all surfaces with cvxpy yielded:\n%s"%str(p))

        if p.status == cvxpy.OPTIMAL:
            seed_point = [ [ float(c.value[0]), float(c.value[1]), float(c.value[2]) ]  for c in self.cvxpy_coordinates ]
        
        if seed_point is None:
            # No possible intersection
            if self.debug: logger.info("Could not find a seed point in the interior of all E-surfaces. There is therefore no intersection.")
            return None


        if self.debug: logger.info("The E-surface equations evaluate as follows for the seed point found x=%s\n%s"%(
            seed_point,
            '\n'.join('E-surface with id %d : %.16e'%(
                E_surf['id'], self.E_surface( seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
            )
        ))

        for i_shift, seed_point_shift in enumerate(self.seed_point_shifts):

            shifted_seed = [ [ seed_comp+seed_point_shift[i_vec][i_comp] for i_comp, seed_comp in enumerate(seed_vec) ] for i_vec, seed_vec in enumerate(seed_point) ]
            scipy_seed = sum(shifted_seed,[])

            if self.debug: logger.info("An intersection is possible, scipy.optimize will now attempt to find an intersection from the seed point with shift #%d: %s."%(i_shift,str(scipy_seed)))

            #pprint(self.intersection_function(scipy_seed))
            #pprint(self.intersection_function_jac(scipy_seed))
            scipy_res = optimize.root(
                self.intersection_function,
                scipy_seed,
                jac=self.intersection_function_jac, 
                method=self.scipy_method,
                options={
                    'xtol'   : 1.49012e-08,
                    'maxfev' : 1000*(len(scipy_seed)+1)
                }
            )

            if self.debug:
                logger.info("The attempt at finding an intersection with scipy.optimize for seed shift #%d yielded:\n%s"%(i_shift,str(scipy_res)))

            #logger.debug(str(scipy_res))
            if not scipy_res.success:
                if self.debug: logger.info("Could not find an intersection with scipy.optimize for seed shift #%d"%i_shift)
                continue
            else:
                #pprint(self.intersection_function(scipy_res.x))
                #pprint(self.intersection_function_jac(scipy_res.x))
                intersection_point = [ scipy_res.x[i:i+3] for i in range(0,len(scipy_res.x),3) ]
                if self.debug: logger.info("The E-surface equations evaluate as follows for the intersection point found using shift #%d x=%s\n%s"%(
                    i_shift,
                    intersection_point,
                    '\n'.join('E-surface with id %d : %.16e'%(
                        E_surf['id'], self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] ) ) for E_surf in self.E_surfaces
                    )
                ))
                return intersection_point

        if self.debug: logger.info("Could not find an intersection with scipy.optimize with %d attempts of shifting the seed point."%len(self.seed_point_shifts))

        return None

    def delta(self, loop_momenta, loop_sig, shift, m_squared):
        
        k = [ sum( l[i]*factor for l, factor in zip(loop_momenta,loop_sig)) for i in range(0,3) ]

        return math.sqrt(
             (k[0]**2+k[1]**2+k[2]**2)+
            +2*(k[0]*shift[0]+k[1]*shift[1]+k[2]*shift[2])
            +shift[0]**2+shift[1]**2+shift[2]**2
            +m_squared
        )
    
    def ddelta(self, loop_momenta, loop_sig, shift, m_squared, loop_index, component_index):

        k_derived = loop_momenta[loop_index][component_index]*loop_sig[loop_index]

        if k_derived == 0.:
            return 0.

        k_not_derived = sum( l[component_index]*factor for i_loop_mom, (l, factor) in enumerate(zip(loop_momenta,loop_sig)) if i_loop_mom!=loop_index )

        return (k_derived + k_not_derived + shift[component_index])/self.delta(loop_momenta, loop_sig, shift, m_squared)

    def E_surface(self, loop_momenta, onshell_propagators, E_surface_shift):

        return sum( self.delta(loop_momenta,osp['loop_sig'],osp['v_shift'],osp['m_squared']) for osp in onshell_propagators) + E_surface_shift

    def dE_surface(self, loop_momenta, onshell_propagators, E_surface_shift, loop_index, component_index):

        return sum( self.ddelta(loop_momenta,osp['loop_sig'],osp['v_shift'],osp['m_squared'],loop_index, component_index) for osp in onshell_propagators)

    def intersection_function(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]

        res = []
        for E_surf in self.E_surfaces:
            res.append( self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift']) )

        # Pad remaining entries with no constraints if there are some left.
        res += [0.]*(len(xs)-len(res))

        return res

    def intersection_function_jac(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]

        res = []
        for E_surf in self.E_surfaces:
            res.append([])
            for loop_index in range(0,len(loop_momenta)):
                for component_index in range(0,3):
                    res[-1].append( self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index, component_index) )

        # Pad remaining entries with no constraints if there are some left.
        res += [ [0.]*(len(loop_momenta)*3) ]*(len(xs)-len(res))

        return res