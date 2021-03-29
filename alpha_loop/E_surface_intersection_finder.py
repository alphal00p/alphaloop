import logging
logger = logging.getLogger('alphaLoop.EsurfaceIntersectionFinder')
import math
import random
import warnings
import numpy
from LTD.vectors import Vector
from pprint import pformat
import itertools

try:
    import cvxpy
    import scipy.optimize as optimize
except:
    logger.critical("Error: could not import package cvxpy and scipy.optimize necessary for processing the E-surfaces for the IR profile. Make sure it is installed.")

from pprint import pprint, pformat

class EsurfaceIntersectionFinder(object):

    def __init__(self, E_surfaces, cvxpy_coordinates, E_cm, seed_point_shifts=None, debug=False, frozen_momenta=None, **opts):
        
        self.E_surfaces = E_surfaces

        self.cvxpy_coordinates = cvxpy_coordinates

        self.debug = debug
        self.E_cm = E_cm
        self.consistency_check_threshold = 1.0e-6
        self.scipy_tolerance = 1.49012e-08

        # Ignore the spatial component of the frozen momenta
        if frozen_momenta is not None:
            self.frozen_momenta = {}
            self.frozen_momenta['in'] = [list(v[1:]) for v in frozen_momenta['in']]
            self.frozen_momenta['out'] = [list(v[1:]) for v in frozen_momenta['out']]
        else:
            self.frozen_momenta = None

        self.n_frozen_momenta = 0 if self.frozen_momenta is None else len(self.frozen_momenta['out'])

        # Assign a truly random seed to cvxpy problem solving
        for cvxpy_coord in self.cvxpy_coordinates:
            cvxpy_coord.value = [ numpy.float64(random.random()*self.E_cm) for i in range(0,3) ]

        self.cvxpy_directions = []
        for i_lop in range(0,len(self.cvxpy_coordinates)):
            direction = []
            for i_comp in range(0,3):
                r = (random.random()-0.5)*2.
                direction.append( ((r-1.)*2.*self.E_cm if r<0. else (r+1.)*2.*self.E_cm) )
            self.cvxpy_directions.append(cvxpy.Constant(direction))
        
        if seed_point_shifts is None or isinstance(seed_point_shifts, int):
            n_shifts = seed_point_shifts if isinstance(seed_point_shifts, int) else 10
            self.seed_point_shifts = [
                [
                    [ random.random()*(self.E_cm/100.)*(10*((i_shift+1)/float(n_shifts//2))) for i_comp in range(3) ]
                    for i_loop in range(len(self.cvxpy_coordinates))
                ] for i_shift in range(0,n_shifts//2)
            ]
            self.seed_point_shifts += [
                [
                    [ random.random()*(self.E_cm/100.)*(200*((i_shift+1)/float(n_shifts-(n_shifts//2)))) for i_comp in range(3) ]
                    for i_loop in range(len(self.cvxpy_coordinates))
                ] for i_shift in range(0,n_shifts-(n_shifts//2))
            ]
        else:
            self.seed_point_shifts = seed_point_shifts

        # Always finish with the no-shift case
        self.seed_point_shifts.append( [ [0.,]*3 ]*len(self.cvxpy_coordinates) )
    
        self.scipy_method = 'hybr'

        self.use_scipy_root = True
        self.use_scipy_fsolve = False
        self.use_scipy_minimize = False

        self.maximal_cvxpy_seed_point_improvement_steps = 10
        self.maximal_cvxpy_find_intersection_steps = 30

    def project_onto_E_surface_closest(self, point, direction, orientations=('plus','minus')):

        #point = [ [pi+random.random()*1.0 for pi in p] for p in point ]

        max_norm = max( abs(n) for n in direction )

        start_point = [ Vector(v) for v in point ]

        if max_norm == 0.:
            return point

        normalized_direction = [normal/max_norm for normal in direction]

        if self.debug: 
            logger.info("Now projecting point x=\n%s\nusing normal=\n%s"%(str(start_point),direction))

        plus_point = [ v + 2.*self.E_cm*normalized_direction[i_loop] for i_loop, v in enumerate(start_point) ]
        minus_point = [ v - 2.*self.E_cm*normalized_direction[i_loop] for i_loop, v in enumerate(start_point) ]
        cvxpy_directions = {
            'plus' : {
                'direction' : plus_point,
                'starting_distance' : sum([ (s-v).square() for s,v in zip(plus_point, start_point) ])
            },
            'minus' : {
                'direction' : minus_point,
                'starting_distance' : sum([ (s-v).square() for s,v in zip(minus_point, start_point) ])
            }
        }
        cvxpy_directions = {k:v for k,v in cvxpy_directions.items() if k in orientations}

        for orientation in cvxpy_directions:
            objective = cvxpy.norm(cvxpy.hstack([cvxpy.Constant(direction)-self.cvxpy_coordinates[i_dir] for i_dir, direction in enumerate(cvxpy_directions[orientation]['direction']) ]), 2)

            p = cvxpy.Problem(
                cvxpy.Minimize( objective ), 
                [ E_surf['cxpy_expression'] <= 0 for E_surf in self.E_surfaces[:] ] )

            cvxpyres=p.solve()
            #print(cvxpyres, p.status)
            try:
                cvxpyres=p.solve()
            except cvxpy.SolverError:
                p.solve(solver=cvxpy.SCS, eps=1e-9)
            except Exception as e:
                return None
            
            cvxpy_directions[orientation]['solution'] = [ Vector([float(c.value[0]), float(c.value[1]), float(c.value[2])]) for c in self.cvxpy_coordinates ]
            cvxpy_directions[orientation]['distance'] = sum([ (s-v).square() for s,v in zip(cvxpy_directions[orientation]['solution'], start_point) ])
            cvxpy_directions[orientation]['relative_distance'] = cvxpy_directions[orientation]['distance']/cvxpy_directions[orientation]['starting_distance']
        
        all_res = sorted(cvxpy_directions.items(), key = lambda s: s[1]['distance'])

        if self.debug: 
            logger.info("Solutions obtained:\n%s\nuReturning:%s"%(pformat(cvxpy_directions),all_res[0][0]))

        return [ list(v) for v in all_res[0][1]['solution'] ]

    def improve_seed_point_directions_per_signature(self, seed_point, maximal_cvxpy_seed_point_improvement_steps, include_push_away_soft=False):

        solved = False
        n_steps = 0
        new_seed_point = seed_point

        dir_norm = 5*self.E_cm
        if not include_push_away_soft:
            base_directions = [
                [1.,0.,0.],[0.,1.,0.],[0.,0.,1.]
            ]
            extra_directions = [
                [-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.]
            ]
        else:
            # Pick the directions below to get less specific points
            base_directions = [
                [1.,0.1,0.2],[0.2,1.,0.1],[0.3,0.1,1.]
            ]
            extra_directions = [
                [-1.,0.1,0.2],[0.1,-1.,0.2],[0.1,0.2,-1.]
            ]

        all_loop_sigs = list(set([ osp['loop_sig'] for E_surf in self.E_surfaces for osp in E_surf['onshell_propagators'] ]))

        # Now add all possible combinations, from less likely to yield a solution to more likely
        all_possible_loop_signatures_directions = list(itertools.combinations_with_replacement(base_directions, len(all_loop_sigs)))
        if len(all_possible_loop_signatures_directions)<maximal_cvxpy_seed_point_improvement_steps:
            for additional_combinations in itertools.combinations_with_replacement(base_directions+extra_directions, len(all_loop_sigs)):
                if additional_combinations not in all_possible_loop_signatures_directions:
                    all_possible_loop_signatures_directions.append(additional_combinations)
        if len(all_possible_loop_signatures_directions)<maximal_cvxpy_seed_point_improvement_steps:
            for additional_combinations in itertools.product(*([base_directions,]*len(all_loop_sigs))):
                if additional_combinations not in all_possible_loop_signatures_directions:
                    all_possible_loop_signatures_directions.append(additional_combinations)
        if len(all_possible_loop_signatures_directions)<maximal_cvxpy_seed_point_improvement_steps:
            for additional_combinations in itertools.product(*([base_directions+extra_directions,]*len(all_loop_sigs))):
                if additional_combinations not in all_possible_loop_signatures_directions:
                    all_possible_loop_signatures_directions.append(additional_combinations)

        best_consistency_test = None
        best_seed_so_far = None
        for signatures_directions in all_possible_loop_signatures_directions:
            
            if solved or n_steps > maximal_cvxpy_seed_point_improvement_steps:
                break
            
            # scramble momenta that are not appearing in any of the E-surfaces selected
            for i_loop in range(0,len(self.cvxpy_coordinates)):
                if not any(sig[i_loop]!=0 for sig in all_loop_sigs):
                    new_seed_point[i_loop] = [random.random()*(self.E_cm/2.) for _ in range(0,3)]

            direction_for_signatures = {s: signatures_directions[i] for i, s in enumerate(all_loop_sigs) }

            solution_sum = sum(abs(self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
            consistency_test = abs(solution_sum/self.E_cm)

            if self.debug:
                logger.info("Step %d of improvement on the seed point\nx=%s\n that evaluates E-surface to:\n%s\nwith score %.4e by testing loop signature directions:\n%s"%(
                    n_steps+1,
                    str(new_seed_point),
                    '\n'.join('E-surface with id %d : %.16e'%(
                        E_surf['id'], self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                    ),
                    abs(solution_sum/self.E_cm),
                    str(direction_for_signatures)
                ))

            if consistency_test < self.consistency_check_threshold:
                if self.debug: logger.info('Current solution satisfactory enough, returning now.')
                solved = True
                return new_seed_point, solved
            
            if best_consistency_test is None or consistency_test<best_consistency_test:
                best_consistency_test = consistency_test
                best_seed_so_far = new_seed_point

            n_steps += 1
            objective = sum( cvxpy.norm(
                     cvxpy.sum([coord*wgt for coord, wgt in zip(self.cvxpy_coordinates, osp['loop_sig'])])- cvxpy.Constant(dir_norm)*cvxpy.Constant(direction_for_signatures[osp['loop_sig']]), 2) 
                     for E_surf in self.E_surfaces for osp in E_surf['onshell_propagators'] )

            if include_push_away_soft:
                objective += sum( cvxpy.norm(coord, 2) for coord in self.E_surfaces for coord in self.cvxpy_coordinates )

            p = cvxpy.Problem(
                cvxpy.Minimize( objective ), 
                [ E_surf['cxpy_expression'] <= 0 for E_surf in self.E_surfaces[:] ] )

            cvxpyres=p.solve()
            #print(cvxpyres, p.status)
            try:
                cvxpyres=p.solve()
            except cvxpy.SolverError:
                p.solve(solver=cvxpy.SCS, eps=1e-9)
            except Exception as e:
                continue

            new_seed_point = [ Vector([float(c.value[0]), float(c.value[1]), float(c.value[2])]) for c in self.cvxpy_coordinates ]

        if best_seed_so_far is None:
            if self.frozen_momenta is not None:
                new_seed_point[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']
            return new_seed_point, solved
        else:
            if self.frozen_momenta is not None:
                best_seed_so_far[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']
            return best_seed_so_far, solved

    def improve_seed_point_cardinal(self, seed_point):

        solved = False
        n_steps = 0
        new_seed_point = seed_point
        while (not solved) and n_steps <= self.maximal_cvxpy_seed_point_improvement_steps:

            for orientation in [1,-1]:
                if solved: break
                for cardinal_direction in [ (i,j) for j in range(0,3) for i in range(0, len(self.cvxpy_coordinates)) ]:

                        solution_sum = sum(abs(self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
                        consistency_test = abs(solution_sum/self.E_cm)

                        if self.debug:
                            logger.info("Step %d of improvement on the seed point\nx=%s\n that evaluates E-surface to:\n%s\nwith score %.4e by projecting onto cardinal direction %s."%(
                                n_steps+1,
                                str(new_seed_point),
                                '\n'.join('E-surface with id %d : %.16e'%(
                                    E_surf['id'], self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                                ),
                                abs(solution_sum/self.E_cm),
                                'component #%d, loop #%d'%(cardinal_direction[0],cardinal_direction[1])
                            ))

                        if consistency_test < self.consistency_check_threshold:
                            if self.debug: logger.info('Current solution satisfactory enough, returning now.')
                            solved = True
                            break
                        
                        n_steps += 1

                        new_seed_point = [ [pi+random.random()*self.E_cm for pi in p] for p in new_seed_point ]

                        for cvxpy_coord in self.cvxpy_coordinates:
                            cvxpy_coord.value = [ numpy.float64(random.random()*self.E_cm) for i in range(0,3) ]

                        cardinal_direction = [ Vector([ (0. if (i,j)!=cardinal_direction else 1.) for i in range(0,3)])*float(orientation) 
                                                for j in range(0, len(self.cvxpy_coordinates)) ]

                        new_seed_point = self.project_onto_E_surface_closest(new_seed_point, cardinal_direction, orientations=('plus',))


        all_loop_sigs = list(set([ osp['loop_sig'] for E_surf in self.E_surfaces for osp in E_surf['onshell_propagators'] ]))
        # scramble momenta that are not appearing in any of the E-surfaces selected
        for i_loop in range(0,len(self.cvxpy_coordinates)):
            if not any(sig[i_loop]!=0 for sig in all_loop_sigs):
                new_seed_point[i_loop] = [random.random()*(self.E_cm/2.) for _ in range(0,3)]

        if self.frozen_momenta is not None:
            new_seed_point[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']

        return new_seed_point, solved

    def improve_seed_point_E_surf_projection(self, seed_point):

        solved = False
        n_steps = 0
        new_seed_point = seed_point
        while (not solved) and n_steps <= self.maximal_cvxpy_seed_point_improvement_steps:

            for i_E_surf, E_surf in enumerate(self.E_surfaces):

                solution_sum = sum(abs(self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
                consistency_test = abs(solution_sum/self.E_cm)

                if self.debug:
                    logger.info("Step %d of improvement on the seed point\nx=%s\n that evaluates E-surface to:\n%s\nwith score %.4e by projecting onto E_surface #%d."%(
                        n_steps+1,
                        str(new_seed_point),
                        '\n'.join('E-surface with id %d : %.16e'%(
                            E_surf['id'], self.E_surface( new_seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                        ),
                        abs(solution_sum/self.E_cm),
                        E_surf['id']
                    ))

                if consistency_test < self.consistency_check_threshold:
                    logger.info('Current solution satisfactory enough, returning now.')
                    solved = True
                    break
                
                n_steps += 1

                E_surf_normal = [ Vector([ 
                            self.dE_surface(new_seed_point, self.E_surfaces[i_E_surf]['onshell_propagators'], self.E_surfaces[i_E_surf]['E_shift'], loop_index, component_index) 
                        for component_index in range(0,3)
                    ]) for loop_index in range(0,len(self.cvxpy_coordinates))
                ]

                new_seed_point = self.project_onto_E_surface_closest(new_seed_point, E_surf_normal, orientations=('plus','minus'))

        if self.frozen_momenta is not None:
            new_seed_point[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']

        return new_seed_point, solved

    def find_intersection(self):

        if self.debug:
            logger.info("EsurfaceIntersectionFinder will now attempt to find the intersection of the following E_surfaces:\n%s"%(
                '\n'.join(pformat(E_surf) for E_surf in self.E_surfaces)
            ))

        # First test with cvxpy if we can find a point in the interior of all specified E-surfaces.
        # This is a prerequisite for the intersection to exist (when including the boudary as part of the interior) and 
        # it is useful to server as a powerful see for the scipy optimize function
        #objective = sum(cvxpy.norm(lmom-self.cvxpy_direction, 2) for lmom in self.cvxpy_coordinates)
        #self.cvxpy_direction = cvxpy.Constant([self.E_cm*2.*i_comp for i_comp in range(0,3)])
        #objective = sum( cvxpy.norm(
        #         cvxpy.sum([coord*wgt for coord, wgt in zip(self.cvxpy_coordinates, osp['loop_sig'])])- self.cvxpy_direction, 2) 
        #         for E_surf in self.E_surfaces for osp in E_surf['onshell_propagators'] )
        objective=1
        #objective = sum( cvxpy.norm(coords-self.cvxpy_directions[i_loop],2) for i_loop, coords in enumerate(self.cvxpy_coordinates) )
        #print(objective)
        #print(self.E_surfaces[0]['cxpy_expression'])
        p = cvxpy.Problem(
            cvxpy.Minimize( objective ), 
            [ E_surf['cxpy_expression'] <= 0 for E_surf in self.E_surfaces[:] ] )

        seed_point = None
        cvxpyres=p.solve()
        #print(cvxpyres, p.status)
        try:
            cvxpyres=p.solve()
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

        if self.frozen_momenta is not None:
            seed_point[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']

        if self.maximal_cvxpy_seed_point_improvement_steps is not None and self.maximal_cvxpy_seed_point_improvement_steps > 0:
            #seed_point, solved = self.improve_seed_point_cardinal(seed_point)
            seed_point, solved = self.improve_seed_point_directions_per_signature(
                seed_point,
                abs(self.maximal_cvxpy_seed_point_improvement_steps),
                include_push_away_soft= (self.maximal_cvxpy_seed_point_improvement_steps>0)
            )

            # Only return at this step if we have set self.maximal_cvxpy_seed_point_improvement_steps negative
            if solved and (self.maximal_cvxpy_seed_point_improvement_steps<0):
                if self.debug: logger.info("The improvement of the seed point already lead to a solution for the intersection:\nx=%s\nfor which the E-surfaces evaluate as follows:\n%s"%(
                    str(seed_point),
                    '\n'.join('E-surface with id %d : %.16e'%(
                        E_surf['id'], self.E_surface( seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                    )
                ))

                return seed_point

        if self.debug: 
            logger.info("The E-surface equations evaluate as follows for the seed point found x=%s\n%s"%(
                seed_point,
                '\n'.join('E-surface with id %d : %.16e'%(
                    E_surf['id'], self.E_surface( seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                )
            ))

        for i_shift, seed_point_shift in enumerate(self.seed_point_shifts):
            
            #seed_point_shift = [[0.,]*3]*4
            shifted_seed = [ [ seed_comp+seed_point_shift[i_vec][i_comp] for i_comp, seed_comp in enumerate(seed_vec) ] for i_vec, seed_vec in enumerate(seed_point) ]

            if self.frozen_momenta is not None:
                shifted_seed[-len(self.frozen_momenta['out']):] = self.frozen_momenta['out']

            if self.debug: 
                logger.info("An intersection is possible, scipy.optimize.root will now attempt to find an intersection from the scipy seed point with shift #%d: %s."%(i_shift,str(shifted_seed)))
                logger.info("The E-surface equations evaluate as follows for this scipy seed point:\n%s"%(
                    '\n'.join('E-surface with id %d : %.16e'%(
                        E_surf['id'], self.E_surface( shifted_seed, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                    )
                ))

            if self.frozen_momenta is not None:
                shifted_seed = shifted_seed[:-len(self.frozen_momenta['out'])]

            scipy_seed = sum(shifted_seed,[])

            if self.use_scipy_root:
                #pprint(self.intersection_function(scipy_seed))
                #pprint(self.intersection_function_jac(scipy_seed))
                scipy_res = optimize.root(
                    self.intersection_function,
                    scipy_seed,
                    jac=self.intersection_function_jac, 
                    method=self.scipy_method,
                    options={
                        'xtol'   : self.scipy_tolerance,
                        'maxfev' : 1000*(len(scipy_seed)+1)
                    }
                )

                if self.debug:
                    logger.info("The attempt at finding an intersection with scipy.optimize.root for seed shift #%d yielded:\n%s"%(i_shift,str(scipy_res)))

                #logger.debug(str(scipy_res))
                if not scipy_res.success:
                    if self.debug: logger.info("Could not find an intersection with scipy.optimize.root for seed shift #%d"%i_shift)
                else:
                    #pprint(self.intersection_function(scipy_res.x))
                    #pprint(self.intersection_function_jac(scipy_res.x))
                    intersection_point = [ scipy_res.x[i:i+3] for i in range(0,len(scipy_res.x),3) ]
                    if self.frozen_momenta is not None:
                        intersection_point += self.frozen_momenta['out']
                    # Double-check the solution
                    solution_sum = sum(abs(self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
                    consistency_test = abs(solution_sum/self.E_cm)
                    if consistency_test > self.consistency_check_threshold:
                        logger.warning("Scipy root found an intersection, however the testing of that solution was not satisfactory. %.3e > %.3e"%(
                            consistency_test, self.consistency_check_threshold))
                    else:
                        if self.debug: logger.info("The E-surface equations evaluate as follows for the intersection point found using shift #%d x=%s\n%s"%(
                            i_shift,
                            intersection_point,
                            '\n'.join('E-surface with id %d : %.16e'%(
                                E_surf['id'], self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] ) ) for E_surf in self.E_surfaces
                            )
                        ))
                        return intersection_point

            if self.use_scipy_fsolve:

                if self.debug: logger.info("scipy.optimize.fsolve will now attempt to find an intersection from the seed point with shift #%d: %s."%(i_shift,str(scipy_seed)))

                #pprint(self.intersection_scalar_function(scipy_seed))
                #pprint(self.intersection_scalar_function_jac(scipy_seed))
                scipy_res = None
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        scipy_res = optimize.fsolve(
                            self.intersection_function,
                            scipy_seed,
                            fprime=self.intersection_function_jac, 
                            **{
                                'xtol'   : self.scipy_tolerance,
                                'maxfev' : 1000*(len(scipy_seed)+1)
                            }
                        )
                    fsolve_success = True
                except Warning as e:
                    fsolve_success = False

                if self.debug:
                    logger.info("The attempt at finding an intersection with scipy.optimize.fsolve for seed shift #%d yielded:\n%s"%(i_shift,str(scipy_res)))

                #logger.debug(str(scipy_res))
                if not fsolve_success:
                    if self.debug: logger.info("Could not find an intersection with scipy.optimize.fsolve for seed shift #%d"%i_shift)
                else:
                    #pprint(self.intersection_function(scipy_res.x))
                    #pprint(self.intersection_function_jac(scipy_res.x))
                    intersection_point = [ scipy_res[i:i+3] for i in range(0,len(scipy_res),3) ]
                    if self.frozen_momenta is not None:
                        intersection_point += self.frozen_momenta['out']
                    # Double-check the solution
                    solution_sum = sum(abs(self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
                    consistency_test = abs(solution_sum/self.E_cm)
                    if consistency_test > self.consistency_check_threshold:
                        logger.warning("Scipy fsolve found an intersection, however the testing of that solution was not satisfactory. %.3e > %.3e"%(
                            consistency_test, self.consistency_check_threshold))
                    else:
                        if self.debug: logger.info("The E-surface equations evaluate as follows for the intersection point found using shift #%d x=%s\n%s"%(
                            i_shift,
                            intersection_point,
                            '\n'.join('E-surface with id %d : %.16e'%(
                                E_surf['id'], self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] ) ) for E_surf in self.E_surfaces
                            )
                        ))
                        return intersection_point

            if self.use_scipy_minimize:

                if self.debug: logger.info("scipy.optimize.minimize will now attempt to find an intersection from the seed point with shift #%d: %s."%(i_shift,str(scipy_seed)))

                #pprint(self.intersection_scalar_function(scipy_seed))
                #pprint(self.intersection_scalar_function_jac(scipy_seed))
                scipy_res = None
                try:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error')
                        scipy_res = optimize.minimize(
                            self.intersection_scalar_function,
                            scipy_seed,
                            tol=self.scipy_tolerance,
                            jac=self.intersection_scalar_function_jac,
                            hess=self.intersection_scalar_function_hessian,
                            method='trust-ncg'#'trust-krylov'#'dogleg'#'trust-ncg'#'Newton-CG'#'BFGS'
                        )
                    fsolve_success = True
                except Warning as e:
                    fsolve_success = False

                if self.debug:
                    logger.info("The attempt at finding an intersection with scipy.optimize.minimize for seed shift #%d yielded:\n%s"%(i_shift,str(scipy_res)))

                #logger.debug(str(scipy_res))
                if not fsolve_success:
                    if self.debug: logger.info("Could not find an intersection with scipy.optimize.minimize for seed shift #%d"%i_shift)
                else:
                    #pprint(self.intersection_function(scipy_res.x))
                    #pprint(self.intersection_function_jac(scipy_res.x))
                    intersection_point = [ scipy_res.x[i:i+3] for i in range(0,len(scipy_res.x),3) ]
                    if self.frozen_momenta is not None:
                        intersection_point += self.frozen_momenta['out']
                    # Double-check the solution
                    solution_sum = sum(abs(self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces)
                    consistency_test = abs(solution_sum/self.E_cm)
                    if consistency_test > self.consistency_check_threshold:
                        #logger.warning("Scipy minimize found an intersection, however the testing of that solution was not satisfactory. %.3e > %.3e"%(
                        #    consistency_test, self.consistency_check_threshold))
                        pass
                    else:
                        if self.debug: logger.info("The E-surface equations evaluate as follows for the intersection point found using shift #%d x=%s\n%s"%(
                            i_shift,
                            intersection_point,
                            '\n'.join('E-surface with id %d : %.16e'%(
                                E_surf['id'], self.E_surface( intersection_point, E_surf['onshell_propagators'], E_surf['E_shift'] ) ) for E_surf in self.E_surfaces
                            )
                        ))
                        return intersection_point


        if self.maximal_cvxpy_find_intersection_steps is not None and self.maximal_cvxpy_find_intersection_steps > 0:
            #seed_point, solved = self.improve_seed_point_cardinal(seed_point)
            seed_point, solved = self.improve_seed_point_directions_per_signature(
                seed_point,
                self.maximal_cvxpy_find_intersection_steps,
                include_push_away_soft=False
            )

            if solved:
                if self.debug: logger.info("The second pass of cvxpy wtihout push away soft found a solution for the intersection:\nx=%s\nfor which the E-surfaces evaluate as follows:\n%s"%(
                    str(seed_point),
                    '\n'.join('E-surface with id %d : %.16e'%(
                        E_surf['id'], self.E_surface( seed_point, E_surf['onshell_propagators'], E_surf['E_shift'] )) for E_surf in self.E_surfaces
                    )
                ))
                return seed_point
    
        if self.debug: logger.info("Could not find an intersection with scipy.optimize with %d attempts of shifting the seed point."%len(self.seed_point_shifts))

        return None

    @classmethod
    def delta(cls, loop_momenta, loop_sig, shift, m_squared):
        
        k = [ sum( l[i]*factor for l, factor in zip(loop_momenta,loop_sig)) for i in range(0,3) ]

        return math.sqrt(
             (k[0]**2+k[1]**2+k[2]**2)+
            +2*(k[0]*shift[0]+k[1]*shift[1]+k[2]*shift[2])
            +shift[0]**2+shift[1]**2+shift[2]**2
            +m_squared
        )
    
    @classmethod
    def ddelta(cls, loop_momenta, loop_sig, shift, m_squared, loop_index, component_index):

        k_derived = loop_momenta[loop_index][component_index]*loop_sig[loop_index]

        if k_derived == 0.:
            return 0.

        k_not_derived = sum( l[component_index]*factor for i_loop_mom, (l, factor) in enumerate(zip(loop_momenta,loop_sig)) if i_loop_mom!=loop_index )

        return loop_sig[loop_index]*(k_derived + k_not_derived + shift[component_index])/cls.delta(loop_momenta, loop_sig, shift, m_squared)

    @classmethod
    def dddelta(cls, loop_momenta, loop_sig, shift, m_squared, loop_indexA, component_indexA,loop_indexB, component_indexB):

        if (loop_indexA, component_indexA)==(loop_indexB, component_indexB):
            return (1.-(cls.ddelta(loop_momenta, loop_sig, shift, m_squared,loop_indexA, component_indexA)**2))/cls.delta(loop_momenta, loop_sig, shift, m_squared)
        else:
            return ((cls.ddelta(loop_momenta, loop_sig, shift, m_squared,loop_indexA, component_indexA)*cls.ddelta(loop_momenta, loop_sig, shift, m_squared,loop_indexB, component_indexB))/
                cls.delta(loop_momenta, loop_sig, shift, m_squared) )

    @classmethod
    def E_surface(cls, loop_momenta, onshell_propagators, E_surface_shift):
        
        # Not necessary any longer as I append frozen momenta upstream now.
        #if len(loop_momenta) < len(self.cvxpy_coordinates):
        #    loop_momenta += self.frozen_momenta['out']

        return sum( cls.delta(loop_momenta,osp['loop_sig'],osp['v_shift'],osp['m_squared']) for osp in onshell_propagators) + E_surface_shift

    @classmethod
    def dE_surface(cls, loop_momenta, onshell_propagators, E_surface_shift, loop_index, component_index):

        return sum( cls.ddelta(loop_momenta,osp['loop_sig'],osp['v_shift'],osp['m_squared'],loop_index, component_index) for osp in onshell_propagators)

    @classmethod
    def ddE_surface(cls, loop_momenta, onshell_propagators, E_surface_shift, loop_indexA, component_indexA,loop_indexB, component_indexB):

        return sum( cls.dddelta(loop_momenta,osp['loop_sig'],osp['v_shift'],osp['m_squared'],loop_indexA, component_indexA,loop_indexB, component_indexB) for osp in onshell_propagators)

    def intersection_function(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]
        if self.frozen_momenta is not None:
            loop_momenta += self.frozen_momenta['out']

        res = []
        for E_surf in self.E_surfaces:
            res.append( self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift']) )

        # Pad remaining entries with no constraints if there are some left.
        res += [0.]*(len(xs)-len(res))

        return res

    def intersection_function_jac(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]
        if self.frozen_momenta is not None:
            loop_momenta += self.frozen_momenta['out']

        res = []
        for E_surf in self.E_surfaces:
            res.append([])
            for loop_index in range(0,len(loop_momenta)-self.n_frozen_momenta):
                for component_index in range(0,3):
                    res[-1].append( self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index, component_index) )

        # Pad remaining entries with no constraints if there are some left.
        res += [ [0.]*len(xs) ]*(len(xs)-len(res))

        return res

    def intersection_scalar_function(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]
        if self.frozen_momenta is not None:
            loop_momenta += self.frozen_momenta['out']

        res = 0.
        for E_surf in self.E_surfaces:
            res += self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'])**2

        return res

    def intersection_scalar_function_jac(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]
        if self.frozen_momenta is not None:
            loop_momenta += self.frozen_momenta['out']

        res = []
        for loop_index in range(0,len(loop_momenta)-self.n_frozen_momenta):
            for component_index in range(0,3):
                res.append(0.)
                for E_surf in self.E_surfaces:
                    res[-1] += ( 2.0*self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'])*
                                   self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index, component_index) )

        return res

    def intersection_scalar_function_hessian(self, xs):

        loop_momenta = [ xs[i:i+3] for i in range(0,len(xs),3) ]
        if self.frozen_momenta is not None:
            loop_momenta += self.frozen_momenta['out']

        res = []
        for loop_index_A in range(0,len(loop_momenta)-self.n_frozen_momenta):
            for component_index_A in range(0,3):
                res.append([])
                for loop_index_B in range(0,len(loop_momenta)-self.n_frozen_momenta):
                    for component_index_B in range(0,3):
                        res[-1].append(0.)
                        for E_surf in self.E_surfaces:
                            if (loop_index_A,component_index_A)==(loop_index_B,component_index_B):
                                res[-1][-1] += (2.*(self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index_A, component_index_A)**2)+
                                                2.*self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'])*
                                                self.ddE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index_A, component_index_A,loop_index_A, component_index_A)
                                            )
                            else:
                                res[-1][-1] += (
                                    2.*self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index_A, component_index_A)*self.dE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index_B, component_index_B)
                                   +2.*self.E_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'])*self.ddE_surface(loop_momenta, E_surf['onshell_propagators'], E_surf['E_shift'], loop_index_A, component_index_A,loop_index_B, component_index_B)
                                )

        return res        
