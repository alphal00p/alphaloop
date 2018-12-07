#!/usr/bin/env python

import sys

import os
pjoin = os.path.join

root_path = os.path.dirname(os.path.realpath( __file__ ))
sys.path.insert(0, pjoin(root_path,os.pardir,os.pardir))
sys.path.insert(0, pjoin(root_path,'rust_backend'))

import math
import time
import random
import yaml
from yaml import Loader, Dumper
from pprint import pprint

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

import madgraph.integrator.vectors as vectors

try:
    import deformation
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'deformation' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory." )

class Propagator(object):

    def __init__(self, id=None, loop_momenta_signs = (1,1), q=vectors.LorentzVector(), m_sq=0.):
        self.id                 = id
        # Propagator is ((\sum_i sgn_loop_momentum_i * l_i) + q)^2 -m^2
        # Notice that the sequence of loop momenta signs can only contain -1, 0 or 1
        self.loop_momenta_signs = loop_momenta_signs
        self.q                  = q
        self.m_sq               = m_sq

    def get_n_contributing_loop_momenta(self):
        return sum(ls != 0 for ls in self.loop_momenta_signs)

    def get_first_contributing_loop_momentum_index(self):
        for i_loop_mom, ls in enumerate(self.loop_momenta_signs):
            if ls != 0:
                return i_loop_mom
        # No loop momentum dependence
        return None

    def put_onshell(self, direction_vectors, offset_vectors, lambdas, solution_sign=1):
        """ Given direction (ld, kd) and offset vectors (L,K), specifying:
                (l,k) (\lambda1, \lambda2) = (ld*\lambda1 + L, kd*\lambda2 + K)
         returns a lambda function for the solution function \lambda2*=F(\lambda1) such that the propgator D satisfies:
                D( l(\lambda1, F(\lambda1)) , k(\lambda1, F(\lambda1)) ) == 0
        This function will substitute the first loop momentum that this propagator depends on with the value that
        puts this propagator on-shell, and return its index in the list as an integer.
        """

        assert(solution_sign in [-1,1])

        # First identify the loop momentum index that will need to be fixed
        # flmi == first_loop_momentum_index
        flmi = self.get_first_contributing_loop_momentum_index()

        # Then cast this polynomial in the form, choosing l to be the first contributing loop momentum index
        # (\lambda * l^\mu + R^\mu)^2 - m_sq == \lambda^2 l\dotl + \lambda * 2 * l\dotR + R\dotR - m_sq
        # Start with the offset from the flmi vector as well as the constant + q
        Rvec = offset_vectors[flmi]*self.loop_momenta_signs[flmi] + self.q
        for i_lambda, lambda_value in enumerate(lambdas):
            if i_lambda==flmi or self.loop_momenta_signs[i_lambda]==0:
                continue
            Rvec += (direction_vectors[i_lambda]*lambda_value + offset_vectors[i_lambda])*self.loop_momenta_signs[i_lambda]

        # Compute the coefficient of the quadratic equation ' A \lambda^2 + B \lambda + C == 0 '
        A = direction_vectors[flmi].dot(direction_vectors[flmi])
        B = 2.0*self.loop_momenta_signs[flmi]*direction_vectors[flmi].dot(Rvec)
        C = Rvec.dot(Rvec) - self.m_sq

        # Now modify the appropriate lambda value with the solution
        discriminant = B**2 - 4.*A*C
        if discriminant < 0.:
            lambdas[flmi] = None
            # No real solution, return None
            return None
        else:
            lambda_onshell = (-B + solution_sign*math.sqrt(discriminant))/(2.*A)
            lambdas[flmi] = lambda_onshell
            return flmi

    def evaluate(self, loop_momenta_real, loop_momenta_complex=None):

        real_part = self.q.get_copy()
        for i_momentum, l_real in enumerate(loop_momenta_real):
            if self.loop_momenta_signs[i_momentum] != 0:
                real_part += l_real*self.loop_momenta_signs[i_momentum]

        imag_part = None
        if loop_momenta_complex is not None:
            imag_part = vectors.LorentzVector()
            for i_momentum, l_imag in enumerate(loop_momenta_complex):
                imag_part += l_imag * self.loop_momenta_signs[i_momentum]

        result_real = real_part.dot(real_part) - self.m_sq
        if imag_part is None:
            return result_real

        result_real -= imag_part.dot(imag_part)
        result_imag = 2.0*real_part.dot(imag_part)

        return complex(result_real, result_imag)

class Topology(object):

    # Add below the new topology IDs you are interested in
    _RUST_TOPOLOGY_MAP = {
        'double-box'                : 1,
        'triangle-box'              : 2,
        'triangle-bo-alternative'   : 3,
        'diagonal-box'              : 4,
        'double-box-SB'             : 5,
    }

    def __init__(self, topology_name='double-box-SB', PS_point=None):
        self.name           = topology_name
        self.PS_point       = PS_point

        if topology_name!='double-box-SB':
            raise BaseException('Only topology double-box-SB supported for now.')
        self.propagators    = self.get_double_box_SB_propagators()

    def get_RUST_topology_ID(self):
        try:
            return self._RUST_TOPOLOGY_MAP[self.name]
        except KeyError:
            raise BaseException("Topology with name '%s' is unknown to the rust backend. "%self.name+
                                "Consider adding it manually in the _RUST_TOPOLOGY_MAP dictionary.")

    def get_masses_sq(self):
        return [p.m_sq for p in self.propagators]

    def get_qs(self):
        return []

    def get_e_cm_sq(self):

        initial_momenta_summed = vectors.LorentzVector()
        for v in self.PS_point:
            if v[0]>0.:
                initial_momenta_summed += v
        res = initial_momenta_summed.dot(initial_momenta_summed)
        assert(res > 0.)
        return res

    def get_double_box_SB_propagators(self):
        return [
            Propagator(
                id                  = (23,2),
                loop_momenta_signs  = (1,0),
                q                   = self.PS_point[2],
                m_sq                = 0.
            ),
            Propagator(
                id                  = (2,1),
                loop_momenta_signs  = (1,0),
                q                   = self.PS_point[2]+self.PS_point[1],
                m_sq                = 0.
            ),
            Propagator(
                id                  = (1,14),
                loop_momenta_signs  = (1,0),
                q                   = self.PS_point[3]*-1.,
                m_sq                = 0.
            ),
            Propagator(
                id                  = (14,4),
                loop_momenta_signs  = (0,1),
                q                   = self.PS_point[3]*-1.,
                m_sq                = 0.
            ),
            Propagator(
                id                  = (4,3),
                loop_momenta_signs  = (0,1),
                q                   = vectors.LorentzVector(),
                m_sq                = 0.
            ),
            Propagator(
                id                  = (3,23),
                loop_momenta_signs  = (0,1),
                q                   = self.PS_point[2],
                m_sq                = 0.
            ),
            Propagator(
                id                  = (23,14),
                loop_momenta_signs  = (-1,1),
                q                   = vectors.LorentzVector(),
                m_sq                = 0.
            )
        ]

class PoleScanner(object):

    @staticmethod
    def map_to_infinity(x):
        return ((1./(1.-x)) - 1./x)

    @staticmethod
    def map_from_infinity(x):
        return -2./(-2.+x-math.sqrt(4.+x**2))

    def __init__(self,
                 log_stream         =   None,
                 topology           =   'double-box-SB',
                 PS_point           =   None,
                 direction_vectors  =   None,
                 offset_vectors     =   None,
                 t_values           =   None,
                 fixed_t_values     =   None,
                 deformation_configuration = None,
                ):
        self.PS_point           = PS_point
        self.topology           = Topology(topology_name=topology, PS_point=self.PS_point)
        self.direction_vectors  = direction_vectors
        self.offset_vectors     = offset_vectors
        # Promote the t-values (in [0,1]) to tuples (t, \lambda(t)) where \lambda(t) is their equivalent
        # on the infinite segment [-inf, inf]
        self.t_values           = [(t, PoleScanner.map_to_infinity(t)) for t in t_values]

        self.fixed_t_values     = [(t, PoleScanner.map_to_infinity(t)) for t in fixed_t_values]

        self.deformation_configuration = deformation_configuration
        # Now save the default yaml configuration to the log stream if specified.
        self.configuration = {
                    'topology'          :   self.topology.name,
                    'PS_point'          :   [[float(v) for v in vec] for vec in self.PS_point],
                    'direction_vectors' :   [[float(v) for v in vec] for vec in self.direction_vectors],
                    'offset_vectors'    :   [[float(v) for v in vec] for vec in self.offset_vectors],
                    'deformation_configuration' : deformation_configuration
        }
        if log_stream is not None:
            log_stream.write(yaml.dump([('configuration', self.configuration ),], Dumper=Dumper))

        # Now setup the deformation
        self.deformation = deformation.Deformation(
            e_cm_sq     = float(self.topology.get_e_cm_sq()),
            mu_sq       = float(self.deformation_configuration['mu_sq']),
            m1_fac      = float(self.deformation_configuration['m1_fac']),
            m2_fac      = float(self.deformation_configuration['m2_fac']),
            m3_fac      = float(self.deformation_configuration['m3_fac']),
            m4_fac      = float(self.deformation_configuration['m4_fac']),
            gamma1      = float(self.deformation_configuration['gamma1']),
            gamma2      = float(self.deformation_configuration['gamma2']),
            soft_fac    = float(self.deformation_configuration['soft_fac']),
            region      = 0,
            # For double-box-SB the q's are computed internally, we only want to set the externals here
            qs_py       = [[0.]*4]*7,
            masses      = [float(math.sqrt(m_sq)) for m_sq in self.topology.get_masses_sq()]
        )
        self.deformation.set_external_momenta([[v for v in vec] for vec in self.PS_point])

    def scan(self, propagators = None):
        all_results = []

        for prop in self.topology.propagators:
            # Loop over both onshell solutions that each lambda value choice can offer.
            for solution_sign in [-1,1]:
                if propagators is not None and (prop.id,solution_sign) not in propagators:
                    # This propagator should be ignored
                    continue

                one_result = {
                    'propagator_id'         :   (prop.id, solution_sign),
                    'pole_imaginary_parts'  :   []
                }

                # Having potentially more than one loop-momentum dependence, this onshellness condition will yield a function
                # of a "varying" lambda being always the *first* loop momentum with sign 0 in the propagator or the *second*
                # with abs(sign) 1, whichever comes first (we could generalise this when doing more than one-loop).
                varying_lambda_index = None
                n_abs_one_count = 0
                for lambda_index, l_sign in enumerate(prop.loop_momenta_signs):
                    if l_sign == 0:
                        varying_lambda_index = lambda_index
                        break
                    else:
                        n_abs_one_count += 1
                        if n_abs_one_count>=2:
                            varying_lambda_index = lambda_index
                            break
                if varying_lambda_index is None:
                    raise BaseException('ERROR: Could not find a suitable lambda to scan over onshell propagator with id: %s'%str(prop.id))

                # Now scan various values
                for t, lambda_value in self.t_values:
                    lambda_onshell_values = [l for (_,l) in self.fixed_t_values]
                    # Adjust the lambda to vary to current value
                    lambda_onshell_values[varying_lambda_index] = lambda_value

                    modified_lambda_index = prop.put_onshell(self.direction_vectors, self.offset_vectors,
                                                                        lambda_onshell_values, solution_sign=solution_sign)
                    if modified_lambda_index is None:
                        # Solution is complex, no data point here:
                        one_result['pole_imaginary_parts'].append( (
                            ( t, tuple( (float(lov) if lov is not None else None) for lov in lambda_onshell_values) ) ,
                                                                                                                 None ))
                        continue

                    # Now compute the real loop momenta prior deformation
                    loop_momenta = [dir_vec*lam + off_vec for dir_vec,lam,off_vec in
                                        zip(self.direction_vectors, lambda_onshell_values, self.offset_vectors)]

                    # Sanity check that the propagator is indeed on-shell. Normalise to some quantity with the same dimension
                    # and always positive
                    normalisation_factor = max(
                        (  sum(sum(abs(le) for le in lm) for lm in loop_momenta) /
                           sum(sum(1. for le in lm) for lm in loop_momenta)   ) ,
                        1.0e-99)
                    threshold = 1.0e-10
                    prop_evaluated = prop.evaluate(loop_momenta)
                    if abs(prop_evaluated/(normalisation_factor**2)) > threshold:
                        print(
                            "WARNING:: Propagator with id '%s' was incorrectly put onshell: abs(onshell) / norm = %.3e / %.3e = %.3e > %.3e"%(
                            prop.id, abs(prop_evaluated), (normalisation_factor**2),
                            abs(prop_evaluated)/(normalisation_factor**2), threshold)
                        )

                    # Now compute the deformed momenta
                    deformed_momenta = [complex(*x) for x in self.deformation.deform_two_loops(
                        self.topology.get_RUST_topology_ID(),
                        [[v for v in loop_momenta[0]], [v for v in loop_momenta[1]]])[0]]
                    deformed_momenta_real = [ vectors.LorentzVector([x.real for x in deformed_momenta[:4]]),
                                              vectors.LorentzVector([x.real for x in deformed_momenta[4:]]) ]
                    deformed_momenta_imag = [ vectors.LorentzVector([x.imag for x in deformed_momenta[:4]]),
                                              vectors.LorentzVector([x.imag for x in deformed_momenta[4:]]) ]

                    # On the plot produced we will register on the x-axis the t-value corresponding to the lambda scanned.
                    # For completeness, we also provide as a second entry the list all lambda values for that scanned point
                    prop_deformed = prop.evaluate(deformed_momenta_real, deformed_momenta_imag)
                    one_result['pole_imaginary_parts'].append( (  ( t, tuple(float(lov) for lov in lambda_onshell_values) ),
                                                                                            float(prop_deformed.imag) ) )
                    #if prop_deformed.imag < 0.:
                    #    print('WARNING: An onshell propagator evaluates to a negative imaginary part!:\n%s'%str(one_result))

                # Save the result into the logstream if specified
                if log_stream is not None:
                    log_stream.write(yaml.dump([( 'propagator_result', one_result ), ], Dumper=Dumper))
                all_results.append(one_result)

        return {
            'configuration' : self.configuration,
            'scan_results'  : all_results
        }

class ResultsAnalyser(object):
    
    def __init__(self, results):
        # Each entry of the scan_results has the following format:
        ###    {
        ###     'propagator_id': (prop.id, solution_sign),
        ###     'pole_imaginary_parts': [
        ###       ( (t_value, all_lambda_values), prop_imag_part )
        ###     ]
        ###    }
        self.scan_results   = results['scan_results']
        self.configuration  = results['configuration']
        self.direction_vectors = [ vectors.LorentzVector(lv) for lv in self.configuration['direction_vectors'] ]
        self.offset_vectors    = [ vectors.LorentzVector(ov) for ov in self.configuration['offset_vectors'] ]
        self.run_time       = results['run_time']

    def test_imaginary_parts_in_scan(self):
        found_problem = 0
        for entry in self.scan_results:
            for (t_value, lambda_values), prop_imag_part in entry['pole_imaginary_parts']:
                if prop_imag_part is None:
                    continue
                if prop_imag_part <= 0.:
                    found_problem += 1
                    loop_momenta = [dir_vec*lam + off_vec for dir_vec,lam,off_vec in
                                        zip(self.direction_vectors, lambda_values, self.offset_vectors)]
                    print('WARNING: The onshell propagator %s evaluates to a negative imaginary part (%s) for the following loop momenta:'%(
                        entry['propagator_id'],prop_imag_part
                    ))
                    for i_l, lm in enumerate(loop_momenta):
                        print('    Loop momentum #%d : %s'%(i_l, lm))
        return found_problem

    def generate_plots(self):

        plt.title('Onshell propagators test')
        plt.ylabel('Imaginary part')
        plt.xlabel('t')
        plt.yscale('log')

        # Generate the line data for each propagator
        propagator_lines = [ ('Prop. %s@%d'%entry['propagator_id'], zip(
              *[(t_value, prop_imag_part) for (t_value, lambda_values), prop_imag_part in entry['pole_imaginary_parts']]
            )) for entry in self.scan_results ]

        for line_name, (x_data, y_data) in propagator_lines:
            plt.plot(x_data, y_data, label=line_name)
        plt.legend()
        plt.show()

    def analyse(self, n_optimal_to_show=5):

        print('='*80)
        print("run time: %.3e [h]"%(self.run_time/3600.0))
        print('-'*80)
        found_problem = self.test_imaginary_parts_in_scan()
        if found_problem == 0:
            print('OK:: All imaginary parts of onshell propagators tested are positive.')
        else:
            print('ERROR:: Some imaginary part of onshell propagators tests were found to be negative. The deformation is likely wrong.')
        print('='*80)
        self.generate_plots()

def load_results_from_yaml(log_file_path):
    """Load a full-fledged scan from a yaml dump"""

    raw_data = yaml.load(open(log_file_path,'r'), Loader=Loader)
    scan_results = []
    processed_data = {}
    for entry_name, entry_value in raw_data:
        if entry_name == 'propagator_result':
            scan_results.append(entry_value)
            continue
        processed_data[entry_name] = entry_value
    processed_data['scan_results'] = scan_results
    return processed_data

if __name__ == '__main__':

    prog_options = ' '.join(sys.argv[1:]).split('--')

    topology            = 'double-box-SB'
    # Default values assume a two-loop topology
    # None means chosen at random
    direction_vectors   = [None,None]
    offset_vectors      = [None,None]
    t_n_points          = 100
    t_range             = [0.,1.0]
    # If Non, it will be assigned to a default later on
    t_values            = None
    fixed_t_values      = [0.5, 0.5]
    # Specifying None for propagators indicated that they should all be considered
    propagators         = None
    # Offers the possibility of looping 'do_loop' times (or at infinity when set to -1) the test with different random
    # direction and offset vectors each time.
    do_loop             = None
    # For now the selected PS Point will have to be specified below
    PS_point            = vectors.LorentzVectorList([
        vectors.LorentzVector([19.6586, -7.15252, -0.206016, 8.96383]),
        vectors.LorentzVector([26.874, 7.04203, -0.0501295, -12.9055]),
        vectors.LorentzVector([43.4674, 0.11049, 0.2561455, 3.94167]),
        vectors.LorentzVector([-90.0, 0.0, 0.0, 0.0]),
    ])
    # If specified to None, then a random PS Point is generated.
    #PS_point            = None #TOIMPLEMENT
    load_results_from   = None
    save_results_to     = pjoin(os.getcwd(),'poles_imaginary_part_scan.yaml')
    config_file_path    = pjoin(os.getcwd(),'rust_backend','config.yaml')

    for prog_option in prog_options:
        if prog_option.strip()=='':
            continue
        try:
            option_name, option_value = prog_option.strip().split('=')
        except ValueError:
            option_name, option_value = prog_option, None

        if option_name in ['topology','t']:
            topology = option_value
        elif option_name in ['config_file_path', 'c']:
            config_file_path = option_value
        elif option_name in ['seed','s']:
            random.seed(eval(option_value))
        elif option_name in ['direction_vectors','d']:
            # A vector set to None means it will be randomly chosen
            direction_vectors = [ (vectors.LorentzVector(vec) if vec is not None else vec) for vec in eval(option_value) ]
        elif option_name in ['offset_vectors','o']:
            # A vector set to None means it will be randomly chosen
            offset_vectors = [ (vectors.LorentzVector(vec) if vec is not None else vec) for vec in eval(option_value) ]
        elif option_name in ['t_values','tv']:
            # The ranges must be iterable of floats or None to select the default iterable.
            t_values = eval(option_value)
        elif option_name in ['t_range','tr']:
            t_range = eval(option_value)
        elif option_name in ['t_n_points','tn']:
            t_n_points = eval(option_value)
        elif option_name in ['fixed_t_values','ft']:
            # The ranges must be iterable of floats or None to select the default iterable.
            fixed_t_values = eval(option_value)
            if fixed_t_values is None:
                if len(direction_vectors)>0:
                    fixed_t_values = [0.5,]*len(direction_vectors)
                else:
                    # Assume two-loop momenta, as in the default fixed_t_values value set above
                    fixed_t_values = [0.5,0.5]
        elif option_name in ['propagators','p']:
            # The propagators option must be an iterable yielding integers or None to select them all
            propagators = eval(option_value)
        elif option_name in ['load_results_from','load']:
            load_results_from = option_value
        elif option_name in ['save_results_to','save']:
            save_results_to = option_value
        elif option_name in ['loop','l']:
            do_loop = eval(option_value)
        else:
            raise BaseException("Parameter '%s' not reckognised."%option_name)

    try:
        deformation_configuration = yaml.load(open(config_file_path, 'r'), Loader=Loader)
    except Exception as e:
        print("ERROR: Could not parse yaml configuration file at '%s'." % config_file_path)
        raise e

    # Generate the list of t_values if None
    if t_values is None:
        t_values = [ ( t_range[0]+ ( ( (t_range[1]-t_range[0]) / float(t_n_points+1)) * i ) ) for i in range(1, t_n_points+1)]

    n_loops_performed = 0
    while True:
        try:
            # post-process the vector definitions by replacing the ones that are None and should be taken random.
            direction_vectors = [
                (vec if vec is not None else vectors.LorentzVector([(random.random()-0.5)*2 for _ in range(4)]))
                for vec in direction_vectors]
            offset_vectors = [
                (vec if vec is not None else vectors.LorentzVector([(random.random() - 0.5) * 2 for _ in range(4)]))
                for vec in offset_vectors]
            # Normalise the direction vectors
            direction_vectors = [vec*(1.0/max(abs(v) for v in vec)) for vec in direction_vectors]

            if load_results_from is None:
                log_stream = open(save_results_to,'w')
                scanner = PoleScanner(
                    log_stream              =   log_stream,
                    topology                =   topology,
                    PS_point                =   PS_point,
                    direction_vectors       =   direction_vectors,
                    offset_vectors          =   offset_vectors,
                    t_values                =   t_values,
                    fixed_t_values          =   fixed_t_values,
                    deformation_configuration = deformation_configuration
                )
                start_time = time.time()
                test_results = scanner.scan(propagators=propagators)
                test_results['run_time'] = time.time()-start_time
                log_stream.write(yaml.dump([('run_time', test_results['run_time']), ], Dumper=Dumper))
                log_stream.close()
            else:
                try:
                    test_results = load_results_from_yaml(load_results_from)
                except Exception as e:
                    print("ERROR: Could not hyperparameters scan results from file '%s'. Error: %s"%(
                        load_results_from, str(e) ))
                    raise e

            #pprint(test_results)
            analyser = ResultsAnalyser(test_results)

            n_loops_performed += 1
            if do_loop is not None:
                analyser.test_imaginary_parts_in_scan()
            else:
                analyser.analyse()

            if (do_loop is None) or (do_loop is not None and do_loop > 0 and n_loops_performed>=do_loop):
                break

            if n_loops_performed%10 == 0:
                print('Number of test scans performed: %d'%n_loops_performed)

        except KeyboardInterrupt:
            print('Run aborted by user. Number of test scans performed so far: %d' % n_loops_performed)
            break