#!/usr/bin/env python3

import os
import argparse
import datetime as dt
from pySecDec.loop_integral import *
from pySecDec.integral_interface import IntegralLibrary
pjoin = os.path.join
import logging
import subprocess
import multiprocessing
import time
import sys
import shutil
import sympy as sp
import math
from pprint import pprint, pformat

root_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

class MyFormatter(logging.Formatter):
    converter=dt.datetime.fromtimestamp
    def formatTime(self, record, datefmt=None):
        ct = self.converter(record.created)
        if datefmt:
            s = '%s%s%s'%('\033[94m',ct.strftime(datefmt),'\033[0m')
        else:
            t = ct.strftime("%Y-%m-%d %H:%M:%S")
            s = "%s%s,%03d%s" % ('\033[94m',t, record.msecs,'\033[0m')
        return s

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

console = logging.StreamHandler()
logger.addHandler(console)

formatter = MyFormatter(fmt='%(asctime)s %(message)s',datefmt='%Y-%m-%d,%H:%M:%S.%f')
console.setFormatter(formatter)

class Silence:
    
    def __init__(self,active=True):
        self.active = active

    def __enter__(self):
        if not self.active: return
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.active: return
        sys.stdout.close()
        sys.stdout = self._original_stdout
        sys.stderr.close()
        sys.stderr = self._original_stderr

def generate(graph_name, quiet=False):
    
    start_time = time.time()
    old_level = logger.level
    if quiet:
        logger.setLevel(logging.CRITICAL)

    if os.path.isfile(pjoin(root_path,graph_name,'%s_pylink.so'%graph_name)):
        logger.info("%s: Generation already performed."%graph_name)
        if quiet:
            logger.setLevel(old_level)
        return

    def draw_graph():
        internal_lines = {drawing_input_internal_lines:s}
        power_list = {drawing_input_power_list:s}
        external_lines = {drawing_input_external_lines:s}
        extension = 'pdf'
        Gstart = 0 # If output looks bad, try changing this value to a non-negative integer
        neato = "neato" # shell command line to 'neato'
        drawing_filename = pjoin(root_path,graph_name,'%s_drawing'%graph_name)
        draw.plot_diagram(internal_lines, external_lines, drawing_filename, 
            powerlist=power_list, neato=neato, extension=extension, Gstart=0)

    logger.info("%s: Defining pySecDec loop integral..."%graph_name)
    propagators = {propagators:s}
    loop_momenta = {loop_momenta:s}
    external_momenta= {external_momenta:s}
    Lorentz_indices = {lorentz_indices:s}
    powerlist = {power_list:s}
    numerator = []
    with open("{numerator_path:s}",'r') as f:
        for line in f.readlines():
            numerator.append(line.strip())
    numerator = ' '.join(numerator)

    replacement_rules = {replacement_rules:s}
    with Silence(active=quiet):
        loop_integral = LoopIntegralFromPropagators(
            propagators, loop_momenta,
            external_momenta = external_momenta,   
            Lorentz_indices = Lorentz_indices,
            numerator = numerator,
            metric_tensor = 'g',
            dimensionality='4-2*eps',
            powerlist=powerlist,
            replacement_rules=replacement_rules
        )

    real_parameters = {real_parameters:s}
    n_loops = {n_loops:s}
    loop_additional_prefactor = {loop_additional_prefactor:s}
    # Note that the factor above is essentially the conversion from pySecDec normalisation conventions to those of alphaLoop, i.e.
    #loop_additional_prefactor = '( (I*(pi**((4-2*eps)/2)))/((2*pi)**(4-2*eps)) )**(%d)'%n_loops
    # Add the correcting MSbar prefactor without the renormalisation scale dependence. 
    # WARNING: this does not automatically give MSbar renormalisation because of nested subdivergences.
    #MSbar_factor = '(1/(4*pi*exp(-EulerGamma)))**eps'
    #loop_additional_prefactor = '((%s)*%s)'%(MSbar_factor,loop_additional_prefactor)

    logger.info("%s: Generating pySecDec package..."%graph_name)
    with Silence(active=quiet):
        loop_package(
            name = graph_name,
            loop_integral = loop_integral,
            real_parameters = real_parameters,
            complex_parameters = [],
            # the highest order of the final epsilon expansion --> change this value to whatever you think is appropriate
            requested_orders = [{max_epsilon_order:d}],
            contour_deformation = {contour_deformation:s},
            # the optimization level to use in FORM (can be 0, 1, 2, 3, 4)
            form_optimization_level = 2,
            form_threads = 2,
            processes = 16, # parallelisation
            # the WorkSpace parameter for FORM
            form_work_space = '100M',
            # the method to be used for the sector decomposition
            # valid values are ``iterative`` or ``geometric`` or ``geometric_ku``
            decomposition_method = 'iterative',
            # if you choose ``geometric[_ku]`` and 'normaliz' is not in your
            # $PATH, you can set the path to the 'normaliz' command-line
            # executable here
            #normaliz_executable='/path/to/normaliz',
            additional_prefactor = loop_additional_prefactor
        )

    logger.info("%s: Drawing graph..."%graph_name)
    draw_graph()

    logger.info("%s: Generating successfully completed in %.2f h."%(graph_name, (time.time()-start_time)/3600. ))

    start_time = time.time()
    logger.info("%s: Compiling package..."%graph_name)
    cmd = ['make', '-j', '%d'%multiprocessing.cpu_count()]
    if quiet:
        subprocess.run(cmd, cwd=pjoin(root_path,graph_name), shell=True, capture_output=False)
    else:
        subprocess.call(cmd, cwd=pjoin(root_path,graph_name), shell=True)

    if os.path.isfile(pjoin(root_path,graph_name,'%s_pylink.so'%graph_name)):
        logger.info("%s: Compilation successfully completed in %.2f h."%(graph_name, (time.time()-start_time)/3600. ))
    else:
        logger.info("%s: Compilation failed."%(graph_name))

    if quiet:
        logger.setLevel(old_level)

def run(graph_name, p_externals, quiet=True, eps_rel = 0.01, max_eval=1000000, min_n=10000, no_picobarns=False, integrator='vegas'):

    start_time = time.time()
    old_level = logger.level
    if quiet:
        logger.setLevel(logging.CRITICAL)

    if not os.path.isfile(pjoin(root_path,graph_name,'%s_pylink.so'%graph_name)):
        logger.critical("%s : Cannot integrate seems pySecDec compiled package is not found at '%s'."%(graph_name, pjoin(root_path,graph_name,'%s_pylink.so'%graph_name)))
        return

    logger.info("%s: Integrating with pySecDec..."%graph_name)

    integral_integrator = IntegralLibrary(pjoin(root_path,graph_name,'%s_pylink.so'%graph_name))

    if integrator == 'vegas':
        integral_integrator.use_Vegas(flags=2, epsrel=eps_rel, epsabs=1e-99, maxeval=max_eval) # ``flags=2``: verbose --> see Cuba manual
    elif integrator == 'qmc':
        integral_integrator.use_Qmc(transform='Korobov3', epsrel=eps_rel, epsabs=1e-99, maxeval=max_eval, minn=min_n)
    else:
        logger.critical("Integrator not recognized: %s"%integrator)

    def dot(v1,v2):
        return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3]
    complex_parameters_input = {complex_parameters_input:s}
    real_parameters_input = {real_parameters_input:s}    
    for i in range(len(p_externals)):
        for j in range(i, len(p_externals)):
            real_parameters_input.append(dot(p_externals[i],p_externals[j]))

    str_integral_without_prefactor, str_prefactor, str_integral_with_prefactor = integral_integrator(
        real_parameters=real_parameters_input, complex_parameters = complex_parameters_input
    )

    str_integral_without_prefactor = str_integral_without_prefactor.replace(',','+I*')
    str_prefactor = str_prefactor.replace(',','+I*')
    str_integral_with_prefactor = str_integral_with_prefactor.replace(',','+I*')

    integral_with_prefactor = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    integral_with_prefactor_err = sp.sympify(str_integral_with_prefactor.replace('+/-','*value+error*'))
    prefactor = sp.sympify(str_prefactor)
    integral_without_prefactor = sp.sympify(str_integral_without_prefactor.replace('+/-','*value+error*'))
    integral_without_prefactor_err = sp.sympify(str_integral_without_prefactor.replace('+/','*value+error*'))
    
    couplings_prefactor = "{couplings_prefactor:s}"
    couplings_values = {couplings_values:s}
    couplings_prefactor_eval = complex(eval(couplings_prefactor, couplings_values))

    if not no_picobarns and len(p_externals)>1:
        couplings_prefactor_eval *= 0.389379304e9

    # Include cross-section and normalisation correction factor
    n_loops = {n_loops:s}
    # Note that this additional_overall_factor is *already* included in the numerator given to FORM, 
    # so it does not need to be added here to couplings_prefactor_eval in order to recover the physical prediction
    additional_overall_factor = {additional_overall_factor:s}
    ptot = [sum(p[i] for p in p_externals) for i in range(4)]
    Ecm = math.sqrt(dot(ptot,ptot))
    if len(p_externals)>1:
        couplings_prefactor_eval *= (1j)*(2 / Ecm)**len(p_externals)

    str_res = []
    complex_res = dict([])
    def cmplx_str(c):
        return '%s%.16e %s %.16e*I'%('+' if c.real>=0 else '-', abs(c.real), '+' if c.imag>=0 else '-', abs(c.imag))
    for pole in range(-10,{max_epsilon_order:d}+1,1):
        if integral_with_prefactor.coeff('eps',pole).coeff('value')==0:
            continue
        # Remove additional_overall_factor as it contains helicity averaging factors that are not relevant for the raw answer
        complex_res[pole] = (
            complex(integral_with_prefactor.coeff('eps',pole).coeff('value'))/additional_overall_factor,
            complex(integral_with_prefactor.coeff('eps',pole).coeff('error'))/additional_overall_factor
        )
        str_res.append(' > %-9s: %-60s +/- (%-60s)'%(
            'eps^(%d) '%pole,
            cmplx_str(complex(couplings_prefactor_eval*integral_with_prefactor.coeff('eps',pole).coeff('value'))),
            cmplx_str(complex(couplings_prefactor_eval*integral_with_prefactor_err.coeff('eps',pole).coeff('error')))
        ))
        # print('Numerical Result')
        # print('eps^-2:', integral_with_prefactor.coeff('eps',-2).coeff('value'), '+/- (',integral_with_prefactor_err.coeff('eps',-2).coeff('error'), ')')
        # print('eps^-1:', integral_with_prefactor.coeff('eps',-1).coeff('value'), '+/- (',integral_with_prefactor_err.coeff('eps',-1).coeff('error'), ')')
        # print('eps^0 :', integral_with_prefactor.coeff('eps',0).coeff('value'), '+/- (',integral_with_prefactor_err.coeff('eps',0).coeff('error'), ')')

    logger.info("%s: Integral result:\n%s"%(graph_name,'\n'.join(str_res)))

    logger.info("%s: integration completed in %.0f h."%(graph_name, (time.time()-start_time)/3600.))
    
    if quiet:
        logger.setLevel(old_level)

    return complex_res

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""pySecDec loop generation and integrations script""")
    requiredNamed = parser.add_argument_group('required named arguments')
    # Required options
    #requiredNamed.add_argument(...)
    # Optional options
    parser.add_argument('--mode', '-m', dest='mode', type=str, default=None, 
        choices=('generate','run','generate_and_run'), help='Specify the run mode of this script (default: %(default)s)')

    parser.add_argument('--dump', '-d', dest='dump', type=str, default=None,
                        help='Dump Python-formatted result into a file named after this option. (default: no dump).')

    parser.add_argument('--integrator', '-i', dest='integrator', type=str, default='vegas', 
        choices=('vegas','qmc'), help='Specify which integrator to pick (default: %(default)s)')

    parser.add_argument('--max_eval', '-n', dest='max_eval', type=float, default=1.0e9, 
        help='Maximum number of evaluations for the integration (default: %(default)s)')

    parser.add_argument('--min_n', '-mn', dest='min_n', type=float, default=1.0e4, 
        help='Minimum number of evaluations for the integration with QMC (default: %(default)s)')

    parser.add_argument('--eps_rel', '-r', dest='eps_rel', type=float, default=0.01, 
        help='Target relative accuracy (default: %(default)s)')

    parser.add_argument(
        "-c", "--clean", action="store_true", dest="clean", default=False,
        help="Remove existing output.")

    parser.add_argument(
        "-np", "--no_picobarns", action="store_true", dest="no_picobarns", default=False,
        help="Never include picobarn conversion.")

    parser.add_argument(
        "-q", "--quiet", action="store_true", dest="quiet", default=False,
        help="Mute pySecDec output.")

    parser.add_argument(
        "-e", "--externals", dest="externals", type=str, default="{default_externals:s}",
        help="Specify numerical values of externals as a string of a Python-formatted list of tuples of four elements. (default: externals onshell and at rest if massive : %(default)s).")

    # parser.add_argument('--mode', '-hf', type=str, metavar='history_file',
    #                     default='rebalancing_history.dat',
    #                     help='Specify the history file to analysis')

    args = parser.parse_args()

    if not args.quiet:
        logger.info("You are running the pySecDec loop generation and integrations script.")

    try:
        parsed_externals = eval(args.externals)
        args.externals = parsed_externals
    except Exception as e:
        logger.critical("The numerical value specified for external momenta ('%s')cannot be parsed by Python. Error: %s"%(args.externals,str(e)))

    graph_name = "{graph_name:s}"

    if args.clean and os.path.isdir(pjoin(root_path,graph_name)):
        shutil.rmtree(pjoin(root_path,graph_name))

    if args.mode is None:
        if not os.path.isdir(pjoin(root_path,graph_name)):
            args.mode = 'generate_and_run'
        else:
            args.mode = 'run'

    final_res = None
    
    if args.mode in ['generate','generate_and_run']:

        generate(graph_name, quiet=args.quiet)
    
    if args.mode in ['run','generate_and_run']:

        final_res = run(graph_name, args.externals, quiet=args.quiet, eps_rel=args.eps_rel, max_eval=int(args.max_eval), no_picobarns=args.no_picobarns, integrator=args.integrator)

    if args.dump is not None:
        if not args.quiet: logger.info("Dumping final result into file '%s'."%args.dump)
        with open(args.dump, 'w') as f:
            f.write(pformat(final_res))