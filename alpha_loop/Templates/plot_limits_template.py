#!/usr/bin/env python

import datetime
import sys
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from matplotlib.backends.backend_pdf import PdfPages

def main(output_path=None, SG_names=None, limits=None, use_latex=None, **opts):
     
    with PdfPages(output_path) as pdf:
    
        plt.rcParams['text.usetex'] = use_latex 
        plt.figure(figsize=(8, 6))
        
%(main_body)s

        d = pdf.infodict()
        d['Title'] = 'Profiling of Local Unitarity integrand limits'
        d['Author'] = 'alphaLoop'
        d['Subject'] = 'Visualisation of the approach of singular limits of Local Unitarity integrands.'
        d['Keywords'] = 'Local Unitarity Limits Singularity Profiling'
        d['CreationDate'] = datetime.datetime(2009, 11, 13)
        d['ModDate'] = datetime.datetime.today()     

%(plot_function_definitions)s

%(data_loading_definitions)s

plotter_parser = ArgumentParser(prog='plotter')
plotter_parser.add_argument('output_path', metavar='output_path', type=str,
    help='The full paths to write the plots into, including the .pdf extension.')
plotter_parser.add_argument("--limits", dest='limits', type=int, nargs='?', default=None,
    help='Specify a list of limits ID to render.')
plotter_parser.add_argument("--SG_name", dest='SG_names', type=str, nargs='?', default=None,
    help='Specify a list of SG names to render.')
plotter_parser.add_argument(
    "-wl", "--with_latex", action="store_true", dest="use_latex", default=False,
    help="Enable latex rendering. Prettier labels but slower generation then.")
plotter_parser.add_argument(
    "-sc", "--show_cuts", action="store_true", dest="show_cuts", default=False,
    help="Also show cuts in the final plots.")

if __name__ == '__main__':
    args = plotter_parser.parse_args()    
    main(**vars(args))
    print("Plots successfully generated in file '{}'".format(args.output_path))
