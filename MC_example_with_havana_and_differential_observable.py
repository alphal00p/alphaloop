#!/usr/bin/env python3
import ltd
import time
import math
import os
from pprint import pprint, pformat
from ltd.havana import Havana, Sample, GridConstructor, ContinuousGridConstructor, DiscreteGridConstructor

_BATCH_SIZE = 100000
_SEED = 1337
_N_ITERATIONS = 10
_N_OBSERVABLE_BINS_PER_DIM = 100
_N_GRID_BINS = 128
_N_MIN_POINTS = 10000
_N_DIMS = 3
_LEARNING_RATE = 1.5
_BIN_REPARTITIONING_SCHEDULE = [1,1,2,1,1,2,1,1,2]
_TRAIN_ON_AVG = False
_OUTPUT_DIR = 'MC_example_output'
_START_GRID = None # Specify a grid to load from, e.g 'havana_grid_iteration_003.yaml'

class Bin(object):

    def __init__(self):
        self.cum_wgt = 0.
        self.cum_wgt_sqr = 0.
        self.n_points = 0
        self.estimates_per_iterations = []
        self.current_estimate = 0.
        self.current_error = 0.

    def accumulate(self, wgt):

        self.cum_wgt += wgt
        self.cum_wgt_sqr += wgt**2
        self.n_points += 1

    def update(self):

        if self.n_points > 0:
            self.estimates_per_iterations.append(
                (
                    self.cum_wgt / self.n_points,
                    math.sqrt( ((self.cum_wgt_sqr / self.n_points) - (self.cum_wgt / self.n_points )**2) / self.n_points ),
                    self.n_points
                )
            )
        else:
            self.estimates_per_iterations.append( (0., 0.) )
        
        # Reset the cumulative quantities now
        self.cum_wgt = 0.
        self.cum_wgt_sqr = 0.
        self.n_points = 0

        if all(estimate[1]==0. for estimate in self.estimates_per_iterations):
            self.current_estimate = 0.
            self.current_error = 0.
            return
    
        # Compute current best estimate by accumulating past results including a weighting according to inverse errors.
        sum_inverse_variances = sum(1./estimate[1] for estimate in self.estimates_per_iterations if estimate[1]!=0.)
        self.current_estimate = sum( estimate[0]/estimate[1] for estimate in self.estimates_per_iterations if estimate[1]!=0. ) / sum_inverse_variances
        
        sum_inverse_variances_squared = sum(1./estimate[1]**2 for estimate in self.estimates_per_iterations if estimate[1]!=0.)
        self.current_error = math.sqrt( sum( 1 for estimate in self.estimates_per_iterations if estimate[1]!=0. ) / sum_inverse_variances_squared )

class Histogram2D(list):

    # We simplify the implementation of this class by assuming input normalised within [0, 1]

    def __init__(self, n_bins, **opts):
        self.n_bins = n_bins
        self[:] = [[ Bin() for _ in range(n_bins) ] for _ in range(n_bins) ]

    def fill(self, x,y,wgt):

        a_bin = self[ math.floor(x*self.n_bins) ][ math.floor(y*self.n_bins) ]
        a_bin.accumulate(wgt)

    def update(self):

        for row in self:
            for a_bin in row:
                a_bin.update()

    def draw(self, filename, format='png', title=None):

        if format=='raw':
            with open(filename, 'w') as f:
                f.write(pformat(
                    [ [ (a_bin.current_estimate, a_bin.current_error) for a_bin in row ] for row in self]
                ))
            return

        import numpy as np
        import matplotlib.pyplot as plt
        import random
        

        x_bins = np.array([i/float(self.n_bins) for i in range(self.n_bins+1)])
        y_bins = np.array([i/float(self.n_bins) for i in range(self.n_bins+1)])

        x = np.array(sum([[coord,]*self.n_bins for coord in x_bins[:-1]], []))
        y = np.array(sum( [list(y_bins[:-1]),]*self.n_bins ,[]))

        fig, ax = plt.subplots(figsize =(10, 7))

        weights = np.array(sum( [ [a_bin.current_estimate for a_bin in row] for row in self ],[]))

        # Creating plot
        plt.hist2d(x, y, [x_bins, y_bins], weights=weights)
        if title:
            plt.title(title)
        
        ax.set_xlabel('x_1')
        ax.set_ylabel('x_2')

        # show plot
        plt.tight_layout()
        #plt.show()
        plt.savefig(filename)

class Observable(dict):

    def __init__(self, *args, **opts):

        # Potentially add additional histograms if so enclined, in this example we pick one only
        # This is what this layer is useful, to accommodate more historgrams

        self.my_histogram2D = Histogram2D(_N_OBSERVABLE_BINS_PER_DIM)

    def fill(self, xs, wgt):
        self.my_histogram2D.fill(xs[0], xs[1], wgt)

    def update(self, *args, **opts):
        self.my_histogram2D.update(*args, **opts)

    def draw(self, *args, **opts):
        self.my_histogram2D.draw(*args, **opts)

def my_function(xs, havana_weight, observable):

    integrand_wgt = math.exp(-sum((x/(1.-x))**2 for x in xs))
    for x in xs:
        integrand_wgt *= 1./(1.-x)**2

    observable.fill(xs, integrand_wgt*havana_weight)
    
    return integrand_wgt

def main():
    
    if not os.path.exists(_OUTPUT_DIR):
        os.mkdir(_OUTPUT_DIR)

    n_grid_bins = _N_GRID_BINS

    if _START_GRID is None:
        grid = GridConstructor(cgc=ContinuousGridConstructor(_N_DIMS, n_grid_bins, _N_MIN_POINTS) )
        havana_sampler = Havana(grid, seed=_SEED)
    else:
        # One must specify the format, either 'bin' or 'yaml', which we take here from the file extension
        with open(_START_GRID,'rb') as f:
            havana_sampler = Havana.load_grid( f.read(), seed=_SEED, format=_START_GRID.split('.')[-1] )

    target_result = (math.sqrt(math.pi)/2.)**_N_DIMS

    my_observable = Observable()

    for i_iteration in range(_N_ITERATIONS):

        havana_sampler.sample(_BATCH_SIZE)
        function_evaluations = []

        for sample in havana_sampler.get_samples():

            havana_weight, xs = sample.continuous_sample[0]
            function_evaluations.append(
                my_function(
                    xs, havana_weight, my_observable
                )
            )
        my_observable.update()

        havana_sampler.add_training_samples(function_evaluations)

        n_grid_bins *= (_BIN_REPARTITIONING_SCHEDULE[i_iteration] if len(_BIN_REPARTITIONING_SCHEDULE)>i_iteration else 1) 

        havana_sampler.update(
            alpha=_LEARNING_RATE, 
            new_bin_length=n_grid_bins,
            train_on_avg=_TRAIN_ON_AVG
        )
        avg, err, chi_sq, max_eval_negative, max_eval_positive, n_evals, n_zero_evals = havana_sampler.get_current_estimate()

        print("Result a iteration %-2d after %-8d evaluations: %.5g +/- %.5g (%.2g%%)"%(
                        i_iteration , n_evals, avg, err, 0. if avg==0. else (abs(err)/avg)*100.))
        print( "vs target                                       : %.5g del %.5g (%.2g%%)"%(
            target_result,avg-target_result,  0. if target_result==0. else ((avg-target_result)/target_result)*100.0
        ))

        # Dump the sampling grid obtained at the end of this iteration
        bytesvec = havana_sampler.dump_grid(format='yaml')
        with open(os.path.join(_OUTPUT_DIR,'havana_grid_iteration_%03d.yaml'%(i_iteration+1)),'wb') as f:
            f.write(bytesvec)

        # Dump the histograms obtained from this iteration
        my_observable.draw(os.path.join(_OUTPUT_DIR,'my_observable_iteration_%03d.py'%(i_iteration+1)), format='raw')
        my_observable.draw(os.path.join(_OUTPUT_DIR,'my_observable_iteration_%03d.png'%(i_iteration+1)), format='png', 
                            title='My differential observable at iteration %d'%(i_iteration+1))

    print("Integration complete")
    print("You can run 'ffmpeg -i my_observable_iteration_%%03d.png my_observable.gif' in %s to animate your observable plot evolution across iterations."%_OUTPUT_DIR)

if __name__ == '__main__':
    main()