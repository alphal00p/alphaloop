#!/usr/bin/env python2
import multiprocessing
import numpy as np
import vegas
import ltd
import signal
import argparse
import yaml
import math

class Integrand(vegas.BatchIntegrand):
    def __init__(self, integrand, n_loops, phase, nproc):
        self.integrand = integrand
        self.n_loops = n_loops
        self.q = []
        self.phase = 0 if phase == "real" else 1
        self.nproc = multiprocessing.cpu_count() if nproc is None else nproc

        self.q_in = multiprocessing.Queue()
        self.q_out = multiprocessing.Queue()

        self.proc = [multiprocessing.Process(target=self.integrand_bridge, args=(
            integrand, self.q_in, self.q_out)) for _ in range(self.nproc)]
        for p in self.proc:
            p.daemon = True
            p.start()

    def integrand_bridge(self, integrand, q_in, q_out, f128=False):
        signal.signal(signal.SIGINT, signal.SIG_IGN)

        loop_momenta = [np.zeros(3, np.double) for _ in range(self.n_loops)]
        while True:
            i, x = q_in.get()
            if i is None:
                break

            # FIXME: phase both is not supported yet
            ans = np.zeros(x.shape[0], np.double)
            for y in range(x.shape[0]):
                if f128:
                    # do the parameterization
                    jac = 1.
                    for l in range(self.n_loops):
                        r = integrand.parameterize(x[y][l*3:(l+1)*3], l, 1.0) # FIXME: e_cm_squared set to 1.0
                        loop_momenta[l] = r[:3]
                        jac *= r[3]

                    ans[y] = integrand.evaluate_f128(loop_momenta)[self.phase] * jac
                else:
                    ans[y] = integrand.evaluate_integrand(x[y])[self.phase]

                if math.isnan(ans[y]):
                    print(loop_momenta, ans[y])
            q_out.put((i, ans))

    def __call__(self, x):
        nx = x.shape[0] // self.nproc + 1

        # distribute the input over the cores and send them to the queue with a label
        # so we can reconstruct the proper order later
        worker_id = 0
        for i in range(self.nproc):
            xr = x[i*nx: (i+1)*nx]
            if len(xr) != 0:
                self.q_in.put((worker_id, xr))
                worker_id += 1

        res = [[]] * worker_id
        for _ in range(worker_id):
            (index, result) = self.q_out.get()
            res[index] = result

        return np.concatenate(res)


def main():
    # seed the random number generator so results reproducible
    np.random.seed((1, 2, 3))

    parser = argparse.ArgumentParser(description='Integrate LTD topologies using pyVegas',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('loops', default="1", type=int,
                        help='Number of loops of the topology')
    parser.add_argument('--phase', default="real", type=str,
                        choices=['real', 'imag'], help='The phase')
    parser.add_argument('-c', default=None, type=int,
                        help='Number of cores. Defaults to maximum available')
    parser.add_argument('-hf', default="LTD/hyperparameters.yaml", type=str,
                        help='Hyperparameters file')
    parser.add_argument('--cross_section', default="LTD/bu_squared.yaml", type=str,
                        help='Cross section file')
    parser.add_argument('--survey_iter', '--si', default=5, type=int,
                        help='number of survey iterations')
    parser.add_argument('--survey_neval', '--sn', default=int(1e6), type=int,
                        help='number of survery evaluations')
    parser.add_argument('--refine_iter', '--ri', default=10, type=int,
                        help='number of refine iterations')
    parser.add_argument('--refine_neval', '--rn', default=int(1e7), type=int,
                        help='number of refine evaluations')
    parser.add_argument('--out', default=None, type=str,
                        help='output file')

    args = parser.parse_args()

    l1 = ltd.CrossSection(args.cross_section, args.hf)

    print('Integrating %s with n_loops=%s %s iterations survey with eval=%s, %s iterations refine eval=%s' % (
        args.cross_section, args.loops, args.survey_iter, args.survey_neval, args.refine_iter, args.refine_neval))

    integ = vegas.Integrator(3 * args.loops * [[0, 1]])
    # adapt
    fparallel = Integrand(l1, args.loops, args.phase, args.c)
    result = integ(fparallel, nitn=args.survey_iter, neval=args.survey_neval)
    print('Done adapting: %s    Q = %.2f' % (result, result.Q))
    # final results
    result = integ(fparallel, nitn=args.refine_iter, neval=args.refine_neval)

    print(result.summary())
    print('result = %s    Q = %.2f' % (result, result.Q))

    if args.out is not None:
        data = {
            'neval': args.survey_iter * args.survey_neval + args.refine_iter * args.refine_neval,
            'fail': 1,
            'result': [result.mean],
            'error': [result.sdev],
            'prob': [result.Q]
        }

        with open(args.out, 'w') as outfile:
            yaml.dump(data, outfile, default_flow_style=False)
            outfile.write('...\n')


if __name__ == '__main__':
    main()
