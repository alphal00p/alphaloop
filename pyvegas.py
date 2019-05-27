#!/usr/bin/env python2
import multiprocessing
import numpy as np
import vegas
import ltd
import signal

class parallelintegrand(vegas.BatchIntegrand):
    """ Convert (batch) integrand into multiprocessor integrand.
    Integrand should return a numpy array.
    """
    def __init__(self, integrand, fcn, nproc=multiprocessing.cpu_count()):
        " Save integrand; create pool of nproc processes. "
        self.integrand = integrand
        self.fcn = fcn
        self.nproc = nproc

        self.q_in = multiprocessing.Queue()
        self.q_out = multiprocessing.Queue()

        proc = [multiprocessing.Process(target=fcn, args=(integrand, self.q_in, self.q_out)) for _ in range(nproc)]
        for p in proc:
            p.daemon = True
            p.start()

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        self.pool = multiprocessing.Pool(processes=nproc)
        signal.signal(signal.SIGINT, original_sigint_handler)
    def __del__(self):
        " Standard cleanup. "
        self.pool.close()
        self.pool.join()
    def __call__(self, x):
        " Divide x into self.nproc chunks, feeding one to each process. "
        nx = x.shape[0] // self.nproc + 1

        #print('x call', x)

        # distribute the input over the cores and send them to the queue with a label
        # so we can reconstruct the proper order later
        worker_id = 0
        for i in range(self.nproc):
            xr = x[i*nx : (i+1)*nx]
            if len(xr) != 0:
                self.q_in.put((worker_id, xr))
                worker_id += 1

        res = [[]] * worker_id
        for _ in range(worker_id):
            (index, result) = self.q_out.get()
            res[index] = result

        return np.concatenate(res)


def integrand_bridge(integrand, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break

        # FIXME: phase both is not supported yet
        ans = np.zeros(x.shape[0], np.double)   
        for y in range(x.shape[0]):
            ans[y] = integrand.evaluate(x[y])[1] # FIXME: hardcoded imaginary part!
        q_out.put((i, ans))

def main():
    # seed the random number generator so results reproducible
    np.random.seed((1, 2, 3))

    l1 = ltd.LTD("LTD/topologies.yaml", "manual_Box_no_ellipse", "LTD/hyperparameters.yaml")

    integ = vegas.Integrator(3 * [[0, 1]])
    # adapt
    fparallel = parallelintegrand(l1, integrand_bridge)
    integ(fparallel, nitn=10, neval=1e6)
    # final results
    result = integ(fparallel, nitn=10, neval=1e6)

    print('result = %s    Q = %.2f' % (result, result.Q))

if __name__ == '__main__':
    main()

