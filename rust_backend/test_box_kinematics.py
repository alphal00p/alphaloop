import ruamel.yaml
from ruamel.yaml.comments import CommentedSeq, CommentedMap
from copy import copy, deepcopy
import numpy as np
import sys
import subprocess
import os
import pandas as pd
import qcdloop as ql
pjoin = os.path.join


def sp(p, q):
    return p[0]*q[0]-sum(p[1:]*q[1:])


def square(p):
    return sp(p, p)


def get_sij(momenta, i, j):
    if i == j:
        return sp(momenta[i-1], momenta[i-1])
    else:
        mom = momenta[i-1]+momenta[j-1]
        return sp(mom, mom)


class Box(object):

    def __init__(self, externals, evaluations=250000000, cores=4, on_shell=False):
        self.on_shell = on_shell
        self.externals = externals
        self.new_externals = deepcopy(externals)

        # Rust
        self.config_file_name = 'config.yaml'
        self.new_config_file_name = 'tmpconfig.yaml'
        self.rust_executable_path = 'target/release/integrator'
        self.get_yaml()  # (self.config_file_name)
        self.box_instances = self.config['topologies']['box']
        self.config['max_eval'] = evaluations
        self.config['cores'] = cores

        # QCDLoop
        self.m = [0.0]*4
        self.mu2 = 1.7**2
        self.QL = ql.QCDLoop()

    def get_yaml(self):
        yaml = ruamel.yaml.YAML()
        yaml.indent(mapping=4, sequence=2, offset=4)
        with open(self.config_file_name) as fp:
            self.config = yaml.load(open(self.config_file_name))

    def store_yaml(self):
        yaml = ruamel.yaml.YAML()
        yaml.indent(mapping=4, sequence=2, offset=4)
        with open(self.new_config_file_name, 'w') as fp:
            yaml.dump(self.config, fp)

    def update_box(self):
        # Write new momenta with the usual format to the yaml file
        for i in range(4):
            self.box_instances['external_momenta'][i] = momentum = CommentedSeq(
                self.new_externals[i].tolist())
            momentum.fa.set_flow_style()
        # Update e_cm_sq
        e_cm_sq = square(self.new_externals[0] + self.new_externals[1])
        self.box_instances['e_cm_sq'] = float(e_cm_sq)

    ':: Approach different PS regions ::'

    def approach_soft(self, factor, typ=24):

        if self.on_shell:
            print("No on_shell routine for approach_soft()")
            exit()

        if factor <= 0:
            update_box(initial_momenta)
        else:
            if typ == 24:  # p2 and p4 going to zero
                self.new_externals = deepcopy(self.externals)
                self.new_externals[1] = box_momenta[1]*float(factor)
                self.new_externals[3] = box_momenta[3]*float(factor)
                self.new_externals[2] = self.new_externals[2] - \
                    sum(self.new_externals)
            elif typ == 23:  # p2 and p3 going to zero
                self.new_externals = deepcopy(self.externals)
                self.new_externals[1] = box_momenta[1]*float(factor)
                self.new_externals[2] = box_momenta[2]*float(factor)
                self.new_externals[3] = self.new_externals[3] - \
                    sum(self.new_externals)
            else:
                print("Unknwn option for approach_soft()")
                exit()

        self.update_box()

    def approach_coll(self, factor, on_shell=False):
        angle = np.pi*(float(factor))
        if self.on_shell:
            self.box_instances['on_shell_flag'] = 0b1111

            self.new_externals[0] = [-0.5, -0.5, 0.0, 0.0]
            self.new_externals[1] = [-0.5, 0.5, 0.0, 0.0]
            self.new_externals[2] = [0.5, -0.5 * np.cos(angle),
                                     -0.5*np.sin(angle), 0.0]
            self.new_externals[3] = [0.5, 0.5 * np.cos(angle),
                                     0.5*np.sin(angle), 0.0]
        else:
            self.new_externals[0] = [-0.5, -0.25, 0.0, 0.0]
            self.new_externals[1] = [-0.5, 0.25, 0.0, 0.0]
            self.new_externals[2] = [0.5, -0.25 * np.cos(angle),
                                     -0.25*np.sin(angle), 0.0]
            self.new_externals[3] = [0.5, 0.25 * np.cos(angle),
                                     0.25*np.sin(angle), 0.0]

        self.update_box()

    ':: The analytic results are computed with qcdloop ::'

    def analytic(self):
        if self.on_shell:
            return self.analytic_sub()
        else:
            return self.analytic_box()

    def analytic_box(self):
        s11 = get_sij(self.new_externals, 1, 1)
        s22 = get_sij(self.new_externals, 2, 2)
        s33 = get_sij(self.new_externals, 3, 3)
        s44 = get_sij(self.new_externals, 4, 4)
        s12 = get_sij(self.new_externals, 1, 2)
        s23 = get_sij(self.new_externals, 2, 3)
        p = [s11, s22, s33, s44, s12, s23]
        return self.QL.integral(self.mu2, self.m, p)[0].real

    def analytic_sub(self):
        s11 = get_sij(self.new_externals, 1, 1)
        s22 = get_sij(self.new_externals, 2, 2)
        s33 = get_sij(self.new_externals, 3, 3)
        s44 = get_sij(self.new_externals, 4, 4)
        s12 = get_sij(self.new_externals, 1, 2)
        s23 = get_sij(self.new_externals, 2, 3)

        p = [s11, s22, s33, s44, s12, s23]
        box_eval = np.array(self.QL.integral(self.mu2, self.m, p))
        tri_s_eval = np.array(self.QL.integral(
            self.mu2, [0.0, 0.0, 0.0], [0.0, 0.0, s12]))
        tri_t_eval = np.array(self.QL.integral(
            self.mu2, [0.0, 0.0, 0.0], [0.0, 0.0, s23]))

        return (box_eval - 2.0*tri_s_eval/s23 - 2.0*tri_t_eval/s12)[0].real

    ':: Numerical integratoin done calling rust ::'

    def run_rust(self):
        """ Steer a run of rust integration with the current yaml config setup and return the result."""

        # Remove some temporary files if necessary
        try:
            os.remove('box_state.dat')
        except:
            pass

        # create new config file
        self.store_yaml()

        # Run rust
        subprocess.call(' '.join([
            self.rust_executable_path, '-t=box', '--config=%s' % self.new_config_file_name,
        ]), shell=True)
        # Now we should be able to return the result of the run
        try:
            yaml = ruamel.yaml.YAML()
            result = yaml.load(open('box_res.dat', 'r'))
            os.remove(pjoin(os.getcwd(), 'box_res.dat'))
            return result
        except Exception as e:
            raise BaseException("ERROR: Could not load rust integration result from file: %s. Error: %s" % (
                'box_res.dat',
                str(e)
            ))


import yaml

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("The function needs one arguments: [coll/soft]")
        exit()

    coll_or_soft = sys.argv[1]

   # box_momenta is used by approach_soft()
    box_momenta = np.array([[-0.47809952420694085, 0, 0, -0.4675244967345596],
                            [-0.5085067895779792, 0, 0, 0.4675244967345596],
                            [0.45782801395958195, 0.13758384614384497,
                             0.0812175730388203, -0.3067260691172595],
                            [0.5287782998253381, -0.13758384614384497, -0.0812175730388203, 0.3067260691172595]])

    # Create Box object
    box = Box(box_momenta, evaluations=1000000, cores=4, on_shell=False)

    factors = []
    analytic_results = []
    rust_results = []
    errors = []
    probs = []

    for factor in np.linspace(0.0001, 1.0, num=40):
        # Create new PS point
        if coll_or_soft == 'coll':
            print("Collinear case with factor %f" % factor)
            box.approach_coll(factor)
        elif coll_or_soft == 'soft':
            print("Soft case with factor %f" % factor)
            box.approach_soft(factor)
        else:
            print("Unknown option %s" % coll_or_soft)
            exit()

        # Evaluate
        rust_result = box.run_rust()
        qcdloop_result = box.analytic()

        # Store
        analytic_results += [qcdloop_result]
        rust_results += [float(rust_result['result'][0])]
        errors += [float(rust_result['error'][0])]
        probs += [float(rust_result['prob'][0])]
        factors += [factor]

    # Save in csv file
    panda_results = pd.DataFrame(
        {'factor': factors, 'analytic': analytic_results, 'result': rust_results, 'error': errors, 'prob': probs})
    panda_results.to_csv('%s.csv' % coll_or_soft,
                         sep='\t',
                         index=False,
                         columns=['factor', 'analytic', 'result', 'error', 'prob'])
