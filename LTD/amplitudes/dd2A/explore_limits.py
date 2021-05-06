#!/usr/bin/env python

import itertools
import math
import os
import sys

# colorful print
if True:
    from pprint import pformat

    from pygments import highlight
    from pygments.formatters import Terminal256Formatter
    from pygments.lexers import Python3Lexer
    from pygments.styles import get_all_styles

    # print(list(get_all_styles()))
    def cprint(out):
        print(color_out(out))

    def color_out(out):
        return highlight(pformat(out),
                         Python3Lexer(),
                         Terminal256Formatter(style='stata-dark'))[:-1]
        # Terminal256Formatter(style='inkpot')))
        # Terminal256Formatter(style='paraiso-dark')))
        # Terminal256Formatter(style='monokai')))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.font_manager import FontProperties
from scipy.optimize import curve_fit

pjoin = os.path.join

# Global paths
MG_PATH = "/home/andrea/Programs/MG5_aMC_v2_7_2_py3/"
AL_PATH = "/home/andrea/BitBucket/alphaloop/"
LIB_PATH = AL_PATH
LTD_PATH = pjoin(AL_PATH, "LTD")


try:
    file_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, pjoin(AL_PATH, "LTD"))
    sys.path.insert(0, AL_PATH)

    # Import Alpha loop bindings
    import ltd_commons
    import vectors
    # Import the rust bindings
    from ltd import CrossSection
except ImportError:
    raise BaseException(
        "Could not import the rust back-end 'ltd' module. Compile it first with:"
        " ./make_lib from within the pyNLoop directory.")


class CrossSectionWorker(object):
    def __init__(self, cross_section_folder, settings_file, SGs=None, cross_section_set=None):
        self.process_title = os.path.basename(
            os.path.normpath(cross_section_folder))
        if cross_section_set == None:
            self.get_SG_names(cross_section_folder)
        else:
            self.get_SG_names(cross_section_folder, cross_section_set)
        self.cross_section_files = {
            sg: pjoin(cross_section_folder, "%s.yaml" % sg) for sg in self.SG_names}
        self.settings_file = settings_file
        self.settings = yaml.load(open(settings_file, 'r'), Loader=yaml.Loader)
        self.cross_sections = {}
        self.workers = {}

        if SGs is None:
            self.SGs = self.SG_names
        else:
            self.SGs = SGs

        for sg in self.SGs:
            if not sg in self.SG_names:
                print("%s is not present in %s. Possible SGs:")
                cprint(self.SG_names)
                sys.exit(1)
            with open(pjoin(cross_section_folder, sg+'.yaml'), 'r') as f:
                squared_topology = yaml.load(f, Loader=yaml.Loader)
                self.cross_sections[sg] = squared_topology
                self.workers[sg] = CrossSection(
                    squared_topology_file=self.cross_section_files[sg],
                    settings_file=self.settings_file)

    @classmethod
    def from_MG_folder(class_object, process_name, settings_file=None, SGs=None):
        cross_section_folder = pjoin(MG_PATH, process_name, 'Rust_inputs')
        if settings_file is None:
            settings_file = pjoin(
                MG_PATH, process_name, 'run_workspace', 'hyperparameters.yaml')
        CSW = class_object(cross_section_folder=cross_section_folder,
                           settings_file=settings_file,
                           SGs=SGs)
        print("Initiate workers for each SuperGraph of process {}:\n{}".format(color_out(process_name),
                                                                               color_out(CSW.workers)))
        return class_object(cross_section_folder=cross_section_folder,
                            settings_file=settings_file,
                            SGs=SGs)

    def get_SG_names(self, cross_section_folder, cross_section_set='all_QG_supergraphs.yaml'):
        all_supergraphs = pjoin(cross_section_folder, cross_section_set)
        self.SG_names = []
        with open(pjoin(cross_section_folder, cross_section_set), 'r') as stream:
            self.cross_section_set = yaml.load(stream, Loader=yaml.Loader)
            for topo in self.cross_section_set['topologies']:
                self.SG_names.append(topo['name'])


class SuperGraphWorker(object):
    def __init__(self, csw: CrossSectionWorker, SG, process_title=None):
        if not SG in csw.SGs:
            raise ValueError("SG unkown.")
        if process_title is None:
            self.process_title = csw.process_title
        else:
            self.process_title = process_title

        self.SG = SG
        self.worker = csw.workers[SG]
        self.settings = csw.settings
        self.squared_topology = csw.cross_sections[SG]
        self.cross_section_set = csw.cross_section_set

    def get_loop_momenta(self):
        return self.squared_topology['loop_momentum_basis']

    def get_name(self):
        return self.process_title


#####################
# Utility functions
#####################


def chunks(long_list, chunk_size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(long_list), chunk_size):
        yield long_list[i:i + chunk_size]


def rel_diff(ef64, ef128):
    rel = []
    for (c1, c2) in zip(ef64, ef128):
        d = abs(c1-c2)
        if abs(c1) == 0 and abs(c2) == 0:
            rel += [0]
        elif abs(c1) == 0 and abs(c2) != 0:
            rel += [None]
        else:
            rel += [d/abs(c2)]
    return np.array(rel)


def rescale(loop_momenta, c, signature):
    return np.array([l*c if s == 1 else l for l, s in zip(loop_momenta, signature)])


def line(x, a, b):
    return a*x+b


def sp_spatial(mom1, mom2):
    return mom1[1]*mom2[1] + mom1[2]*mom2[2] + mom1[3]*mom2[3]


def make_onshell(mom):
    mom[0] = np.sqrt(sp_spatial(mom, mom))


def sp(mom1, mom2):
    return mom1[0]*mom2[0] - sp_spatial(mom1, mom2)


def sq(mom):
    return sp(mom, mom)


def collinear_limit(mom, delta, shift=np.array([0, 0, 0, 0]), x=1.0):
    p1 = mom
    p2 = mom * [1, -1, -1, -1]

    q = [0., 1.2, 12., -1.45, 0.03]  # vectors.LorentzVector(np.random.rand(4))
    l_perp = -np.array([0.,
                        p1[2]*q[3] - p1[3]*q[2],
                        p1[3]*q[1] - p1[1]*q[3],
                        p1[1]*q[2] - p1[2]*q[1]])

    coll_mom = (shift + x*p1 + delta*p2 + np.sqrt(delta)*l_perp)
    return coll_mom  # , rust_instance.inv_parameterize_f128(c_mom[1:], 0)[:3]


def check_instability(rust_instance, loop_momenta, min_log=-5, max_log=5, n_points=100, rescale_which=None):
    if rescale_which is None:
        rescale_which = [1]*len(loop_momenta)
    # Rescale loop momenta
    ts = np.logspace(min_log, max_log, n_points)
    ls = [rescale(loop_momenta, x, rescale_which) for x in ts]

    # Evaluate
    e_f64 = [complex(*rust_instance.evaluate(l)) for l in ls]
    e_f128 = [complex(*rust_instance.evaluate_f128(l)) for l in ls]

    # Parametrization Jacobian
    paramerization_jacobians = []
    for momenta in ls:
        jac = 1.0
        for n, l in enumerate(momenta):
            jac *= rust_instance.inv_parameterize_f128(l, n, 125**2)[3]
        paramerization_jacobians += [1/jac]

    # Relative Difference
    rdiff = rel_diff(e_f64, e_f128)

    return (ts, rdiff, e_f64, e_f128, paramerization_jacobians)


def fit_exp_curve(xs, ys, x_min=None, x_max=None):
    if x_min is None:
        x_min = min(xs)
    if x_max is None:
        x_max = max(xs)

    Xs = []
    Ys = []
    for x, y in zip(xs, ys):
        if x <= x_max and x >= x_min:
            Xs += [x]
            Ys += [y]
    return curve_fit(line, np.log(Xs), np.log(Ys))


#########################
# Setting Files
#########################
settings_file = pjoin(AL_PATH, 'LTD', 'hyperparameters.yaml')

#############################
# Import individual diagrams
#############################
qqbar2a_diags = []
diags_titles = ['D1', 'D2', 'D3', 'D4', 'IR',
                'D2UV', 'D3UV', 'D4UV', 'IRUV', 'D4UV2']
for diag_n in range(1, 10+1):
    process = "qqbar2a_%d" % diag_n
    cross_section_folder = pjoin(
        AL_PATH, 'LTD/amplitudes/dd2A', process, 'color_struc_0/Rust_inputs')
    cross_section_set = '%s_color_0.yaml' % process
    csw = CrossSectionWorker(cross_section_folder, settings_file=settings_file,
                             cross_section_set=cross_section_set)
    qqbar2a_diags += [SuperGraphWorker(csw, csw.SG_names[0],
                                       process_title=diags_titles[diag_n-1])]
    qqbar2a_diags[-1].worker.externals_from_set(
        pjoin(cross_section_folder, cross_section_set))


#########################
# Import process
#########################
process = "qqbar2a"
cross_section_folder = pjoin(
    AL_PATH, 'LTD/amplitudes/dd2A', process, 'color_struc_0/Rust_inputs')
cross_section_set = '%s_color_9.yaml' % process
csw = CrossSectionWorker(cross_section_folder, settings_file=settings_file,
                         cross_section_set=cross_section_set)

# Print list of supergraphs
print("\nSupergraphs {}:\n\t- {}".format(
    color_out(process), color_out(csw.SG_names)))

#############################
# Select SG - RUST INSTANCE
############################
qqbar2a = SuperGraphWorker(
    csw, '%s_color_9_0' % process,
    process_title='q q~ -> a a @ NLO')
qqbar2a.worker.externals_from_set(
    pjoin(cross_section_folder, cross_section_set))
# Print incoming momentum
print("\nSettings -> incoming momenta:\n\t- {}".format(
    color_out(qqbar2a.settings['CrossSection']['incoming_momenta'])))

# Print loop momenta
for sgw in [qqbar2a]:
    print("\nLoop Momentum Basis of {} {}:\n\t- {}".format(
        color_out(sgw.get_name()), sgw.SG, color_out(sgw.get_loop_momenta())))


######################################
# TEST Evaluation
######################################
np.random.seed(1)
#rust_instance = qqbar2a.worker
#external_data = qqbar2a.cross_section_set['external_data']
rust_instance = qqbar2a_diags[1].worker
external_data = qqbar2a_diags[1].cross_section_set['external_data']
n_loops = 1+1
#
# Print functions
#print("Callable functions for CrossSection:")
#cprint(list(filter(lambda x: not "__" in x, rust_instance.__dir__())))
# sys.exit(1)
#
xs = np.random.rand(n_loops*3)
ls = []
overall_jacobian = 1.
for n, x in enumerate(chunks(xs, 3)):
    out = rust_instance.parameterize(x, n, 125**2)
    ls.append([0.]+list(out[:3]))
    # make_onshell(ls[-1])
    overall_jacobian *= out[3]

ls = [[0.0, 0.008068658704146322,  -0.006698226849113513, 0.7561329592139684]]\
    + external_data['out_momenta']

scaling = list(filter(lambda x: x[0] > 0, rust_instance.get_scaling(ls, 0)))[0]
ls_scaled = [np.array(l)*scaling[0] for l in ls]
scaling2 = list(filter(lambda x: x[0] > 0,
                       rust_instance.get_scaling(ls_scaled, 0)))[0]
#scaling = (1, 1)
evaluate_cut = complex(*rust_instance.evaluate_cut(ls, 0, *scaling))
print("EVALUATE CUT:\t{}".format(color_out(evaluate_cut)))
evaluate = complex(*rust_instance.evaluate_f128(ls))
print("EVALUATE:\t{}".format(color_out(evaluate)))
#xs = []

e_cm_sq = sq(sum(np.array(x) for x in external_data['in_momenta']))
cprint(e_cm_sq)

# sys.exit(0)

print("Scaling {}".format(color_out(scaling)))
print("Scaling2 {}".format(color_out(scaling2)))
print("Scaled Momenta:\n{}".format(
    color_out([list(np.array(l)*scaling[0]) for l in ls])))
ls_scaled = [np.array(l)*scaling[0] for l in ls]
evaluate_integrand = complex(*rust_instance.evaluate_integrand(xs))
evaluate_cut = complex(*rust_instance.evaluate_cut(ls, 0, *scaling))
evaluate = complex(*rust_instance.evaluate_f128(ls))
print("EVALUATE INTEGRAND:\t{}".format(color_out(evaluate_integrand)))
print("EVALUATE CUT * JAC:\t{}".format(color_out(evaluate_cut*overall_jacobian)))
print("EVALUATE CUT:\t{}".format(color_out(evaluate_cut)))
print("EVALUATE:\t{}".format(color_out(evaluate)))


##############################
# EXPLORATION
##############################

to_explore = []
to_explore.append(qqbar2a)
# to_explore.append(
#    next(filter(lambda x: x.process_title == 'D1', qqbar2a_diags)))
# to_explore.append(
#    next(filter(lambda x: x.process_title == 'D2', qqbar2a_diags)))
# to_explore.append(
#    next(filter(lambda x: x.process_title == 'D3', qqbar2a_diags)))
to_explore.append(
    next(filter(lambda x: x.process_title == 'D4', qqbar2a_diags)))
# to_explore.append(
#    next(filter(lambda x: x.process_title == 'IR', qqbar2a_diags)))
to_explore.append(
    next(filter(lambda x: x.process_title == 'D4UV', qqbar2a_diags)))
to_explore.append(
    next(filter(lambda x: x.process_title == 'D4UV2', qqbar2a_diags)))
# to_explore.append(
#   next(filter(lambda x: x.process_title == 'IRUV', qqbar2a_diags)))
# to_explore.append(
#   next(filter(lambda x: x.process_title == 'IR', qqbar2a_diags)))
#to_explore = qqbar2a_diags

# Momentum q
p1 = np.array(qqbar2a.settings['CrossSection']['incoming_momenta'][0])
# Momentum qbar
p2 = np.array(qqbar2a.settings['CrossSection']['incoming_momenta'][1])

# loop_momenta
loops_momenta = [(np.random.rand(4)-0.5)*1e0 for _ in range(
    n_loops-len(external_data['out_momenta']))]+external_data['out_momenta']

# Storage
modes = ['coll_p1', 'coll_p2', 'uv']
ts = np.logspace(-12, 1, 100)
storage_f64 = [pd.DataFrame(columns=modes, index=ts)
               for _ in range(len(to_explore))]
storage_f128 = [pd.DataFrame(columns=modes, index=ts)
                for _ in range(len(to_explore))]
# Single point
uv_test = 1e10
l_test = [list(uv_test*np.array(loops_momenta[0]))] + loops_momenta[1:]
for t, d in zip(diags_titles, to_explore):
    print("{}\t: {}".format(t, color_out(
        abs(uv_test**2*complex(*d.worker.evaluate(l_test))))))
res = 0
for t, d in zip(diags_titles, to_explore):
    res += uv_test**2*complex(*d.worker.evaluate(l_test))
cprint(abs(res))
cprint(uv_test**2*abs(complex(*qqbar2a.worker.evaluate(l_test))))
# cprint(uv_test**2*complex(*to_explore[1].worker.evaluate(l_test)) +
#       uv_test**2*complex(*to_explore[2].worker.evaluate(l_test)))
# sys.exit(1)

# Plot grid
n_explore = len(to_explore)
n_modes = len(modes)
fig, axes = plt.subplots(nrows=n_explore, ncols=n_modes, figsize=(30, 15))
if n_explore == 1:
    axes = np.array([axes])
for axs, mode in zip(axes[0], modes):
    axs.set_title(mode, size='xx-large')
for axs, sgw in zip(axes[:, 0], to_explore):
    axs.set_ylabel(sgw.process_title, size='large')

# Plot
for plot_i, sgw in enumerate(to_explore):
    title = sgw.process_title
    rust_instance = sgw.worker
    n_loops = len(sgw.get_loop_momenta())

    for axs, mode in zip(axes[plot_i], modes):
        if mode == 'coll_p1':
            # Collinear explore
            ls = [[list(collinear_limit(p1, t, x=0.5))]
                  + loops_momenta[1:] for t in ts]
            rescaling_jacobians = np.array([t**1 for t in ts])
            lpsq = [sq(l[0]+p1) for l in ls]
            lpsq_label = 'lp1sq'
            if plot_i == n_explore - 1:
                axs.set_xlabel(r'$\lambda$', size='xx-large')

        elif mode == 'coll_p2':
            # Collinear explore
            ls = [[list(collinear_limit(p2, t, x=-0.5))]
                  + loops_momenta[1:] for t in ts]
            print(ls[0])
            rescaling_jacobians = np.array([t**1 for t in ts])
            lpsq = [sq(l[0]+p2) for l in ls]
            lpsq_label = 'lp2sq'
            if plot_i == n_explore - 1:
                axs.set_xlabel(r'$\lambda$', size='xx-large')

        elif mode == 'uv':
            # UV explore
            ls = [[list(1/t*np.array(loops_momenta[1]))] +
                  loops_momenta[1:] for t in ts]
            rescaling_jacobians = np.array([1/t**3 for t in ts])
            lpsq = [sp_spatial(l[0], l[0]) for l in ls]
            lpsq_label = 'llsq'
            if plot_i == n_explore - 1:
                axs.set_xlabel(r'$\frac{1}{\lambda}$', size='xx-large')
        else:
            raise KeyError(mode)

        # Evaluate
        e_f64 = [complex(*rust_instance.evaluate(l)) for l in ls]
        e_f128 = [complex(*rust_instance.evaluate_f128(l))
                  for l in ls]
        storage_f64[plot_i][mode] = e_f64
        storage_f128[plot_i][mode] = e_f128

        # Parametrization Jacobian
        par_jacobians = []
        for t in ts:
            jac = 1.0
            for n, l in enumerate(external_data['out_momenta']):
                jac *= rust_instance.inv_parameterize_f128(l, n, e_cm_sq)[3]
            par_jacobians += [1.0/jac]
        par_jacobians *= rescaling_jacobians

        # Plot the asymptotic behaviour
        try:
            popt, pcov = fit_exp_curve(ts,
                                       np.abs(e_f128)*par_jacobians,
                                       x_min=1e-7, x_max=1e-5)
            axs.plot(ts, np.exp(line(np.log(ts), *popt)),
                     '-m', linewidth=5, alpha=0.2, label=r"$\lambda^{%f}$" % popt[0])
            print("{} scaling\t: {}".format(mode, color_out(popt[0])))
        except ValueError as err:
            print("ValueError: ", err)
            print(" => Skipping fitting the curve")

        # Plot f128 and rescale by the Jacobian to match evaluate_integrand
        axs.plot(ts, [abs(v.real) for x, v in zip(
            ts, np.real(e_f128)*par_jacobians)], '-r', label="evaluate f128")
        axs.plot(ts, [abs(v.real) for x, v in zip(
            ts, np.imag(e_f128)*par_jacobians)], '--r')
        # Plot f64 and rescale by the Jacobian to match evaluate_integrand
        axs.plot(ts, [abs(v.real) for x, v in zip(
            ts, np.real(e_f64)*par_jacobians)], '-g', label="evaluate f64")
        axs.plot(ts, [abs(v.real) for x, v in zip(
            ts, np.imag(e_f64)*par_jacobians)], '--g')
        # Plot collinearity of the final states
        axs2 = axs.twinx()  # share same x-axis
        axs2.plot(ts, [abs(v) for v in lpsq], linestyle='dashed',
                  color='tab:blue', label=lpsq_label)
        axs2.set_yscale('log')
        axs2.tick_params(axis='y', labelcolor="tab:blue", labelsize='large')

        axs.set_yscale('log')
        axs.set_xscale('log')

        lines, labels = axs.get_legend_handles_labels()
        lines2, labels2 = axs2.get_legend_handles_labels()
        # axs2.legend(lines + lines2, labels + labels2, loc='lower center',
        #       bbox_to_anchor=(0.5, -.25), ncol=3, fancybox=True)
        axs2.legend(lines + lines2, labels + labels2, loc='lower center',
                    ncol=2, fancybox=True)
#
#    # ask matplotlib for the plotted objects and their labels
#    lines, labels = axs_bottom.get_legend_handles_labels()
#    lines2, labels2 = ax2.get_legend_handles_labels()
#    ax2.legend(lines + lines2, labels + labels2, loc='lower center',
#               bbox_to_anchor=(0.5, -.25), ncol=3, fancybox=True)
#    axs_top.set_title(title)
#    axs_top.set_yscale('log')
#    axs_top.set_xscale('log')
#    axs_top.set_xlabel(r'$\lambda$')
#    axs_top.legend()
#    axs_bottom.set_yscale('log')
#    axs_bottom.set_xscale('log')
#    axs_bottom.set_xlabel(r'$\lambda$')
#    axs_bottom.tick_params(axis='y', labelcolor="green")
# axs_bottom.legend()
#
# for n,s in enumerate(storage_f64):
# s.to_csv(to_explore_COLL[n][0].replace(">","_").replace("@","_").replace(" ","")+"_f64.csv")
# for n,s in enumerate(storage_f128):
# s.to_csv(to_explore_COLL[n][0].replace(">","_").replace("@","_").replace(" ","")+"_f128.csv")
# print("DONE")
# plt.title(limit)
# plt.legend(bbox_to_anchor=(0.75, 0.5))
fig.tight_layout()
plt.savefig("qqbar2a_regions.pdf")
# plt.show()
#
# ps = topology['external_kinematics']
# for p in ps:
# print(p)
