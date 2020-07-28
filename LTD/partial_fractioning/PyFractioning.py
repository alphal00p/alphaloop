import numpy as np
import itertools
import time
import os
import sys
import copy

class bcolors:
    HEADER = '\033[1;35;48m'
    OKBLUE = '\033[1;34;48m'
    OKGREEN = '\033[1;32;48m'
    WARNING = '\033[1;33;48m'
    FAIL = '\033[1;31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        self.BOLD = ''
        self.UNDERLINE = ''


###############################################################################
#
#    DECORATORS
#
def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        elapsed = (te-ts) * 1000

        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = elapsed
        else:
            if elapsed > 1000:
                print("\033[K{} PF elements generated in {}s"
                      .format(len(result), elapsed/1000))
            else:
                print("\033[K{} PF elements generated in {}ms"
                      .format(len(result), elapsed))
        return result
    return timed


def print_info(method):
    def pretty_print_topology(*args, **kw):
        n_props = args[0]
        signatures = args[1]
        n_loops = len(args[1][0])
        name = kw.get('name')

        # Print title
        blocks = 60
        if name == None:
            title = "{}-LOOP".format(n_loops)
        else:
            title = "{}-LOOP - {}".format(n_loops, name)
        side_blocks = (blocks-len(title))//2
        blocks -= len(title) % 2
        print("\n{0:}{2:}{1:}"
              .format(bcolors.BOLD, bcolors.ENDC, "="*blocks))
        print("{0:}={2:}{3:}{2:}{1:}="
              .format(bcolors.BOLD, bcolors.ENDC, " "*(side_blocks-1), title))
        print("{0:}{2:}{1:}\n{0:} Components{1:}:"
              .format(bcolors.BOLD, bcolors.ENDC, "="*blocks))

        idx = 0
        for (s, n) in zip(signatures, n_props):
            print("\t> signature {}".format(s))
            if idx == idx+n-1:
                print("\t  i: {}".format(idx))
            else:
                print("\t  i: {},..,{}".format(idx, idx+n-1))
            idx += n

            # Propagator with 4-momenta
            den = "\t  props: 1/(("
            for ki, l in enumerate(s):
                if l == 0:
                    continue
                elif l == 1:
                    den += "+{0:}k{2:}{1:}".format(bcolors.BOLD,
                                                   bcolors.ENDC, ki)
                elif l == -1:
                    den += "-{0:}k{2:}{1:}".format(bcolors.BOLD,
                                                   bcolors.ENDC, ki)
                else:
                    den += "{2:+}{0:}k{3:}{1:}".format(
                        bcolors.BOLD, bcolors.ENDC, l, ki)

            den += "+{0:}p{2:}{1:})^2-m{2:}^2)".format(
                bcolors.BOLD, bcolors.ENDC, "i")

            # Propagator with energies
            den += "\n\t\t = 1/(("
            for ki, l in enumerate(s):
                if l == 0:
                    continue
                elif l == 1:
                    den += "+{0:}k{2:}{1:}"\
                        .format(bcolors.OKBLUE, bcolors.ENDC, ki)
                elif l == -1:
                    den += "-{0:}k{2:}{1:}"\
                        .format(bcolors.OKBLUE, bcolors.ENDC, ki)
                else:
                    den += "{2::+}{0:}k{3:}{1:}"\
                        .format(bcolors.OKBLUE, bcolors.ENDC, l, ki)

            den += "+{0:}p{2:}{1:})^2-{0:}E{2:}{1:}^2)"\
                .format(bcolors.OKGREEN, bcolors.ENDC, "i")
            print(den)

        # Notation
        print(" {}Notation{}:".format(bcolors.BOLD, bcolors.ENDC))
        print(
            "\t{0:}ki{1:}, {0:}pi{1:}   : 4-vectors".format(bcolors.BOLD, bcolors.ENDC))
        print("\t{0:}ki{2:}, {1:}pi{2:}   : Energy components".format(
            bcolors.OKBLUE, bcolors.OKGREEN, bcolors.ENDC))
        print("\tvki, vpi : Spatial part of 4-vector")
        loop_energies_vector = '('
        for i in range(n_loops):
            if i+1 == n_loops:
                loop_energies_vector += "vk{})".format(i)
            else:
                loop_energies_vector += "vk{},".format(i)
        print("\t{}Ei{}^2     : ({}.signature)^2+vpi^2+mi^2\n".format(
            bcolors.OKGREEN, bcolors.ENDC, loop_energies_vector))

        # Prefactor
        print(" {}Prefactor{}:".format(bcolors.BOLD, bcolors.ENDC))
        energies_product = ""
        size = 0
        for i in range(sum(n_props)):
            size += len("2E{}".format(i))
            if i+1 == sum(n_props):
                energies_product += "2{}E{}{}".format(
                    bcolors.OKGREEN, i, bcolors.ENDC)
            else:
                energies_product += "2{}E{}{}*".format(
                    bcolors.OKGREEN, i, bcolors.ENDC)
                size += 1
        size += len("(2πi)^{}".format(n_loops))
        print("{4:}{3:}(2πi)^{0:}\n{4:}{2:}\n{4:}   {1:}\n".format(
            n_loops, energies_product, '-'*size, ' '*((size-5)//2), '\t'))

        return method(*args, **kw)
    return pretty_print_topology


def store_report(method):
    def store_it(*args, **kw):
        result = method(*args, **kw)
        n_loops = len(args[1][0])
        output_type = kw.get('output_type')
        name = kw.get('name')
        if name is None:
            name = ""

        if output_type == 'mathematica':
            file_name = 'pf_{}l_{}.m'.format(
                n_loops, name.lower().replace(' ', '_'))
            store_report_mathematica(result, file_name)
        elif output_type == 'yaml':
            file_name = 'pf_{}l_{}.yaml'.format(
                n_loops, name.lower().replace(' ', '_'))
            store_report_yaml(result, file_name)
        elif output_type == 'FORM':
            file_name = 'pf_{}l_{}.frm'.format(
                n_loops, name.lower().replace(' ', '_'))
            store_report_form(result, file_name)
        elif output_type == 'pickle':
            file_name = 'pf_{}l_{}.pkl'.format(
                n_loops, name.lower().replace(' ', '_'))
            store_report_pkl(result, file_name)
        else:
            print(
                "Skip storing. Valid output_type: ['mathematica','yaml','pickle','FORM']")
        return result

    def store_report_mathematica(res, file_name):
        ss = "{\n"
        t0 = time.time()
        size = len(res)

        for pf_id, (fact, prod, num) in enumerate(res):
            dt = int((time.time()-t0)*2) % 4
            if dt == 0:
                print("Storing as mathematica   ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 1:
                print("Storing as mathematica.  ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 2:
                print("Storing as mathematica.. ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 3:
                print("Storing as mathematica...({}/{})".format(pf_id+1, size), end="\r")

            # store factor
            ss += " {"
            ss += str(fact) + ", {"

            # store denominators
            if prod == []:
                ss += "1}, {"
            else:
                for n, r in enumerate(prod):
                    for (j, (e, p)) in enumerate(zip(r['energies'], r['shifts'])):
                        if e != 0 and p != 0:
                            ss += "{:+}*E{}{:+}*p{}".format(e, j, p, j)
                    if n+1 == len(prod):
                        ss += "}, {"
                    else:
                        ss += ", "

            # store numerator steps
            for n, zs in enumerate(num):
                ss += "{"
                for i, num_step in enumerate(zs):
                    for j, l in enumerate(num_step['lambdas']):
                        if l != 0:
                            ss += "{:+}*k{}".format(l, j)
                    for (j, (e, p)) in enumerate(zip(num_step['energies'], num_step['shifts'])):
                        if e != 0 and p != 0:
                            ss += "{:+}*E{}{:+}*p{}".format(e, j, p, j)

                    if i+1 == len(zs) and n+1 == len(num):
                        ss += "} }"
                    elif i+1 == len(zs):
                        ss += "}, "
                    else:
                        ss += ", "
            if pf_id+1 == len(res):
                ss += "}"
            else:
                ss += "},\n"
        ss += "\n}"
        with open("{}".format(file_name), 'w') as f:
            f.write(ss)
            f.close()
        print("\033[KStored in {}  ".format(file_name))

    def store_report_form(res, file_name):
        # For the moment convert coefficients to int
        def cast_int(f):
            if int(f)-f == 0:
                return int(f)
            else:
                raise ValueError("Cannont safe cast %f into int."%f)

        ss = "L F = \n"
        t0 = time.time()
        size = len(res)

        den_library = []
        for pf_id, (fact, prod, num) in enumerate(res):
            dt = int((time.time()-t0)*2) % 4
            if dt == 0:
                print("Storing as FORM   ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 1:
                print("Storing as FORM.  ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 2:
                print("Storing as FORM.. ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 3:
                print("Storing as FORM...({}/{})".format(pf_id+1, size), end="\r")

            # store factor
            ss += "%+d" % cast_int(fact)

            # store denominators
            if prod != []:
                den = ""
                for n, r in enumerate(prod):
                    den = ""
                    for (j, (e, p)) in enumerate(zip(r['energies'], r['shifts'])):
                        if e != 0 and p != 0:
                            den += "{:+}*E{}{:+}*p{}".format(cast_int(e), j, cast_int(p), j)
                    try:
                        d_idx = den_library.index(den)
                    except ValueError:
                        d_idx = len(den_library)
                        den_library += [den]
                    ss += "*invd%d" % d_idx

            # store numerator steps
            ss += "*num(ncmd("
            for n, zs in enumerate(num):
                for i, num_step in enumerate(zs):
                    for j, l in enumerate(num_step['lambdas']):
                        if l != 0:
                            ss += "{:+}*k{}".format(cast_int(l), j)
                    for (j, (e, p)) in enumerate(zip(num_step['energies'], num_step['shifts'])):
                        if e != 0 and p != 0:
                            ss += "{:+}*E{}{:+}*p{}".format(cast_int(e), j, cast_int(p), j)
                    if i+1 == len(zs) and n+1 == len(num):
                        ss += "))"
                    elif i+1 == len(zs):
                        ss += "), ncmd("
                    else:
                        ss += ", "
            ss += "\n"
        

        with open("{}".format(file_name.replace(".frm", ".den")), 'w') as f:
            ss_den = ""
            for n, den in enumerate(den_library):
                ss_den += "d%d = %s;\n" % (n, den)
            f.write(ss_den)
            f.close()
        with open("{}".format(file_name), 'w') as f:
            f.write("%s;" % ss)
            f.close()
        print("\033[KStored in {}  ".format(file_name))
 
    def store_report_yaml(res, file_name):
        try:
            import yaml
        except ModuleNotFoundError:
            print("Missing module 'yaml'. Skip save.")
            return
        t0 = time.time()
        size = len(res)

        with open(file_name, "w") as stream:

            for pf_id, x in enumerate(res):
                dt = int((time.time()-t0)*2) % 4
                if dt == 0:
                    print("Storing as yaml   ({}/{})".format(pf_id+1, size), end="\r")
                elif dt == 1:
                    print("Storing as yaml.  ({}/{})".format(pf_id+1, size), end="\r")
                elif dt == 2:
                    print("Storing as yaml.. ({}/{})".format(pf_id+1, size), end="\r")
                elif dt == 3:
                    print("Storing as yaml...({}/{})".format(pf_id+1, size), end="\r")
                res[pf_id] = {
                    'factor': float(x[0]),
                    'denominators': [{key: [float(v) for v in val] for key, val in y.items()} for y in x[1]],
                    'numerators': {"step{}".format(n): [{key: [float(v) for v in val] for key, val in yy.items()} for yy in y] for n, y in enumerate(x[2])}
                }
                # with open(file_name, "w") as stream:
                yaml.dump([res[pf_id]], stream, default_flow_style=None)
        print("\033[KStored in {}  ".format(file_name))

    def store_report_pkl(res, file_name):
        try:
            import pickle
        except ModuleNotFoundError:
            print("Missing module 'pickle'. Skip save.")
            return
        t0 = time.time()
        size = len(res)

        for pf_id, x in enumerate(res):
            dt = int((time.time()-t0)*2) % 4
            if dt == 0:
                print("Storing as pickle   ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 1:
                print("Storing as pickle.  ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 2:
                print("Storing as pickle.. ({}/{})".format(pf_id+1, size), end="\r")
            elif dt == 3:
                print("Storing as pickle...({}/{})".format(pf_id+1, size), end="\r")
            res[pf_id] = {
                'factor': float(x[0]),
                'denominators': [{key: [float(v) for v in val] for key, val in y.items()} for y in x[1]],
                'numerators': {"step{}".format(n): [{key: [float(v) for v in val] for key, val in yy.items()} for yy in y] for n, y in enumerate(x[2])}
            }
        # with open(file_name, "w") as stream:
        with open(file_name, "wb") as stream:
            pickle.dump(res, stream)
        print("\033[KStored in {}  ".format(file_name))
    return store_it


def print_pretty_report(method):
    def print_report(*args, **kw):
        result = method(*args, **kw)
        n_loops = len(args[1][0])

        if kw.get('verbose') is None or not kw.get('verbose'):
            return result

        for (pf_id, (fact, prod, num)) in enumerate(result):
            print("\t"+"-"*40)
            print("\t{}ID{}: {}".format(
                bcolors.BOLD + bcolors.WARNING, bcolors.ENDC, pf_id))
            print("\t{}Factor{}: {}".format(
                bcolors.BOLD + bcolors.WARNING, bcolors.ENDC, fact))

            # Denominators
            den_f = ""
            if len(prod) == 0:
                den_f = "1"
            for n in range(len(prod)):
                if n == 0:
                    den_f += "d{}".format(n)
                else:
                    den_f += "*d{}".format(n)

            print("\t{}Denominators{}: {}".format(
                bcolors.BOLD + bcolors.WARNING, bcolors.ENDC, den_f))
            
            for n, r in enumerate(prod):
                print("\t\td{} = ".format(n), end='')
                for (j, (e, p)) in enumerate(zip(r['energies'], r['shifts'])):
                    if e != 0 and p != 0:
                        print("{0:+}*{4:}E{1:}{5:}{2:+}*{4:}p{3:}{5:}"
                              .format(e, j, p, j, bcolors.OKGREEN, bcolors.ENDC),
                              end='')
                print("")

            # Numerator
            num_f = "\tN{}() from N0(".format(n_loops)
            for n in range(n_loops):
                if n+1 != n_loops:
                    num_f += "k{}, ".format(n)
                else:
                    num_f += "k{})".format(n)

            print("\t{}Numerator{}: {}".format(
                bcolors.BOLD + bcolors.WARNING, bcolors.ENDC, num_f))
            for n, zs in enumerate(num):
                f_to = "N{}(".format(n+1)
                for i in range(n+1, n_loops):
                    if i+1 != n_loops:
                        f_to += "k{}, ".format(i)
                    else:
                        f_to += "k{}".format(i)
                f_to += ")"
                f_from = "\tN{}({},".format(
                    n, ["z{}".format(i) for i in range(len(zs))])
                for i in range(n+1, n_loops):
                    if i+1 != n_loops:
                        f_from += "k{}, ".format(i)
                    else:
                        f_from += "k{}".format(i)
                f_from += ")"
                print('\t >{3:}Step {0:}: {1:} -> {2:}{4:}'
                      .format(n + 1, f_from, f_to, bcolors.HEADER, bcolors.ENDC))
                for i, num in enumerate(zs):
                    print("\t\tz{} = ".format(i), end='')
                    for j, l in enumerate(num['lambdas']):
                        if l != 0:
                            print("{:+}*{}k{}{}"
                                  .format(l, bcolors.OKBLUE, j, bcolors.ENDC),
                                  end='')
                    for (j, (e, p)) in enumerate(zip(num['energies'], num['shifts'])):
                        if e != 0 and p != 0:
                            print("{0:+}*{4:}E{1:}{5:}{2:+}*{4:}p{3:}{5:}"
                                  .format(e, j, p, j, bcolors.OKGREEN, bcolors.ENDC),
                                  end='')
                    print("")
        return result
    return print_report
#
#    END DECORATORS
#
###############################################################################


# Loop through all subset of lenth n of the set of
# all integers from 0 to set_size-1
def get_next_subset(subset, set_size, bottom_up=True):
    n = len(subset)
    # Move subset forward
    if bottom_up:
        for i in range(n):
            if i+1 == n or subset[i]+1 < subset[i+1]:
                subset[i] += 1
                for j in range(i):
                    subset[j] = j
                return subset[-1] < set_size

    # Move subset backward
    else:
        for i in range(n):
            if subset[i] > i:
                subset[i] -= 1
                for j in range(1, i+1):
                    subset[i-j] = subset[i] - j
                return subset[-1] < set_size
        return False


def get_next_set_pair(subset1, subset2, set_size):
    # Move subset1 forward and subset2 backward
    # if subset1 U subset2 == set
    # then this holds true also after the iteration
    return get_next_subset(subset1, set_size, True) and get_next_subset(subset2, set_size, False)


# Perform parital fractionig to remove the poles from the h_index
def partial_fractioning(h_index, e_index, num=[], num_rank=100):
    n_h = len(h_index)-1
    n_e = len(e_index)
    k = h_index[0]
    num += [k]
    if n_e == 0:
        if num_rank == 0 and len(h_index) != 1:
            return []
        else:
            return [[[], h_index.copy()]]
    if n_h == 0:
        return [[[(k, i) for i in e_index], num]]

    res = []
    splits = [n-1 for n in range(n_h)]
    sub_exp = [[] for n in range(n_e+n_h)]
    while get_next_subset(splits, n_e+n_h-1):
        sub_exp[0:splits[0]+1] = [(k, e_index[j]) for j in range(splits[0]+1)]
        for n in range(n_h):
            if n+1 != n_h:
                sub_exp[splits[n]+1: splits[n+1]+1] = [
                    (h_index[n+1], e_index[j-n]) for j in range(splits[n], splits[n+1])]
            else:
                sub_exp[splits[n]+1:] = [(h_index[n+1], e_index[j-n])
                                         for j in range(splits[n], n_e+n_h-1)]
        res += [[sub_exp.copy(), num.copy()]]

    if num_rank >= 1:
        res += partial_fractioning(h_index[1:], e_index, num, num_rank-1)

    return res


# Map the to new triplet based on which denominator was integrated over
def den_mapper(den_giver, den_receiver, residue_n):
    den_mapped = {}  # 'numerator': den_receiver['numerator']}
    factor = den_receiver['lambdas'][residue_n]/den_giver['lambdas'][residue_n]

    den_mapped['energies'] = den_receiver['energies'] - \
        factor*den_giver['energies']
    den_mapped['shifts'] = den_receiver['shifts'] - \
        factor*den_giver['shifts']
    den_mapped['lambdas'] = den_receiver['lambdas'] - \
        factor*den_giver['lambdas']
    return den_mapped


# Store the instruction for the evaluation of the numerator
def num_mapper(product, indices, residue_n):
    num_mapped = []

    for idx in indices:
        den_giver = product[idx]
        factor = -1/den_giver['lambdas'][residue_n]
        num_mapped += [
            {'lambdas': factor*den_giver['lambdas'],
             'shifts': factor*den_giver['shifts'],
             'energies': factor*den_giver['energies']}
        ]
        num_mapped[-1]['lambdas'][residue_n] = 0
    return num_mapped


# Compute all the terms that appear by applying partial fractioning to a
# topology that contains a set of loop lines defined by their signature
# and multiplicity
#
# Example (Double Box):
#  - loop_line 1:
#      - signature = [1,0] (aka k)
#      - n_props   = 3
#  - loop_line 2:
#      - signature = [1,-1] (aka k-l)
#      - n_props   = 1
#  - loop_line 3:
#      - signature = [0,1] (aka l)
#      - n_props   = 3
#
# Additional options:
#  - verbose: True/False
#       print all the contributions
#  - output_type: mathematica/yaml/pickle
#       choose the format in which to save the result
#  - name:
#       specify a name to be used when saving the file
#  - log_time: {}
#       store in log_time the elapsed time in ms
#
@store_report
@print_pretty_report
@timeit
@print_info
def integrate_energies(ll_n_props, signatures, **kwargs):
    # ll_n_props: number of propagator per loop line
    # signatures: signature of each loop line
    n_loops = len(signatures[0])

    # Check if the signatures are correclty defined
    if any(len(s) != n_loops for s in signatures):
        raise ValueError(
            "All signatures must have the same length. (number of loop variables)")

    # Check that there are no tree elements
    if any(all(v == 0 for v in s) for s in signatures):
        raise ValueError(
            "All signatures must contain at least one loop momentum dependence (tree element factorize)")

    # Check that ll_n_props has enough entires
    if len(ll_n_props) != len(signatures):
        raise ValueError(
            "Need to specify the number of propagators belonging to each signature")

    id_to_ll = []
    for n, n_props in enumerate(ll_n_props):
        id_to_ll += [n for _ in range(n_props)]
    n_props = sum(ll_n_props)
    dens_plus = [{'energies': np.array([0 if n != i else 1 for i in range(n_props)]),
                  'shifts': np.array([0 if n != i else 1 for i in range(n_props)]),
                  'lambdas': np.array(signatures[id_to_ll[n]].copy())} for n in range(n_props)]
    dens_minus = [{'energies': np.array([0 if n != i else -1 for i in range(n_props)]),
                   'shifts': np.array([0 if n != i else 1 for i in range(n_props)]),
                   'lambdas': np.array(signatures[id_to_ll[n]].copy())} for n in range(n_props)]

    # Take all combination of plus and minus energies from the first E_fractioning
    #    0: positive energy
    #    1: negative energy
    pf_res = []
    for choose in itertools.product([0, 1], repeat=n_props):
        product = [[], [[] for _ in range(n_loops)]]
        for dpm, which in zip(zip(dens_plus, dens_minus), choose):
            product[0] += [copy.deepcopy(dpm[which])]

        # factor coming from the initial partial fractioning into positive and negative energies
        global_factor = (-1)**choose.count(0)
        pf_res += pf_product(product, n_loops, global_factor=global_factor)

        print("Select: {} (Total: {})".format(choose, len(pf_res)), end='\r')
    return pf_res


# Perform partial fractioning on a single product for an arbitrary number of loops
def pf_product(product, n_loops, r=0, global_factor=1.0, numerator=[]):
    if r == n_loops:
        # Return the result with the corresponding factor in the case that all
        # the momenta have been removed by the procedure
        if not all([all(den['lambdas'] == 0) for den in product[0]]):
            raise("There are still some loop momenta left")
        return [(global_factor, [{k: v for k, v in x.items() if k != 'lambdas'} for x in product[0]], product[1])]

    indices = []
    h_index = []
    e_index = []
    left_product = []
    for n, den in enumerate(product[0]):
        factor = den['lambdas'][r]
        if factor != 0:  # Element depneds on the loop momenta
            # Extract the factor in front of the loop momenta in order to
            # to have a consistent unit factor before taking the residue
            global_factor /= factor
            den['energies'] = den['energies'] / factor
            den['shifts'] = den['shifts'] / factor
            den['lambdas'] = den['lambdas'] / factor

            indices += [n]
            if all(den['energies'] >= 0):
                e_index += [n]
            else:
                h_index += [n]
        else:
            left_product += [den]

    # factor coming from applying partial fractioning to remove hyperboloids
    global_factor *= (-1)**len(h_index)
    if len(h_index) == 0:  # no pole
        return []
    else:  # apply residue and partial fractioni)g
        res = partial_fractioning(h_index, e_index, num=[], num_rank=1000)
    if res == []:
        return res

    result = []
    for mapping in res:
        new_product = [copy.deepcopy(left_product), copy.deepcopy(product[1])]
        for pair in mapping[0]:
            new_product[0] += [den_mapper(product[0][pair[0]],
                                          product[0][pair[1]],
                                          r)]
        new_product[1][r] = num_mapper(product[0], mapping[1], r)

        num = numerator + [mapping[1]]
        result += pf_product(new_product, n_loops, r+1,
                             global_factor=global_factor, numerator=num)
    return result


if __name__ == '__main__':
    bcolors = bcolors()
    # Uncomment next line to disable colors
    # bcolors.disable()

    #########################################
    #            1-LOOP - Box               #
    #########################################
    #              _  _____  _
    #                |     |
    #                |     |
    #              _ '_____' _
    #
    signatures = [[1]]
    n_props = [4]

    res = integrate_energies(n_props, signatures,
                             verbose=False,
                             name='Box',
                             output_type='mathematica')
    #########################################
    #            2-LOOP - Sunrise           #
    #########################################
    #                  ___
    #             _  /_____\ _
    #                \_____/
    #
    signatures = [[1, 0], [1, -1], [0, 1]]
    n_props = [1, 1, 1]

    res = integrate_energies(n_props, signatures,
                             verbose=False,
                             name='Sunrise',
                             output_type='mathematica')

    #########################################
    #            2-LOOP - PentaBox          #
    #########################################
    #
    #                 |___  _
    #          _  ___/    |
    #            |   |    |
    #          _ |___|____| _
    #
    signatures = [[1, 0], [1, -1], [0, 1]]
    n_props = [1, 3, 4]

    res = integrate_energies(n_props, signatures,
                             verbose=False,
                             name="PentaBox",
                             output_type='mathematica')

    #########################################
    #            4-LOOP - Banana            #
    #########################################
    #              _____
    #             / ___ \
    #          _ (/_____\)_
    #            (\_____/)
    #             \_____/
    #
    signatures = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [-1, -1, -1, 1]]
    n_props = [1, 1, 1, 1, 1]

    res = integrate_energies(n_props, signatures,
                             verbose=False,
                             name="banana",
                             output_type='mathematica')

    #########################################
    #            6-LOOP - Diamond           #
    #########################################
    #              _________
    #           _ /_|_____|_\ _
    #             '. \   / .'
    #               '.\ /.'
    #                 '.'
    #

    signatures = [[1, 0, 0, 0, 0, 0],
                  [1, -1, 0, 0, 0, 0],
                  [0, 1, 0, 0, 0, 0],
                  [0, 1, -1, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0],
                  [0, 0, 1, -1, 0, 0],
                  [0, 0, 0, 1, 0, 0],
                  [0, 0, 0, 1, -1, 0],
                  [0, 0, 0, 0, 1, 0],
                  [-1, 0, 0, 0, 0, 1],
                  [0, 0, -1, 0, 0, 1],
                  [0, 0, 0, 0, -1, 1]]
    n_props = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    res = integrate_energies(n_props, signatures,
                             verbose=False,
                             name="diamond",
                             output_type='mathematica')

    #########################################
    #         4-LOOP -  2x2 FISHNET         #
    #########################################
    #           _ ._________. _
    #             |    |    |
    #             |____|____|
    #             |    |    |
    #           _ '____|____' _
    #
    signatures = [[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [1, -1, 0, 0],
                  [-1, 0, -1, 0],
                  [0, -1, 0, -1],
                  [0, 0, 1, 0],
                  [0, 0, -1, 1],
                  [0, 0, 0, -1]]
    n_props = [2, 2, 1, 1, 1, 2, 1, 2]

    res = integrate_energies(n_props, signatures,
                             name="2x2 Fishnet",
                             output_type='pickle')
