import numpy as np
import itertools
import time
import os
import copy
import sys
import progressbar


# For the moment convert coefficients to int
def cast_int(f):
    if int(f)-f == 0:
        return int(f)
    else:
        raise ValueError("Cannont safe cast %f into int." % f)


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
                      .format(len(result.pf_res), elapsed/1000))
            else:
                print("\033[K{} PF elements generated in {}ms"
                      .format(len(result.pf_res), elapsed))
        return result
    return timed


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
# @timeit
class PartialFractioning:
    def __init__(self, ll_n_props, signatures, **kwargs):
        # ll_n_props: number of propagator per loop line
        # signatures: signature of each loop line
        self.signatures = signatures
        self.ll_n_props = ll_n_props
        self.n_loops = len(signatures[0])
        self.name = kwargs.pop("name", "topo").lower().replace(' ', '_')
        self.shift_map = kwargs.pop("shift_map", None)
        self.n_cuts = kwargs.pop("n_cuts", 0)
        self.den_library = []
        self.esq_library = [[], []]
        self.energies_map = []
        self.file_name = 'pf_{}l_{}.frm'.format(self.n_loops, self.name)

        # Check if the signatures are correclty defined
        if any(len(s) != self.n_loops for s in signatures):
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

        self.id_to_ll = []
        for n, n_props in enumerate(ll_n_props):
            self.id_to_ll += [(n, i) for i in range(n_props)]
        n_props = sum(ll_n_props)
        dens_plus = [{'energies': np.array([0 if n != i else 1 for i in range(n_props)]),
                      'shifts': np.array([0 if n != i else 1 for i in range(n_props)]),
                      'lambdas': np.array(signatures[self.id_to_ll[n][0]].copy())} for n in range(n_props)]
        dens_minus = [{'energies': np.array([0 if n != i else -1 for i in range(n_props)]),
                       'shifts': np.array([0 if n != i else 1 for i in range(n_props)]),
                       'lambdas': np.array(signatures[self.id_to_ll[n][0]].copy())} for n in range(n_props)]
        # Take all combination of plus and minus energies from the first E_fractioning
        #    0: positive energy
        #    1: negative energy
        den_library = []
        self.pf_res = []

        with progressbar.ProgressBar(prefix='PF starting energies: {variables.choose} : ', max_value=2**n_props, variables={'choose': 'N/A'}) as bar:
            for n_choose, choose in enumerate(itertools.product([0, 1], repeat=n_props)):
                bar.update(n_choose)
                bar.update(choose="".join([str(c) for c in choose]))
                product = [[], [[] for _ in range(self.n_loops)]]
                for dpm, which in zip(zip(dens_plus, dens_minus), choose):
                    product[0] += [dpm[which]]

                # factor coming from the initial partial fractioning into positive and negative energies
                global_factor = (-1)**choose.count(0)
                self.pf_res += self.pf_product(product,
                                               global_factor=global_factor)

    # Map the to new triplet based on which denominator was integrated over
    def den_mapper(self, den_giver, den_receiver, residue_n):
        den_mapped = {}  # 'numerator': den_receiver['numerator']}
        factor = den_receiver['lambdas'][residue_n] / \
            den_giver['lambdas'][residue_n]

        den_mapped['energies'] = den_receiver['energies'] - \
            factor*den_giver['energies']
        den_mapped['shifts'] = den_receiver['shifts'] - \
            factor*den_giver['shifts']
        den_mapped['lambdas'] = den_receiver['lambdas'] - \
            factor*den_giver['lambdas']
        return den_mapped

    # Store the instruction for the evaluation of the numerator
    def num_mapper(self, product, indices, residue_n):
        num_mapped = []

        for idx in indices:
            den_giver = product[idx]
            factor = 1/den_giver['lambdas'][residue_n]
            num_mapped += [
                {'lambdas': - factor*den_giver['lambdas'],
                 'shifts': factor*den_giver['shifts'],
                 'energies': factor*den_giver['energies']}
            ]
            num_mapped[-1]['lambdas'][residue_n] = 0
        return num_mapped

    # Perform partial fractioning on a single product for an arbitrary number of loops

    def pf_product(self, product, r=0, global_factor=1.0):
        if r == self.n_loops:
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
            res = self.partial_fractioning(
                h_index, e_index, num=[], num_rank=1000)
        if res == []:
            return res

        result = []
        for mapping in res:
            new_product = [copy.deepcopy(
                left_product), copy.deepcopy(product[1])]
            for pair in mapping[0]:
                new_product[0] += [self.den_mapper(product[0][pair[0]],
                                                   product[0][pair[1]],
                                                   r)]
            new_product[1][r] = self.num_mapper(product[0], mapping[1], r)

            result += self.pf_product(new_product, r+1,
                                      global_factor=global_factor)
        return result

    # Perform parital fractionig to remove the poles from the h_index
    def partial_fractioning(self, h_index, e_index, num=[], num_rank=100):
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
            sub_exp[0:splits[0]+1] = [(k, e_index[j])
                                      for j in range(splits[0]+1)]
            for n in range(n_h):
                if n+1 != n_h:
                    sub_exp[splits[n]+1: splits[n+1]+1] = [
                        (h_index[n+1], e_index[j-n]) for j in range(splits[n], splits[n+1])]
                else:
                    sub_exp[splits[n]+1:] = [(h_index[n+1], e_index[j-n])
                                             for j in range(splits[n], n_e+n_h-1)]
            res += [[sub_exp.copy(), num.copy()]]

        if num_rank >= 1:
            res += self.partial_fractioning(h_index[1:],
                                            e_index, num, num_rank-1)

        return res

    def to_FORM(self):
        ss = ""
        for pf_id, (fact, prod, num) in enumerate(self.pf_res):
            # store factor
            ss += "%+d" % cast_int(fact)
            # store denominators
            if prod != []:
                den = ""
                for n, r in enumerate(prod):
                    den = ""
                    for (j, e) in enumerate(r['energies']):
                        if e != 0:
                            den += "{:+}*E{}".format(cast_int(e), j)
                    shifts = ""
                    for (j, p) in enumerate(r['shifts']):
                        if p != 0:
                            if j < self.n_cuts:
                                shifts += "{:+}*c{}".format(cast_int(p), j+1)
                            else:
                                shifts += "{:+}*p{}".format(cast_int(p),
                                                            j-self.n_cuts+1)
                    if shifts != "":
                        den += "+energies(%s)" % shifts
                    try:
                        d_idx = self.den_library.index(den)
                    except ValueError:
                        d_idx = len(self.den_library)
                        self.den_library.append(den)
                    ss += "*invd%d" % d_idx
            # store numerator steps
            ss += "*num(ncmd("
            for n, zs in enumerate(num):
                for i, num_step in enumerate(zs):
                    for j, l in enumerate(num_step['lambdas']):
                        if l != 0:
                            ss += "{:+}*k{}".format(cast_int(l), j)
                    for (j, e) in enumerate(num_step['energies']):
                        if e != 0:
                            ss += "{:+}*E{}".format(cast_int(e), j)
                    shifts = ""
                    for (j, p) in enumerate(num_step['shifts']):
                        if p != 0:
                            if j < self.n_cuts:
                                shifts += "{:+}*c{}".format(cast_int(p), j+1)
                            else:
                                shifts += "{:+}*p{}".format(cast_int(p),
                                                            j-self.n_cuts+1)
                    if shifts != "":
                        ss += "+energies(%s)" % shifts
                    if i+1 == len(zs) and n+1 == len(num):
                        ss += ") )"
                    elif i+1 == len(zs):
                        ss += "), ncmd("
                    else:
                        ss += ", "
            ss += "\n"
        return ss

    def shifts_to_externals(self):
        if self.shift_map is None:
            return
        # Remap all the shifts
        if not isinstance(self.shift_map, np.ndarray):
            self.shift_map = np.array(self.shift_map)
        for f, dens, num in self.pf_res:
            for d in dens:
                d['shifts'] = self.shift_map.dot(d['shifts'])
            for ncmd in num:
                for lim in ncmd:
                    lim['shifts'] = self.shift_map.dot(lim['shifts'])
        
        # Check for degenerate energies
        for idx, (sign, shift) in enumerate(self.id_to_ll):
            esq = {"signature": self.signatures[sign], 'shifts': self.shift_map.transpose()[
                idx]}
            try:
                e_idx = next(i for i, v in enumerate(self.esq_library[0])
                             if all(s1 == s2 for s1, s2 in zip(v['signature'], esq['signature']))
                             and all(s1 == s2 for s1, s2 in zip(v['shifts'], esq['shifts'])))
                self.energies_map += [[0 if i!= self.esq_library[1][e_idx] else 1 for i in range(sum(self.ll_n_props))]]
            except StopIteration:
                self.energies_map += [[0 if i!= idx else 1 for i in range(sum(self.ll_n_props))]]
                self.esq_library[0].append(esq)
                self.esq_library[1].append(idx)
        self.energies_map = np.array(self.energies_map).transpose()
        # Remap energies
        if len(self.esq_library[1]) == sum(self.ll_n_props):
            return
        print(self.energies_map)
        for f, dens, num in self.pf_res:
            for d in dens:
                d['energies'] = self.energies_map.dot(d['energies'])
            for ncmd in num:
                for lim in ncmd:
                    lim['energies'] = self.energies_map.dot(lim['energies'])
 

if __name__ == '__main__':
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
    shift_map = [[0, 0, 0, -1],
                 [0, 1, 1, 0],
                 [0, 0, 1, 0]]

    # TEST vacuum 
    shift_map = [[0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0]]
    pf = PartialFractioning(n_props, signatures,
                            name='Box', shift_map=shift_map, n_cuts=0)
    pf.shifts_to_externals()
    print("1L Result:")
    print(pf.to_FORM())
    for n, den in enumerate(pf.den_library):
        print("d%d = %s;" % (n, den))
    #sys.exit()
    #    #########################################
    #    #            2-LOOP - Sunrise           #
    #    #########################################
    #    #                  ___
    #    #             _  /_____\ _
    #    #                \_____/
    #    #
    #    signatures = [[1, 0], [1, -1], [0, 1]]
    #    n_props = [1, 1, 1]
    #
    #    res = integrate_energies(n_props, signatures,
    #                             verbose=False,
    #                             name='Sunrise',
    #                             output_type='mathematica')
    #
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
    shift_map = [[0, -1, 0, 0, 0, 0, 0, -1],
                 [0, 0, 0, 1, 1, 1, 1, 0],
                 [0, 0, 0, 0, 0, 1, 1, 0],
                 [0, 0, 0, 0, 0, 0, 1, 0]]

    pf = PartialFractioning(n_props, signatures,
                            name='PentaBox', shift_map=shift_map, n_cuts=1)
    pf.shifts_to_externals()
    output = "L F =\n"
    output += pf.to_FORM()
    output += ";\n"
    for idx, s in enumerate(pf.signatures):
        print("ll {}: {}".format(idx, s))
    print()
    for idx, id_to_ll in enumerate(pf.id_to_ll):
        print("id {} :> ll {}".format(idx, id_to_ll))
    print()
    for n, den in enumerate(pf.den_library):
        print("d%d = %s;" % (n, den))
