import random
import numpy as np
import itertools

random.seed(6)
num_max_rank = 1
n_denominators = 5
n_loops = 3
ll_ids = [[0, 0, 1], [0, 1, 0], [0, 1, -1],
          [1, 0, 0], [-1, 0, 1], [1, 1, 0], [1, -1, 2]]
#ll_ids = [[0, 1], [1, 0], [1, -1]]
#ll_ids = [[1]]

tensor_struct = [list(t) for t in itertools.combinations_with_replacement(
    range(num_max_rank+1), n_loops)]

print(tensor_struct)

# Define Inputs
denominator = [
    {
        'name': 'p%d' % i,
        'signature': ll_ids[random.randint(0, 2**n_loops-2)],
        # 'shift': [20*(random.random()-0.5) for j in range(4)],
        'shift': [random.randint(-100, 100) for j in range(2)],
        'pow': random.randint(1, 3),
    }
    for i in range(n_denominators)]
# numerator = [(c_id, random.random()) for c_id in tensor_struct]
numerator = [(c_id, random.randint(-100, 100)) for c_id in tensor_struct]

# Call this function to cast the (denominator, numerator) pair into a mathematica readable expression
# numerator = [([0, 1], 1)]


def print_expression(denominator, numerator):
    ss = '('
    for key, c in numerator:
        if key == [0]*n_loops:
            ss += '{:+f}'.format(c)
        else:
            ss += '{:+f}*{:s}'.format(c, ''.join(
                ['(k{}[0]^{})'.format(i, p) for i, p in enumerate(key) if p != 0]))
    ss += ')'
    for den in denominator:
        # print(den['signature'])
        k = [''.join(['{:+f}*k{}[{}]'.format(s, n, i)
                      for n, s in enumerate(den['signature'])]) for i in range(2)]
        x = '({}{:+f})^2'.format(k[0], den['shift'][0])
        x0 = ''.join(['-({}{:+f})^2'.format(k[i], den['shift'][i])
                      for i in range(1, 2)])
        ss += '/({}{})^{}'.format(x, x0, den['pow'])
    return ss


def higher_order_derivative(order: list, pden: list, pnum: list, ll_id=0, basis=None):
    # print("::", order, ll_id, pnum)
    res = ''

    # Define basis and transformation matrix
    if basis is None:
        ll_basis = [ll for n, ll in enumerate(ll_ids) if order[n] != 0]
        try:
            map_basis = np.matrix(ll_basis).I.transpose()
        except np.linalg.LinAlgError:
            raise Exception(
                "Basis for higher oder residue not invertible: {}".format(ll_basis))

        basis = (ll_basis, map_basis)
        #print("Higher order residue basis: {}".format(ll_basis))

    # Check if end of recursion
    if order[ll_id] <= 0:
        if len(order)-1 == ll_id:
            return print_expression(pden, [(pnum, 1)])
        else:
            return higher_order_derivative(order, pden, pnum, ll_id + 1, basis)

    order[ll_id] -= 1
    c_n = next(n for n, value in enumerate(
        np.dot(basis[1], ll_ids[ll_id]).flat) if value != 0)
    # Check if numerator is differentiable
    for l_n, k_to_c in enumerate(basis[1].transpose()):
        # print(k_to_c)
        v = k_to_c.take(c_n)[0, 0]
        # print(v)
        if pnum[l_n] != 0:
            pnum[l_n] -= 1
            res += '+%f' % (v * (pnum[l_n]+1)) + '*(' +\
                higher_order_derivative(
                order, pden, pnum, ll_id, basis)+')'
            # Reset
            pnum[l_n] += 1

    # Take derivative of denominator
    for den_n in range(len(pden)):
        # get diff factor
        d_coeff = - 2 * pden[den_n]['pow']\
            * np.dot(basis[1][c_n], pden[den_n]['signature'])[0, 0]

        pden[den_n]['pow'] += 1  # rise denominator

        # Sum over all new monomials
        res += '+ %f' % (d_coeff * pden[den_n]['shift'][0]) + '*(' + \
            higher_order_derivative(order, pden, pnum, ll_id, basis)+')'
        # print(pden[den_n]['signature'])
        for l_n, s in enumerate(pden[den_n]['signature']):
            pnum[l_n] += 1
            res += '+%f' % (d_coeff * s) + '*(' + \
                higher_order_derivative(
                order, pden, pnum, ll_id, basis)+')'
            # Reset
            pnum[l_n] -= 1
        pden[den_n]['pow'] -= 1
    order[ll_id] += 1
    return res


if __name__ == '__main__':
    # We cut over two propagators k^2 and l^2 rised to some power n1 and n2
    loop_momenta = [[random.random() for i in range(4)]
                    for n in range(n_loops)]
    for k in loop_momenta:
        k[0] = np.sqrt(sum(np.array(k[1:])**2))

    # Set the number of derivative for each loop_line
    # FIXME: no linear comination of loop momenta supported yet
    order_n = [0, 0, 1, 0, 1, 0, 1]
    #order_n = [1, 2, 0]
    #order_n = [3]
    print([den['signature'] for den in denominator])
    print('fun = {};'.format(print_expression(denominator, numerator)))

    result = ''
    for pnum, c in numerator:
        result += "+%d" % (c) + "(" + \
            higher_order_derivative(order_n, denominator, pnum)+")"
    print('test = {};'.format(result))
