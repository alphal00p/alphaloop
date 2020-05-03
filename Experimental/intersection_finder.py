import math
import numpy as np
from scipy import linalg
import random

def square(d):
    return d.dot(d)

def f_df(k, foci, q0):
    """ Evaluate the surface and its gradient at point `k`"""
    f = -q0
    df = [np.array([0., 0., 0.]) for _ in k]
    for (signature, shift, mass) in foci:
        d = np.array(signature).dot(k)
        norm = np.sqrt(square(d + shift) + mass**2)
        f += norm

        mompart = (d + shift) / norm
        for li in range(len(k)):
            if signature[li] != 0:
                for i in range(3):
                    df[li][i] += mompart[i]
    return f, np.array(df)

def find_t(k, dir, foci, q0):
    """Use Newton's method. The projection can only have two solutions, so we are guaranteed that there
    are no local minima."""
    t = 200 # find the positive t solution
    for i in range(20):
        # evaluate the surface and the derivative at k = dir * t + k
        f = -q0
        df = 0.
        for (signature, shift, mass) in foci:
            d = np.array(signature).dot(dir)
            k_in_shift = np.array(signature).dot(k) + shift
            s = np.sqrt(square(d * t + k_in_shift) + mass**2)
            f += s
            df += (t * square(d) + d.dot(k_in_shift)) / s

        new_t = t - f / df

        #print(new_t, f, df)

        if abs(f) < 1e-14: # TODO: scale with e_cm
            #print('Found %f in %d steps:' % (t, i))
            return t
        t = new_t

    print("No convergence")
    return None

def get_random_point(foci, q0):
    # construct a point on the inside of the e-surface
    rows = [f[0] for f in foci[:-1]]
    nullspace = linalg.null_space(rows)
    if nullspace.shape[1] > 0:
        rows = np.vstack([rows, list(nullspace.T)])
    mat = np.linalg.inv(rows)

    shifts = [-f[1] for f in foci[:-1]] + [np.array([0., 0., 0.,])] * (len(foci[0][0]) - len(foci) + 1)
    k_in = mat.dot(np.array(shifts))

    dir = [np.array([1., 1., 1.]) for _ in shifts]

    t = find_t(k_in, dir, foci, q0)
    return k_in + np.array(dir) * t

def intersection_finder(surfaces):
    # start with a random point on the surface
    x = get_random_point(surfaces[0][:-1], surfaces[0][-1])
    y = get_random_point(surfaces[1][:-1], surfaces[1][-1])

    for _ in range(20):
        print('current distance', np.linalg.norm(x-y))

        # find the halfway point
        f_x, df_x = f_df(x, surfaces[0][:-1], surfaces[0][-1])
        f_y, df_y = f_df(y, surfaces[1][:-1], surfaces[1][-1])
        assert(abs(f_x) < 1e-13 and abs(f_y) < 1e-13)
        #print('funcs', f_x, f_y, df_x, df_y)

        # now find a point inside the volume
        t = find_t(x, -df_x, surfaces[0][:-1], surfaces[0][-1])
        assert(t > 0.)
        center_x = x - df_x * t * 0.1

        t = find_t(y, -df_y, surfaces[1][:-1], surfaces[1][-1])
        assert(t > 0.)
        center_y = y - df_y * t * 0.1

        #print('centers', center_x, center_y)

        # project the center line to 2 new points
        t_x = find_t(x, center_y - center_x, surfaces[0][:-1], surfaces[0][-1])
        t_y = find_t(y, center_x - center_y, surfaces[1][:-1], surfaces[1][-1])

        #print('tx', t_x, t_y)

        if t_x >= 0.5 and t_y >= 0.5:
            # we are crossing! 
            print('Intersecting!')
            break

        # set the new points
        x = x + (center_y - center_x) * t_x
        y = y + (center_x - center_y) * t_y



# example problem
q0 = 20.
q1 = 10.
p1 = np.array([0.0, 0., 0.])
p2 = np.array([1.0, 2., 0.])
p3 = np.array([0.0, 1., 0.])
p4 = np.array([2.0, 3., 0.])
p5 = np.array([4.0, 0., 2.])
m1 = 1.0
m2 = 1.0
m3 = 1.0
m4 = 1.0
m5 = 1.0

surfaces = [[([1, 0], p1, m1), ([0, 1], p2, m2), ([1, 1], p3, m3), q0], 
            [([1, 0], p4, m4), ([1, 0], p5, m5), q1]]

#surfaces = [[([1], p1, m1), ([1], p2, m2), q0], 
#            [([1], p4, m4), ([1], p5, m5), q0]]


intersection_finder(surfaces)