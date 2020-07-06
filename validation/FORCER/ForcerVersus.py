import sys
import os
from pathlib import Path
import re
from scipy.special import zeta
import numpy as np
import pandas as pd

if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, pjoin(root_path, os.path.pardir))
    sys.path.insert(0, pjoin(root_path, os.path.pardir, os.path.pardir))
    from LTD.squared_topologies import SquaredTopologyGenerator

ST_FORCER_LIBRARY = {}

ST_FORCER_2P_3L = []
ST_FORCER_LIBRARY[3] = ST_FORCER_2P_3L
ST_FORCER_2P_3L.append(([(0, 3), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4), (4, 5)], "+ 2*ep^-1*z3"))
ST_FORCER_2P_3L.append(([(0, 1), (1, 2), (1, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (5, 6)], "- 2*ep^-1*z3"))
ST_FORCER_2P_3L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 5), (4, 6), (5, 6), (2, 7)], "- 2*ep^-1*z3"))

ST_FORCER_2P_4L = []
ST_FORCER_LIBRARY[4] = ST_FORCER_2P_4L
ST_FORCER_2P_4L.append(([(0, 3), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5), (3, 5), (4, 5), (4, 6)], "+5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 3), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5), (4, 6)], "+5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 6), (3, 4), (3, 6), (4, 5), (5, 6), (6, 7)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 5), (3, 4), (3, 6), (4, 6), (5, 6), (6, 7)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 4), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 7)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 5), (2, 3), (2, 6), (3, 4), (3, 6), (4, 5), (4, 6), (5, 6), (5, 7)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 6), (6, 7), (7, 8)], "+5/4*ep^-1*z3"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (7, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 7), (4, 6), (5, 6), (5, 7), (6, 7), (7, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 6), (4, 5), (4, 6), (5, 7), (6, 7), (2, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 4), (3, 5), (4, 6), (5, 6), (5, 7), (6, 7), (2, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (2, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (2, 8)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 3), (2, 7), (3, 4), (3, 7), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7), (2, 8)], "+5/4*ep^-1*z3"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 8), (6, 7), (6, 8), (2, 9)], "-5*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 8), (6, 7), (7, 8), (2, 9)], "-10*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 7), (5, 8), (6, 8), (7, 8), (2, 9)], "-10*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-5/2*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 7), (5, 8), (6, 8), (7, 8), (2, 9)], "-5/2*ep^-1*z5"))
ST_FORCER_2P_4L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "+5/12*ep^-1*z3"))

ST_FORCER_2P_5L = []
ST_FORCER_LIBRARY[5] = ST_FORCER_2P_5L
ST_FORCER_2P_5L.append(([(0, 3), (1, 3), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (4, 7)], "+36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 3), (1, 5), (1, 6), (2, 4), (2, 5), (2, 6), (3, 4), (3, 5), (4, 6), (5, 6), (4, 7)], "+441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 3), (1, 5), (1, 6), (2, 3), (2, 5), (2, 6), (3, 4), (4, 5), (4, 6), (5, 6), (4, 7)], "+36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 3), (1, 4), (1, 5), (2, 3), (2, 5), (2, 6), (3, 6), (4, 5), (4, 6), (5, 6), (4, 7)], "+441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 2), (1, 5), (1, 6), (2, 5), (2, 6), (3, 4), (3, 5), (3, 6), (4, 5), (4, 6), (4, 7)], "+36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 2), (1, 3), (1, 5), (2, 5), (2, 6), (3, 4), (3, 6), (4, 5), (4, 6), (5, 6), (4, 7)], "+441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 2), (1, 3), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6), (4, 7)], "+441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 3), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 5), (4, 6), (5, 6), (4, 7)], "+36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 5), (1, 7), (2, 6), (2, 7), (3, 4), (3, 5), (3, 7), (4, 6), (4, 7), (5, 6), (6, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 5), (1, 7), (2, 5), (2, 7), (3, 4), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (6, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 7), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 7), (3, 6), (3, 7), (4, 5), (4, 7), (5, 6), (6, 7), (6, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 4), (2, 7), (3, 5), (3, 7), (4, 6), (4, 7), (5, 6), (6, 7), (6, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 7), (6, 7), (6, 8)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 5), (1, 2), (1, 3), (1, 4), (2, 5), (2, 7), (3, 5), (3, 7), (4, 6), (4, 7), (5, 6), (6, 7), (6, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 6), (1, 7), (2, 3), (2, 4), (2, 5), (3, 4), (3, 6), (4, 7), (5, 6), (5, 7), (6, 7), (5, 8)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 6), (2, 7), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 7), (5, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 5), (2, 7), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 4), (2, 7), (3, 5), (3, 6), (4, 6), (4, 7), (5, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 4), (2, 6), (3, 5), (3, 7), (4, 6), (4, 7), (5, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 5), (1, 6), (2, 3), (2, 4), (2, 5), (3, 6), (3, 7), (4, 6), (4, 7), (5, 7), (6, 7), (5, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 6), (2, 7), (3, 5), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (6, 7), (5, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 6), (2, 7), (3, 4), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (5, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 6), (2, 7), (3, 4), (3, 5), (3, 6), (4, 6), (4, 7), (5, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 3), (2, 6), (3, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 3), (2, 6), (3, 4), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (5, 8)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 5), (2, 3), (2, 4), (3, 6), (3, 7), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (5, 8)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (5, 8)], "+12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7), (2, 8)], "+12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 7), (2, 3), (2, 7), (3, 8), (4, 5), (4, 6), (4, 7), (5, 6), (5, 8), (6, 8), (8, 9)], "-24*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 7), (3, 6), (3, 8), (4, 7), (4, 8), (5, 6), (5, 8), (6, 7), (8, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 7), (3, 6), (3, 8), (4, 7), (4, 8), (5, 6), (5, 7), (6, 8), (8, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 7), (3, 6), (3, 7), (4, 7), (4, 8), (5, 6), (5, 8), (6, 8), (8, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 7), (3, 5), (3, 8), (4, 6), (4, 7), (5, 6), (6, 8), (7, 8), (8, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 8), (6, 7), (7, 8), (8, 9)], "+6*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 7), (5, 8), (6, 8), (7, 8), (8, 9)], "+6*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 6), (4, 7), (4, 8), (5, 7), (6, 8), (7, 8), (8, 9)], "+3/2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 7), (1, 2), (1, 3), (1, 4), (2, 3), (2, 5), (3, 7), (4, 5), (4, 6), (5, 8), (6, 7), (6, 8), (7, 8), (8, 9)], "+3/2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 7), (1, 8), (2, 3), (2, 4), (2, 5), (3, 6), (3, 8), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 7), (1, 8), (2, 3), (2, 4), (2, 5), (3, 6), (3, 7), (4, 6), (4, 8), (5, 7), (5, 8), (6, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 7), (1, 8), (2, 3), (2, 4), (2, 5), (3, 4), (3, 6), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (7, 9)], "-3/4*ep^-1*z3"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 8), (3, 7), (3, 8), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 8), (7, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 7), (3, 4), (3, 8), (4, 5), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 7), (3, 4), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 7), (3, 4), (3, 8), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 7), (3, 4), (3, 5), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 4), (3, 5), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (6, 7), (7, 8), (7, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 4), (3, 5), (3, 7), (4, 6), (4, 8), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 4), (3, 5), (3, 7), (4, 6), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (7, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 4), (3, 5), (3, 7), (4, 6), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8), (7, 9)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 8), (2, 3), (2, 4), (3, 4), (3, 7), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 7), (2, 8), (3, 4), (3, 5), (3, 6), (4, 5), (4, 8), (5, 8), (6, 7), (6, 8), (7, 9)], "-3/4*ep^-1*z3"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 8), (3, 4), (3, 7), (4, 5), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 8), (3, 4), (3, 5), (4, 6), (4, 8), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 8), (3, 4), (3, 5), (4, 6), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 7), (3, 4), (3, 8), (4, 5), (4, 6), (5, 6), (5, 8), (6, 8), (7, 8), (7, 9)], "-3/4*ep^-1*z3"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (3, 8), (4, 7), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 3), (2, 4), (3, 5), (3, 7), (4, 6), (4, 8), (5, 6), (5, 8), (6, 8), (7, 8), (7, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 8), (3, 7), (3, 8), (4, 5), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 8), (3, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 8), (3, 7), (3, 8), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 5), (3, 8), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 5), (3, 8), (4, 6), (4, 8), (5, 6), (5, 8), (6, 7), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 5), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 5), (3, 8), (4, 5), (4, 8), (5, 6), (6, 7), (6, 8), (7, 8), (7, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 5), (3, 8), (4, 5), (4, 6), (5, 8), (6, 7), (6, 8), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 4), (3, 8), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 7), (3, 4), (3, 8), (4, 5), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 7), (3, 8), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 8), (4, 5), (4, 7), (5, 8), (6, 7), (6, 8), (7, 8), (7, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7), (4, 6), (4, 8), (5, 7), (5, 8), (6, 8), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 7), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (7, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 3), (2, 8), (3, 8), (4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 8), (7, 8), (7, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (7, 9)], "+12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 7), (2, 8), (3, 4), (3, 8), (4, 5), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 7), (2, 8), (3, 4), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 7), (6, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 7), (2, 8), (3, 4), (3, 8), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 7), (2, 8), (3, 4), (3, 5), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 8), (3, 5), (3, 6), (4, 7), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (2, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 8), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 8), (4, 6), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 8), (4, 5), (4, 8), (5, 6), (6, 7), (6, 8), (7, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 5), (3, 6), (4, 5), (4, 8), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 4), (3, 8), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 4), (2, 7), (3, 4), (3, 5), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 3), (2, 8), (3, 4), (4, 7), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 3), (2, 8), (3, 4), (4, 5), (4, 7), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 7), (2, 3), (2, 8), (3, 4), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 7), (2, 8), (3, 5), (3, 6), (4, 7), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8), (2, 9)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 8), (4, 7), (4, 8), (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 8), (4, 5), (4, 8), (5, 6), (5, 7), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 8), (4, 5), (4, 7), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 8), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 5), (4, 7), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 7), (3, 5), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 7), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 8), (3, 4), (3, 5), (3, 7), (4, 6), (4, 8), (5, 6), (5, 8), (6, 7), (7, 8), (2, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 7), (2, 8), (3, 4), (3, 5), (3, 7), (4, 6), (4, 7), (5, 6), (5, 8), (6, 8), (7, 8), (2, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 7), (3, 7), (3, 8), (4, 5), (4, 6), (4, 8), (5, 7), (5, 8), (6, 7), (6, 8), (2, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 7), (3, 4), (3, 7), (4, 5), (4, 8), (5, 6), (5, 8), (6, 7), (6, 8), (7, 8), (2, 9)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 9), (3, 4), (3, 5), (4, 6), (4, 7), (5, 6), (5, 8), (6, 7), (7, 8), (8, 9), (9, 10)], "-1/12*ep^-1*z3"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 7), (4, 9), (5, 7), (5, 8), (6, 7), (6, 8), (8, 9), (9, 10)], "+1/2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 7), (4, 9), (5, 6), (5, 7), (6, 8), (7, 8), (8, 9), (9, 10)], "+1/2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 9), (8, 9), (9, 10)], "+2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 9), (5, 7), (6, 7), (6, 8), (7, 8), (8, 9), (9, 10)], "+1/2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 8), (6, 9), (7, 8), (7, 9), (9, 10)], "+2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 9), (2, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 7), (6, 8), (7, 9), (8, 9), (9, 10)], "+ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 4), (3, 9), (4, 6), (5, 6), (5, 7), (6, 8), (7, 8), (7, 9), (8, 9), (9, 10)], "-7*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (4, 5), (5, 6), (5, 9), (6, 7), (6, 8), (7, 8), (7, 9), (8, 9), (9, 10)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 9), (2, 4), (2, 9), (3, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 8), (6, 7), (6, 9), (8, 9), (2, 10)], "-7*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 9), (2, 4), (2, 9), (3, 5), (3, 6), (4, 5), (4, 7), (5, 8), (6, 8), (6, 9), (7, 8), (7, 9), (2, 10)], "+2*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 9), (3, 6), (3, 7), (4, 6), (4, 9), (5, 6), (5, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-14*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 9), (3, 6), (3, 7), (4, 6), (4, 8), (5, 7), (5, 9), (6, 8), (7, 9), (8, 9), (2, 10)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 9), (3, 5), (3, 6), (4, 6), (4, 9), (5, 7), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-1001/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 9), (3, 5), (3, 6), (4, 6), (4, 7), (5, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 9), (3, 4), (3, 6), (4, 6), (5, 7), (5, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 8), (4, 9), (5, 8), (6, 7), (6, 9), (7, 9), (8, 9), (2, 10)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 9), (5, 9), (6, 8), (6, 9), (7, 8), (8, 9), (2, 10)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 9), (5, 8), (6, 8), (6, 9), (7, 9), (8, 9), (2, 10)], "-36/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 7), (4, 8), (5, 9), (6, 8), (6, 9), (7, 9), (8, 9), (2, 10)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 5), (3, 7), (4, 6), (4, 8), (5, 9), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/20*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 5), (2, 6), (3, 4), (3, 9), (4, 9), (5, 7), (5, 8), (6, 7), (6, 9), (7, 8), (8, 9), (2, 10)], "-12*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 9), (3, 5), (4, 6), (4, 7), (5, 8), (5, 9), (6, 8), (6, 9), (7, 8), (7, 9), (2, 10)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 9), (3, 5), (4, 6), (4, 7), (5, 6), (5, 9), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/80*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 9), (3, 5), (4, 6), (4, 7), (5, 6), (5, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-7*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 9), (3, 5), (4, 5), (4, 6), (5, 7), (6, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-7*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 9), (4, 6), (4, 7), (5, 8), (5, 9), (6, 8), (6, 9), (7, 8), (7, 9), (2, 10)], "+3*ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 9), (4, 6), (4, 7), (5, 6), (5, 9), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/80*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 9), (4, 6), (4, 7), (5, 6), (5, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-161/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 9), (5, 8), (5, 9), (6, 7), (6, 9), (7, 8), (8, 9), (2, 10)], "-18/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 9), (5, 8), (5, 9), (6, 7), (6, 8), (7, 9), (8, 9), (2, 10)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 9), (5, 7), (5, 9), (6, 8), (6, 9), (7, 8), (8, 9), (2, 10)], "-441/40*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 7), (4, 8), (5, 7), (5, 9), (6, 8), (6, 9), (7, 9), (8, 9), (2, 10)], "-18/5*ep^-1*z3^2"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 9), (5, 7), (5, 9), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/80*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 3), (1, 4), (2, 3), (2, 5), (3, 6), (4, 6), (4, 7), (5, 8), (5, 9), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "-441/80*ep^-1*z7"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 9), (3, 4), (3, 9), (4, 5), (4, 6), (5, 7), (5, 8), (6, 7), (6, 9), (7, 8), (8, 9), (2, 10)], "-1/4*ep^-1*z3"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 9), (3, 4), (3, 5), (4, 6), (4, 7), (5, 8), (5, 9), (6, 8), (6, 9), (7, 8), (7, 9), (2, 10)], "+ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 9), (3, 4), (3, 5), (4, 6), (4, 7), (5, 6), (5, 9), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "+ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 9), (5, 7), (6, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "+ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (3, 5), (3, 6), (4, 5), (4, 7), (5, 9), (6, 8), (6, 9), (7, 8), (7, 9), (8, 9), (2, 10)], "+ep^-1*z5"))
ST_FORCER_2P_5L.append(([(0, 1), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5), (4, 9), (5, 6), (5, 9), (6, 7), (6, 8), (7, 8), (7, 9), (8, 9), (2, 10)], "-1/4*ep^-1*z3"))

def analytic_result(pole_string, n_loops):
    r = eval(re.sub(r'z(\d+)',r'zeta(\1)', pole_string.replace("*ep^-1","").replace("^","**")))
    return r * n_loops * np.pi/(16.0 * np.pi**2)**n_loops

if __name__ == '__main__':
    Path(pjoin(root_path,'Rust_inputs')).mkdir(parents=True, exist_ok=True)
    analytic_results = {'name':[], 'result':[]}
    for n_loops in range(3,5+1):
        for n, (edges_raw, pole_string) in enumerate(ST_FORCER_LIBRARY[n_loops]):
            name = "STF_%dL_%d"%(n_loops, n)
            print("Creating %s"%name, end='\r')
            edges = [('p%d'%i , *e) for i, e in enumerate(edges_raw)]
            edges[0] = ('q1',*edges[0][1:])
            edges[-1] = ('q2',*edges[-1][1:])
            incoming_momenta_names = ['q1']
            n_jets=0
            analytic_results['name'] += [name]
            analytic_results['result'] += [analytic_result(pole_string, n_loops)]
            external_momenta = {'q1': [1.0, 0.0, 0.0 ,0.0], 'q2':[1.0, 0.0, 0.0, 0.0]}
            stf_2P = SquaredTopologyGenerator(edges, 
                                              name, 
                                              incoming_momenta_names,
                                              n_jets,
                                              external_momenta,
                                              overall_numerator=1.0,
    )
            stf_2P.export(pjoin(root_path,'Rust_inputs','%s.yaml'%name))
    print("\x1b[2K\rInput YAML files for cross_section can be found in Rust_inputs")
    
    pd.DataFrame(analytic_results).to_csv("analytic_results.csv")
    
