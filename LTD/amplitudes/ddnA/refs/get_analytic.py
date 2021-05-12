#!/bin/python
import sys

import numpy as np
import pandas as pd

if len(sys.argv) != 2:
    print("Must give an amplitude id")
    sys.exit(1)

amp_id = int(sys.argv[1])

df_ct = pd.read_csv('./ddnA_ct0.csv', index_col=0, dtype=np.float64)
ct = df_ct[df_ct.index == amp_id]

df_res = pd.read_csv('./analytic_results.csv', index_col=0, dtype=np.float64)
res = df_res[df_res.index == amp_id]


analytic_results = complex(res['ep0_re'], res['ep0_im'])
analytic_ct = complex(ct['ct0_re'], ct['ct0_im'])
print("analytic resutl : {:+.16e}".format(analytic_results))
print("analytic CT     : {:+.16e}".format(analytic_ct))
print("\nanalytic + ct   : {:+.16e}".format(analytic_results + analytic_ct))
