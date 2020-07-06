import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if True:
    pjoin = os.path.join
    root_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.insert(0, pjoin(root_path, os.path.pardir))
    from validate import get_result

pd.options.display.float_format = '{:+.5e}'.format
if __name__ == '__main__':
    results = pd.read_csv('analytic_results.csv')
    results['neval'] = ['N/A' for _ in range(len(results))]
    results['rust'] = ['N/A' for _ in range(len(results))]
    results['err'] = ['N/A' for _ in range(len(results))]
    results['err(%)'] = ['N/A' for _ in range(len(results))]
    results['rel'] = ['N/A' for _ in range(len(results))]
    results['\u03C3 away'] = ['N/A' for _ in range(len(results))]
    for (idx, topo) in results.iterrows():
        try:
            cuba_output = get_result(pjoin('results', "%s_res.dat"%topo['name']))
            if cuba_output[1] == 0.0:
                res = cuba_output[3]
                err = cuba_output[4]
            else:
                res = cuba_output[1]
                err = cuba_output[2]
            results.loc[idx, 'neval'] = cuba_output[0]
            results.loc[idx, 'rust'] = "{:+.5e}".format(res)
            results.loc[idx, 'err'] = "{:.2e}".format(err)
            results.loc[idx, 'err(%)'] = "{:.2%}".format(abs(err/res))
            results.loc[idx, '\u03C3 away'] = "{:.2f}".format(abs(res-topo['result'])/ err)
            results.loc[idx, 'rel'] = "{:.2}".format(abs((res-topo['result'])/topo['result']))
        except FileNotFoundError:
            print("no result found for %s"%topo['name'])
    
    try:
        print(results.to_markdown())
    except AttributeError:
        print(results)
