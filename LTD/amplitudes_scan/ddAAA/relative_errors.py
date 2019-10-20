#!/usr/bin/python

import numpy as np
import pandas as pd

if __name__ == "__main__":
    #input_filename = "ddAAA_cuhre_2M_results.csv" 
    #output_filename ="ddAAA_cuhre_2M.csv"
    
    input_filename = "ddAAA_vegas_100M_results.csv" 
    output_filename ="ddAAA_vegas_100M.csv"

    dr = pd.read_csv (input_filename)
    da = pd.read_csv (r'analytic_results.csv')
    #Write seed id
    relative = [np.array(dr['seed'].to_numpy())]
    
    #Write relative real error
    relative += [[ x/y for (x,y) in zip(dr["real"].to_numpy(),da['ep0_re'].to_numpy())]]
    relative += [[ abs(x/y) for (x,y) in zip(dr["real_err"].to_numpy(),dr['real'].to_numpy())]]
    #Write relative imag error
    relative += [[ x/y for (x,y) in zip(dr["imag"].to_numpy(),da['ep0_im'].to_numpy())]]
    relative += [[ abs(x/y) for (x,y) in zip(dr["imag_err"].to_numpy(),dr['imag'].to_numpy())]]

    relative = np.array(relative).transpose()
    data = pd.DataFrame(relative, index=None, columns=['seed','re_rel','re_err_rel','im_rel','im_err_rel'])
    data.to_csv(output_filename, encoding='utf-8', index=False)

