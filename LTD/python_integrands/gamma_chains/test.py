#!/usr/bin/env python

import gamma_chain
from numpy import array

MAX_LENGTH = int(gamma_chain.gamma_chain_max_length.maxlength)

vbar    =   [3.0,3.1,3.2,3.3]
u       =   [4.0,4.1,4.2,4.3]
indices =   [1,-1,-2,2,-1,-2,-3,-4,-5,-6,-5,-4,-6,-3]
vectors =   [
    [1.0,1.1,1.2,1.3],
    [2.0,2.1,2.2,2.3],
]

print("Testing benchmark chain of 14 gamma matrices:")
result = gamma_chain.compute_chain( vbar, u, indices, vectors)
target_result = complex(-78354.8416+1312.256j)
print(result)
print("Target benchmark result:")
print(target_result)
