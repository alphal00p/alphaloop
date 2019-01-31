#!/usr/bin/env python2
import ltd

l1 = ltd.LTD("P7")
print l1.evaluate([0.1, 0.2, 0.3])
print l1.evaluate([0.2, 0.3, 0.4])
print l1.deform(loop_momentum=[2.0, 3.0, 4.0])

l2 = ltd.LTD("double-triangle")
print l2.evaluate_cut(x=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7], cut_structure=[1,1,0], cut_indices=[0,0,-1])