import ltd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


a = ltd.CrossSection("LTD/bu_squared.yaml", "LTD/hyperparameters.yaml")

data = [a.evaluate_f128([[1., 1., 1.], [-1.6, -5.6, -2.6 + (i + 1) * 1e-10], [2., 3., 4.]], [[1.,5.,2.], [1.,5.,2.]])[0] for i in range(10000)]
#data = [a.evaluate_f128([[2, 0.5, 2. + (i + 1) * 1e-6], [3.0, 0., 4.]], [[1,1,0], [1,1,0]])[1] for i in range(10000)]
#data = [a.evaluate_f128([[2, 0.5, 2.], [1.0, -0.5, 2. + (i + 1) * 1e-6]], [[1,1,0], [1,1,0]])[1] * ((i + 1) * 1e-8)**2 for i in range(10000)]

fig, ax = plt.subplots()
#ax.set_yscale('log')
ax.plot(data)
ax.grid()
plt.show()

