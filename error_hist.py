#!/usr/bin/env python
import sys
import numpy as np
import matplotlib as mat
mat.use('agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x = np.fromfile(sys.argv[1])

mu, sigma = -0.942, 7.34

# the histogram of the data
n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel(r'$ \frac{\hat{x}-x}{x10^{-6}}$')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ Error:}\ \mu=-0.942,\ \sigma=7.34$')
plt.axis([-50, 50, 0, 0.1])
plt.grid(True)

plt.savefig(sys.argv[2],bbox_inches='tight')
