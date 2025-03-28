# This code works for python 3.
# This codes is a Python notebook file and may be used with Jupyter notebook.
# If you use python GUI, can run from there
# Make sure you have installed numpy, scipy and matplotlib
# Installation guide of matplotlib can be found from https://matplotlib.org/


import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import math

# see https://matplotlib.org/ for documentation
# use Google to search the function you need.


filename1 = 'bb84v3alg2exploret.mat'
figname = 'BB84tplot.eps'
Figdata1 = scipy.io.loadmat(filename1)
blockList1 = Figdata1["n_range"][0]

keyrate1 = Figdata1["maxKeyList"]

# Adjust the size of the figure below. If you want to make it a square,
# use the same number for two inputs;
# If the value is too small, your legend may have overlaps with your curves.

plt.figure(figsize=(15, 9))
plt.rcParams.update({'font.size': 25})

ax = plt.subplot(111)


# The first symbol is for the color of the curve.
# r : red, k : black, g : green, c : cyan, y : yellow, b : blue
# See https://matplotlib.org/3.1.1/api/colors_api.html for more colors.

# The second symbol is for the style of the marker. Some of the styles are similar to MATLAB.
# See https://matplotlib.org/3.1.1/api/markers_api.html for more descriptions.
#
# The third symbol is for the style of the line.
# See https://matplotlib.org/gallery/lines_bars_and_markers/line_styles_reference.html for details.
#
plt.ticklabel_format(style='sci', axis='y')

curve1 = ax.semilogx(blockList1, keyrate1[0], 'ro', linewidth=2, markersize=8)
curve2 = ax.semilogx(blockList1, keyrate1[1], color=(1, 0.5, 0), marker='*', linestyle='none', linewidth=2,
                     markersize=8)
curve3 = ax.semilogx(blockList1, keyrate1[2], 'c^', linewidth=2, markersize=8)
curve4 = ax.semilogx(blockList1, keyrate1[3], 'bs', linewidth=2, markersize=7)
curve5 = ax.semilogx(blockList1, keyrate1[4], 'g>', linewidth=2, markersize=8)

plt.xlabel(r'Number of signals')
plt.ylabel('Key rate per signal')


# Set the legend of curves.

curves = curve1 + curve2 + curve3 + curve4 + curve5
ax.legend(curves, ['$\zeta_{t}$=0', '$\zeta_{t}$=0.005', '$\zeta_{t}$=0.01', '$\zeta_{t}$=0.015', '$\zeta_{t}$=0.02'],
          frameon=False, loc='best')

# Remove the grid in the plot. If you need the grid, can set it to be True.
plt.grid(False)


plt.tight_layout()

# save the figure using the eps file format
# Modify this part as you need.

plt.savefig(figname, format='eps', dpi=1200)

# plt.show()
