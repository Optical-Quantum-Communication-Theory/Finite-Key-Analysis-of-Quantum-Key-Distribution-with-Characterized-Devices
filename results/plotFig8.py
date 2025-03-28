# This code works for python 3.
# In command line, type python3 plotFig.py
# If you use python GUI, can run from there
# Make sure you have installed numpy, scipy and matplotlib
# Installation guide of matplotlib can be found from https://matplotlib.org/


import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import math


# see https://matplotlib.org/ for documentation
# use Google to search the function you need.


filename1 = 'opBB84v3L0t.mat'
figname = 'opticalBB84-revised.eps'
Figdata1 = scipy.io.loadmat(filename1)
blockList1 = Figdata1["n_range"]


keyrate1 = Figdata1["maxKeyList"]

filename2 = 'opBB84v3L10t.mat'

Figdata2 = scipy.io.loadmat(filename2)


keyrate2 = Figdata2["maxKeyList"]


filename3 = 'opBB84v3L20t.mat'
Figdata3 = scipy.io.loadmat(filename3)


keyrate3 = Figdata3["maxKeyList"]


# Adjust the size of the figure below.
# If the value is too small, the legend may overlap with curves.
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
curve1 = ax.semilogx(blockList1[0], keyrate1[0],
                     'ro', linewidth=2, markersize=8)
curve2 = ax.semilogx(blockList1[0], keyrate2[0],
                     'g*', linewidth=2, markersize=8)
curve3 = ax.semilogx(blockList1[0], keyrate3[0],
                     'cd', linewidth=2, markersize=8)
plt.xlabel(r'Number of signals')
plt.ylabel('Key rate per signal')


# Set the legend of curves.

curves = curve1+curve2+curve3
ax.legend(curves, [r'L = 0 km', r'L = 10 km', r'L = 20 km'],
          frameon=False, loc='best')

# Remove the grid in the plot. If you need the grid, can set it to be True.
plt.grid(False)


plt.tight_layout()

# Save the figure using the eps file format
# Modify this part as you need.

plt.savefig(figname, format='eps', dpi=1200)
