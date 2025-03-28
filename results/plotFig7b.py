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


filename1 = 'hd5dv3.mat'
figname = 'figure7b.eps'
Figdata1 = scipy.io.loadmat(filename1)
blockList1 = Figdata1["n_range"]
n_points = len(blockList1[0])

qber = Figdata1["Q"]

keyrate1 = np.array(Figdata1["maxKeyList"][0])

theory1 = np.log2(5) + 2 * (1-qber) * np.log2(1-qber) + 2*qber*np.log2(qber/4)
theory1 = theory1 * np.ones((n_points, 1))

filename2 = 'hd7dv3.mat'

Figdata2 = scipy.io.loadmat(filename2)
blockList2 = Figdata2["n_range"]
n_points = len(blockList2[0])

theory2 = np.log2(7) + 2 * (1-qber) * np.log2(1-qber) + 2*qber*np.log2(qber/6)
theory2 = theory2 * np.ones((n_points, 1))

keyrate2 = np.array(Figdata2["maxKeyList"][0])


filename3 = 'ps.mat'

Figdata3 = scipy.io.loadmat(filename3)
blockList3 = Figdata3["n_range"]


keyrate3 = np.array(Figdata3["d5key"][0])


blockList4 = Figdata3["n_range"]


keyrate4 = np.array(Figdata3["d7key"][0])

# Adjust the size of the figure below.
# If the value is too small, the legend may have overlaps with curves.
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
curve1 = ax.semilogx(blockList1[0], np.maximum(
    keyrate1, 0), 'ro', linewidth=2, markersize=8)
curve2 = ax.semilogx(blockList3[0], np.maximum(
    keyrate3, 0), 'gs', linewidth=2, markersize=8)
curve3 = ax.semilogx(blockList1[0], theory1, 'y', linewidth=2, markersize=8)
curve4 = ax.semilogx(blockList2[0], np.maximum(
    keyrate2, 0), 'b*', linewidth=2, markersize=8)
curve5 = ax.semilogx(blockList4[0], np.maximum(
    keyrate4, 0), 'cx', linewidth=2, markersize=8)
curve6 = ax.semilogx(blockList2[0], theory2,
                     'purple', linewidth=2, markersize=8)
plt.xlabel(r'Number of signals')
plt.ylabel('Key rate per signal')


# Set the legend of curves.

curves = curve1+curve2+curve3+curve4+curve5+curve6
ax.legend(curves, [r'd = 5, EAT', r'd = 5, postselection', r'd = 5, asymptotic',
          r'd = 7, EAT', r'd = 7, postselection', r'd = 7, asymptotic'], frameon=False)


# Remove the grid in the plot.
plt.grid(False)


plt.tight_layout()

# Save the figure using the eps file format
# Modify this part as you need.

plt.savefig(figname, format='eps', dpi=1200)
