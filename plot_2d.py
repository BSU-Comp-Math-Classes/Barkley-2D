import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate
from mpl_toolkits.mplot3d.axes3d import Axes3D


filename = str(sys.argv[1])
with open(filename, 'r') as f:
    next(f)
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
    u = [float(line.split()[2]) for line in lines]


nx = int(math.sqrt(len(x)))

# Set up a regular grid of interpolation points
xi, yi = np.linspace(min(x), max(x), nx), np.linspace(min(y), max(y), nx)
xi, yi = np.meshgrid(xi, yi)

# Interpolate; there's also method='cubic' for 2-D data such as here
ui = scipy.interpolate.griddata((x, y), u, (xi, yi), method='linear')


fig = plt.figure()

plt.contourf(xi,yi,ui,extent=[min(x), max(x), min(y), max(y)])

plt.show()
