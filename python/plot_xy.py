import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as color_map

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection    

def plot(xy, **kwargs):
  xy = np.asarray(xy)
  plt.plot(xy[:,0],xy[:,1])
  plt.grid(True)
  print 'Field range: [%g, %g]'%(min(xy[:,1]),max(xy[:,1]))
  plt.show()

def open_xy_data(file):
  xy = []
  with open(file) as f:
    n = map(int,f.readline().split())[0]
    while n > 0 :
      xy.append(map(float, f.readline().split()))
      n -= 1

  return xy

if __name__ == "__main__":
  if len(sys.argv) > 1:
    xy = open_xy_data(sys.argv[1])
    plot(xy)
