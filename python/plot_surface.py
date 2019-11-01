import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as color_map

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection    

def plot(mesh,verts,cpts,**kwargs):
  cmap = kwargs.get('cmap',color_map.get_cmap('rainbow'))
  size = kwargs.get('size', (10,10))
  transparancy = kwargs.get('alpha', 1.0)
  # get bounds
  def get_bounds(data,idx):
    mn,mx = 1e8,-1e8
    for d in data:
      arr = np.asarray(d)
      mn, mx = min(mn,arr[idx]), max(mx,arr[idx])
    return mn, mx

  fig = plt.figure(figsize=size)
  ax = Axes3D(fig)

  polys = []
  poly_verts = []
  k = 0
  color = [0,0,0,transparancy]
  color[:3] = cmap(0.5)[:3]
  for el in mesh:
    vert = []
    for i in el:
      vert.append(verts[i])
    poly_verts.append(vert)
  polys.append(Poly3DCollection(poly_verts,edgecolors=('black',),facecolors=(color,)))

  xmin,xmax = get_bounds(verts,0)
  ymin,ymax = get_bounds(verts,1)
  zmin,zmax = get_bounds(verts,2)
  padx = (xmax-xmin)*0.01
  pady = (ymax-ymin)*0.01
  padz = (zmax-zmin)*0.01
    # Create cubic bounding box to simulate equal aspect ratio
  max_range = np.array([xmax-xmin + 2*padx, ymax-ymin + 2*pady, zmax-zmin + 2*padz]).max()
  Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(xmax+xmin)
  Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(ymax+ymin)
  Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(zmax+zmin)

  # Comment or uncomment following both lines to test the fake bounding box:
  for xb, yb, zb in zip(Xb, Yb, Zb):
    ax.plot([xb], [yb], [zb], 'w')

  # coll = Poly3DCollection(verts, edgecolors=('black',))
  for poly in polys:
    ax.add_collection3d(poly,zs='z')
  # ax.add_collection3d(coll,zs='z')

  for p in cpts:
    ax.scatter(p[0],p[1],p[2],color='red',marker='s')
  ax.set_xlabel('x-axis')
  ax.set_ylabel('y-axis')
  ax.set_zlabel('z-axis')
  plt.show()
  return fig,ax


def open_file(file):
  mesh = []
  x = []
  cpts = []
  with open(file) as f:
    n = map(int,f.readline().split())[0]
    while n > 0 :
      mesh.append(map(int, f.readline().split()))
      n -= 1

    n = map(int,f.readline().split())[0]
    while n > 0 :
      x.append(map(float, f.readline().split()))
      n -= 1

    n = map(int,f.readline().split())[0]
    while n > 0 :
      cpts.append(map(float, f.readline().split()))
      n -= 1
  return mesh, x, cpts

if __name__ == "__main__":
  if len(sys.argv) > 1:
    m,x,cpt = open_file(sys.argv[1])
    plot(m,x,cpt)