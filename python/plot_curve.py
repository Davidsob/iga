import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as color_map

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection    

def plot3(mesh, verts, cpts, **kwargs):
  size = kwargs.get('size', (10,10))
  fig = plt.figure(figsize=size)
  ax = Axes3D(fig)

  for el in mesh:
    a1,a2,a3 = verts[el[0]]
    b1,b2,b3 = verts[el[1]]
    ax.plot3D([a1,b1],[a2,b2],[a3,b3], 'k')

  for p in cpts:
    x,y,z = p
    ax.scatter(x,y,z,color='red',marker='s')

  vec_data = kwargs.get("vector_data", None)
  if vec_data:
    print "Have vector data"
    pts,vecs = vec_data
    for p,v in zip(pts,vecs):
      x,y,z = p
      dx,dy,dz = v
      ax.plot3D([x,x+dx],[y,y+dy],[z,z+dz],color='red')
      ax.scatter(x,y,color='red')

  ax.set_xlabel('x-axis')
  ax.set_ylabel('y-axis')
  ax.set_zlabel('z-axis')

  plt.show()

def plot(mesh, verts, cpts, **kwargs):
  for el in mesh:
    a = verts[el[0]]
    b = verts[el[1]]
    plt.plot([a[0], b[0]],[a[1], b[1]], 'k')

  for p in cpts:
    plt.scatter(p[0],p[1],color='red',marker='s')

  vec_data = kwargs.get("vector_data", None)
  if vec_data:
    print "Have vector data"
    pts,vecs = vec_data
    for p,v in zip(pts,vecs):
      x,y = p
      dx,dy = v
      plt.arrow(x,y,dx,dy,color='green',head_width=0.05, head_length=0.05)
      plt.scatter(x,y,color='green')

  plt.grid(True)
  plt.axis('equal')
  plt.show()


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

def open_vector_data(file):
  x = []
  vecs = []
  with open(file) as f:
    n = map(int,f.readline().split())[0]
    while n > 0 :
      x.append(map(float, f.readline().split()))
      n -= 1
      
    n = map(int,f.readline().split())[0]
    while n > 0 :
      vecs.append(map(float, f.readline().split()))
      n -= 1
  return x, vecs

if __name__ == "__main__":
  if len(sys.argv) > 1:
    m,x,cpt = open_file(sys.argv[1])
    dim = len(x[0])
    plotter = plot if dim == 2 else plot3
    if len(sys.argv) > 2:
      vec_data = open_vector_data(sys.argv[2])
      plotter(m,x,cpt,vector_data=vec_data)
    else:
      plotter(m,x,cpt)