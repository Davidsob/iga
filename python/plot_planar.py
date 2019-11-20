import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as color_map

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection    

def plot(mesh,verts,cpts,**kwargs):
  # get the arguments
  cmap = kwargs.get('cmap',color_map.get_cmap('rainbow'))
  size = kwargs.get('size', (10,10))
  has_field = kwargs.get('has_field',False)
  transparancy = kwargs.get('alpha', 0.3)

  # get bounds utility function
  def get_bounds(data,idx):
    mn,mx = 1e8,-1e8
    for d in data:
      arr = np.asarray(d)
      mn, mx = min(mn,arr[idx]), max(mx,arr[idx])
    return mn, mx

  fig,ax = plt.subplots()

# set up the mesh representation required
# by tricontourf
  polys = []
  verts = np.asarray(verts)
  mesh = np.asarray(mesh)
  tris = np.asarray([[0,1,2],[2,3,0]])
  for el in mesh:
    for t in tris:
      polys.append(el[t])

  # get the field if one is provided
  if has_field == False:
    field = np.ones(verts.shape[0])*0.5
  else:
    field = verts[:,-1]
    a,b = min(field),max(field)
    print 'Field range: [%g, %g]'%(a,b)

  # plot the surface
  ax.tricontourf(verts[:,0],verts[:,1],polys,field,10,cmap=cmap)

  # plot the facet edges 
  for el in mesh:
    x = verts[el,:]
    for i,j in ((0,1),(1,2),(2,3),(3,0)):
      ax.plot([x[i,0],x[j,0]],[x[i,1],x[j,1]], color='black')

  # plot vector data if available
  vec_data = kwargs.get("vector_data", None)
  if vec_data:
    print "Have vector data"
    pts,uvec,vvec = vec_data
    for p,u,v in zip(pts,uvec,vvec):
      print p
      x,y = p
      dux,duy = u
      dvx,dvy = v

      ax.plot([x,x+dux],[y,y+duy],color='red')
      ax.plot([x,x+dvx],[y,y+dvy],color='green')

## plot the controlpoints
  for p in cpts:
    plt.plot(p[0],p[1],color='red',marker='s',markersize=1)

## set the bounds and aspect ratio of the plot
  xmin,xmax = get_bounds(verts,0)
  ymin,ymax = get_bounds(verts,1)
  padx = (xmax-xmin)*0.1
  pady = (ymax-ymin)*0.1

  ax.set_xlim(xmin - padx, xmax + padx)
  ax.set_ylim(ymin - pady, ymax + pady)
  ax.set_aspect('equal')
  ax.set_xlabel('x-axis')
  ax.set_ylabel('y-axis')
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
    has_field = (len(x[0]) > 2)
    if len(sys.argv) == 4:
      p,u = open_vector_data(sys.argv[2])
      p,v = open_vector_data(sys.argv[3])
      plot(m,x,cpt,vector_data=(p,u,v),has_field=has_field)
    else:
      plot(m,x,cpt,has_field=has_field)