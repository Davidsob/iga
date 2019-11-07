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
  has_field = kwargs.get('has_field',False)
  transparancy = kwargs.get('alpha', 0.3)
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
  field = []
  for el in mesh:
    vert = []
    tmp = []
    for i in el:
      vert.append(verts[i][:3])
      if has_field:
        tmp.append(verts[i][3])
    poly_verts.append(vert)
    field.append(np.mean(tmp))

  if has_field == False:
    color = [0,0,0,transparancy]
    color[:3] = cmap(0.5)[:3]
    polys.append(Poly3DCollection(poly_verts,edgecolors=('black',),facecolors=(color,)))
  else:
    # scale field
    a,b = min(field),max(field)
    colors = map(lambda c: cmap((c-a)/(b-a)), field)
    polys.append(Poly3DCollection(poly_verts,edgecolors=('black',),facecolors=colors))


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

  for poly in polys:
    ax.add_collection3d(poly)
    # ax.add_collection3d(poly,zs='z')

  vec_data = kwargs.get("vector_data", None)
  if vec_data:
    print "Have vector data"
    pts,uvec,vvec = vec_data
    for p,u,v in zip(pts,uvec,vvec):
      w = np.cross(u,v)
      w *= np.linalg.norm(u)/np.linalg.norm(w) 
      x,y,z = p
      dux,duy,duz = u
      dvx,dvy,dvz = v
      dwx,dwy,dwz = w

      ax.plot3D([x,x+dux],[y,y+duy],[z,z+duz],color='red')
      ax.plot3D([x,x+dvx],[y,y+dvy],[z,z+dvz],color='green')
      ax.plot3D([x,x+dwx],[y,y+dwy],[z,z+dwz],color='blue')

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
    has_field = (len(x[0]) > 3)
    if len(sys.argv) == 4:
      p,u = open_vector_data(sys.argv[2])
      p,v = open_vector_data(sys.argv[3])
      plot(m,x,cpt,vector_data=(p,u,v),has_field=has_field)
    else:
      plot(m,x,cpt,has_field=has_field)