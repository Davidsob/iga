#!/u/xeons08/people/bdavidson/anaconda2/bin
import csv
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as color_map

from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection    
from matplotlib.animation import FuncAnimation as animate

def plot(mesh,verts,**kwargs):
  cmap = kwargs.get('cmap',color_map.get_cmap('rainbow'))
  transparancy = kwargs.get('alpha', 0.3)
  field_range = kwargs.get('field_range',(0,1));
  size = kwargs.get('size', (10,10))
  fig = kwargs.get('figure',None)
  ax  = kwargs.get('axis', None)
  # get bounds
  def get_bounds(data,idx):
    mn,mx = 1e8,-1e8
    for d in data:
      arr = np.asarray(d)
      mn, mx = min(mn,arr[idx]), max(mx,arr[idx])
    return mn, mx

  if not fig:
    fig = plt.figure(figsize=size)
  if not ax:
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
      tmp.append(verts[i][3])
    poly_verts.append(vert)
    field.append(np.mean(tmp))

  # scale field
  a,b = min(field),max(field)
  fmn,fmx = field_range
  colors = map(lambda c: cmap((c-fmn)/(fmx-fmn)), field)
  polys.append(Poly3DCollection(poly_verts,edgecolors=('black',),facecolors=colors))
  print 'Field range: [%g, %g]'%(a,b)

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

  ax.set_xlabel('x-axis')
  ax.set_ylabel('y-axis')
  ax.set_zlabel('z-axis')
  # plt.show()
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

    return mesh, x ## no control points for now

if __name__ == "__main__":
  assert len(sys.argv) >= 4, "InputError:: %d inputs, expected 3. (directory, base_name, variable tag)"%(len(sys.argv)-1)
  directory = sys.argv[1]
  base_name = sys.argv[2]
  var_name  = sys.argv[3]
  # color     = sys.argv[4] if len(sys.argv) >= 5 else 'rainbow'

  # get the files
  files = [f for f in filter(lambda f: base_name in f, os.listdir(directory))]
  time_file = [f for f in filter(lambda f: f.endswith("_time.txt"), files)][0]
  field_files = [f for f in filter(lambda f: "_" + var_name + "_" in f , files)]

  field_files.sort(key=lambda x: int(re.findall("\d+",x)[0]))
  field_files.sort(key=lambda x: int(re.findall("\d+",x)[-1]))

  ## open all data
  mesh = []
  data = []
  fmn,fmx = 1e16, -1e16
  for f in field_files:
    m,x = open_file(directory + '/' + f)
    mesh.append(m)
    data.append(x)
    field = np.asarray(x)
    fmn = min(fmn,np.min(field[:,-1]))
    fmx = max(fmx,np.max(field[:,-1]))

  # plot the first file
  field_range = (fmn, fmx)
  print 'plotting step ', 0
  fig,ax = plot(mesh[0],data[0],field_range=field_range)
  last_step = 0
  class Updater:
    def __init__(self):
      self.last_step = 0

    def update(self,step):
      if self.last_step != step:
        print 'plotting step ', step
        ax.clear()
        plot(mesh[step],data[step],axis=ax,figure=fig,field_range=field_range)
        self.last_step = step

  up = Updater()
  ani = animate(plt.gcf(), up.update, np.arange(1, len(field_files)), interval=50, repeat=True) 
  plt.show()