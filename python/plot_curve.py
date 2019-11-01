import sys
import matplotlib.pyplot as plt

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
      plt.scatter(p[0],p[1],color='green')

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
    if len(sys.argv) > 2:
      vec_data = open_vector_data(sys.argv[2])
      plot(m,x,cpt,vector_data=vec_data)
    else:
      plot(m,x,cpt)