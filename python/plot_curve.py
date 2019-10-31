import sys
import matplotlib.pyplot as plt

def plot(mesh, verts, cpts):
  for el in mesh:
    a = verts[el[0]]
    b = verts[el[1]]
    plt.plot([a[0], b[0]],[a[1], b[1]], 'k')

  for p in cpts:
    plt.scatter(p[0],p[1],color='red',marker='s')

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

if __name__ == "__main__":
  if len(sys.argv) > 1:
    m,x,cpt = open_file(sys.argv[1])
    plot(m,x,cpt)