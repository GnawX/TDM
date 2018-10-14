import numpy as np
from ase.io import read
from ase.io.cube import read_cube_data
from timeit import default_timer as timer

start = timer()

homo = read_cube_data('homo.cube')[0]
lumo = read_cube_data('lumo.cube')[0]


with open('homo.cube','r') as f:
     f.readline()
     f.readline()
     f.readline()
     line1 = f.readline().split()
     line2 = f.readline().split()
     line3 = f.readline().split()

n1,n2,n3 = int(line1[0]),int(line2[0]),int(line3[0])
x1,y1,z1 = [float(s) for s in line1[1:]]
x2,y2,z2 = [float(s) for s in line2[1:]]
x3,y3,z3 = [float(s) for s in line3[1:]]
vec = np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
da = np.linalg.norm(vec[0])
db = np.linalg.norm(vec[1])
dc = np.linalg.norm(vec[2])


a,b,c = np.meshgrid(range(n1),range(n2),range(n3),indexing='ij')
r = np.vstack((a.flatten(),b.flatten(),c.flatten())).T
r = np.dot(r,vec)

homo = homo.flatten()
lumo = lumo.flatten()
temp = np.multiply(homo,lumo)
dipole = np.dot(temp,r)*da*db*dc
print dipole,'e Bohr'

end = timer()
print 'elapsed time:', end - start, 'seconds.'
