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
a = np.array([x1,y1,z1])
b = np.array([x2,y2,z2])
c = np.array([x3,y3,z3])
da = np.linalg.norm(a)
db = np.linalg.norm(b)
dc = np.linalg.norm(c)

homo = homo.flatten()
lumo = lumo.flatten()
overlap = np.dot(homo,lumo)*da*db*dc
print overlap

end = timer()
print 'elapsed time:', end - start, 'seconds.'
