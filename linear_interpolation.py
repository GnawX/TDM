from ase import atoms
from ase import io

nimgs = 7

initial = io.read('start.xyz')
final = io.read('end.xyz')

images = [initial]
images += [initial.copy() for i in range(nimgs-2)]
images += [final]
images += [final.copy() for i in range(nimgs-1)]

pos0 = initial.get_positions()
pos1 = final.get_positions()
dist = (pos1 - pos0)/(nimgs - 1)

for i in range(nimgs*2-1):
    images[i].set_positions(pos0 + i*dist)
    fil = 'image'+ str(i+1) + '.xyz'
    io.write(fil, images[i])
