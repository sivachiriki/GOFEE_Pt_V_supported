from ase.io import read, write
from ase.visualize import view

data=read('TiO2_101surface_3by1_2layer.traj@:')
image =data[0]
cell=image.get_cell()
print(image.get_cell())
write('TiO2_101surface_3by1_2layer.xyz',image,vec_cell=True)
