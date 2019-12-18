from ase.io import read, write
from ase.visualize import view

data=read('anataseTi24O48_101surface_optPBEesben_1by2.traj@:')
image =data[0]
image = image * (3,3,1)
cell=image.get_cell()
print(image.get_cell())
write('newimage.xyz',image,vec_cell=True)
