from ase.io import read, write
from ase.visualize import view

data=read('2VO3Htrans_initial_TiO21013by1cell2layer_DFTrelxd.traj@:')
image =data[0]
cell=image.get_cell()
print(image.get_cell())
write('2VO3Htrans_final_TiO21013by1cell2layer.xyz',image,vec_cell=True)
