from ase.io import read, write,Trajectory
from ase.visualize import view
from ase.io.trajectory import TrajectoryWriter

data=read('2VO3Htrans_starting_image_neb_cell4by1.traj@:')
image =data[0]
cell=image.get_cell()
data=read('modifiednebimages.traj@:')
traj = Trajectory('2VO3Htransmodifiednebimages.traj','w')
for i in range(len(data)):
    data[i].set_cell(cell)
    data[i].center()
    print(data[i].get_cell())
    traj.write(data[i])
