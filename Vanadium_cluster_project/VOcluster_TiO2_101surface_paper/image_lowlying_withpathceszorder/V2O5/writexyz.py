from ase.io import read, write,Trajectory
data=read('V2O5_TiO2_101_DFTlowlyingsorted.traj@:')
structures =[data[0],data[8],data[10],data[16]]
numbers=[0,8,10,16]
for i,a in enumerate(structures):
    p=numbers[i]
    name='V2O5_TiO2_101_DFT_{}.xyz'.format(p+1)
    write(name,a)
