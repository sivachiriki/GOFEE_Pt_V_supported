from ase.io import read, write,Trajectory
data=read('V2O5_H2O_TiO2_101_DFTlowlyingsorted.traj@:')
structures =[data[0],data[1],data[3],data[4],data[9]]
numbers=[0,1,3,4,9]
for i,a in enumerate(structures):
    p=numbers[i]
    name='V2O5_H2O_TiO2_101_DFT_{}.xyz'.format(p+1)
    write(name,a)
