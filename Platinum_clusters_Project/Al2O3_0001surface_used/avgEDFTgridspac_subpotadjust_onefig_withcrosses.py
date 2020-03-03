from ase.ga.data import DataConnection
from optparse import OptionParser
from ase.io.trajectory import Trajectory
from ase.io import read,write
from ase import Atoms
import numpy as np
import os
import sys
import matplotlib.gridspec as gridspec
import matplotlib as plt
from matplotlib import pyplot
from matplotlib.pyplot import *
from ase.ga.generate_feavec_correplot_raddes_Ti import drawcorrelationplot
from matplotlib.patches import ConnectionPatch
from ase.data import covalent_radii as radii
from ase.data.colors import jmol_colors
from matplotlib.patches import Circle
import matplotlib.patches as patches
from ase.visualize import view

plt.rcParams['savefig.facecolor'] = "1.0"
#######################################
fig = pyplot.figure(figsize=(8, 4))
data = read(sys.argv[1]+'@:')
gs1 = gridspec.GridSpec(1, 3)
ax1 = fig.add_subplot(gs1[0])
n_images=len(data)
atoms=data[0]*(1,3,1)
#write('131cell.traj',atoms)
atoms1=atoms
atoms.rotate('y',np.radians(-90),rotate_cell=True)
atoms.rotate('z',np.radians(-90),rotate_cell=True)

#del atoms[atoms.positions[:,0] <=8.0]
#del atoms[atoms.positions[:,0] >=28.0]
# Add the atoms to the plot as circles.
positions =atoms.get_positions()
for i,atom in enumerate(atoms):
    color = jmol_colors[atom.number]
    if (positions[i,0]<=12.32 or positions[i,0]>=24.73):
       color =[255/255, 255/255, 255/255]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax1.add_patch(circle)
ax1.set_xlim([8.0,28.0])
ax1.set_ylim([5.0,22.0])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set(aspect=1)
print(ax1.get_xlim(),ax1.get_ylim())
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax1.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax1.plot(box_x, box_y, color='black',linewidth=0.5)
#------------------ blue box ----------------------------------------------#
print(box_x[0],box_y[0],xlim[1]-xlim[0],ylim[1]-ylim[0])
xlim =[11.3,24.5]
ylim =[11.8,20.5]
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax1.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=0.8      # remove background
     ) )
print(box_x[0],box_y[0],xlim[1]-xlim[0],ylim[1]-ylim[0])
ax1.plot(box_x, box_y, color='blue',linewidth=2.0)
ax1.text(8.1,22.5,'a)',fontsize=14)
#------------------------------------------------------------------------#
ax1 = fig.add_subplot(gs1[1])
data = read(sys.argv[2]+'@:')
n_images=len(data)
atoms=data[0]*(1,3,1)
write('131cell.traj',atoms)
atoms1=atoms
atoms.rotate('y',np.radians(-90),rotate_cell=True)
atoms.rotate('z',np.radians(-90),rotate_cell=True)

#del atoms[atoms.positions[:,0] <=8.0]
#del atoms[atoms.positions[:,0] >=28.0]
# Add the atoms to the plot as circles.
positions =atoms.get_positions()
for i,atom in enumerate(atoms):
    color = jmol_colors[atom.number]
    if (atom.number ==78):
       color =[0.1, 0.6, 0.6]
    if (positions[i,0]<=12.32 or positions[i,0]>=24.73):
       color =[255/255, 255/255, 255/255]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax1.add_patch(circle)
#-------------------------- make crosses on top of the circle -------------------#
print(positions[237,0],positions[237,1],positions[237,2])
colorcode_cross =[96/255,96/255,96/255]
for i,atom in enumerate(atoms):
    if (positions[i,1]<8.7 and atom.number ==13):
        mymarker = ax1.scatter(positions[i,0],positions[i,1], s=60, c=colorcode_cross,linewidth=0.1, marker='x', zorder=2.3,clip_on=True)
        ax1.add_artist(mymarker)
    if (positions[i,1]<8.7 and atom.number ==8):
        mymarker = ax1.scatter(positions[i,0],positions[i,1], s=10, c=colorcode_cross,linewidth=0.05, marker='x', zorder=2.3,clip_on=True)
        ax1.add_artist(mymarker)
#--------------------------------------------------------------------------------#
ax1.set_xlim([8.0,28.0])
ax1.set_ylim([5.0,22.0])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set(aspect=1)
print(ax1.get_xlim(),ax1.get_ylim())
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax1.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax1.plot(box_x, box_y, color='black',linewidth=0.5)
ax1.text(8.1,22.5,'b)',fontsize=14)
#-------------------------------------------------------------#
ax1 = fig.add_subplot(gs1[2])
data = read(sys.argv[2]+'@:')
n_images=len(data)
atoms=data[0]*(3,3,1)
atoms1=atoms
#atoms.rotate('y',np.radians(-90),rotate_cell=True)
#atoms.rotate('z',np.radians(-90),rotate_cell=True)

#del atoms[atoms.positions[:,0] <=8.0]
#del atoms[atoms.positions[:,0] >=28.0]
#write('131cell.traj',atoms)
#view(atoms)
# Add the atoms to the plot as circles.
positions =atoms.get_positions()
for i,atom in enumerate(atoms):
    color = jmol_colors[atom.number]
    if (atom.number ==78):
       color =[0.1, 0.6, 0.6]
    if (i<567 or i>709):
       color =[255/255, 255/255, 255/255]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax1.add_patch(circle)
ax1.set_xlim([-2.0,22.0])
ax1.set_ylim([7.0,28.0])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set(aspect=1)
print(ax1.get_xlim(),ax1.get_ylim())
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax1.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax1.plot(box_x, box_y, color='black',linewidth=0.5)
ax1.text(-2.0,28.50,'c)',fontsize=14)
#------------------------------------------------------------#
pyplot.show()
fig.savefig('methodology_surface_Al2O3_withcrosses.pdf')
