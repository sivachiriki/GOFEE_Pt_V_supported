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

plt.rcParams['savefig.facecolor'] = "1.0"
#######################################
fig = pyplot.figure(figsize=(6, 5))
data = read(sys.argv[1]+'@:')
gs1 = gridspec.GridSpec(1, 1)
ax1 = fig.add_subplot(gs1[0])
n_images=len(data)
atoms=data[0]
atoms1=atoms
atoms.rotate('y',np.radians(-90),rotate_cell=True)
atoms.rotate('z',np.radians(-90),rotate_cell=True)
# Add the atoms to the plot as circles.
for atom in atoms:
    color = jmol_colors[atom.number]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax1.add_patch(circle)
ax1.axis('equal')
ax1.set_xticks([])
ax1.set_yticks([])
ax1.axis('off')

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
gs1.tight_layout(fig, rect=[0.0, 0.2, 0.4, 0.5]) # (left, bottom, right, top)
##########################################################################################
data = read(sys.argv[2]+'@:')
gs2 = gridspec.GridSpec(1, 1)
ax2 = fig.add_subplot(gs2[0])
n_images=len(data)
atoms=data[0]
atoms1=atoms
atoms.rotate('y',np.radians(-90),rotate_cell=True)
atoms.rotate('z',np.radians(-90),rotate_cell=True)
# Add the atoms to the plot as circles.
for atom in atoms:
    color = jmol_colors[atom.number]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax2.add_patch(circle)
ax2.axis('equal')
ax2.set_xticks([])
ax2.set_yticks([])
ax2.axis('off')

xlim = ax2.get_xlim()
ylim = ax2.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax2.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax2.plot(box_x, box_y, color='black',linewidth=0.5)
gs2.tight_layout(fig, rect=[0.31, 0.1, 0.60, 0.6]) # (left, bottom, right, top)
#############################################################################################
data = read(sys.argv[3]+'@:')
gs3 = gridspec.GridSpec(1, 1)
ax3 = fig.add_subplot(gs3[0])
n_images=len(data)
atoms=data[0]
atoms1=atoms
atoms.rotate('y',np.radians(-90),rotate_cell=True)
atoms.rotate('z',np.radians(-90),rotate_cell=True)
# Add the atoms to the plot as circles.
for atom in atoms:
    color = jmol_colors[atom.number]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax3.add_patch(circle)
ax3.axis('equal')
ax3.set_xticks([])
ax3.set_yticks([])
ax3.axis('off')

xlim = ax3.get_xlim()
ylim = ax3.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax3.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax3.plot(box_x, box_y, color='black',linewidth=0.5)
gs3.tight_layout(fig, rect=[0.61, 0.1, 0.91, 0.6],h_pad=0.2,w_pad=0.2) # (left, bottom, right, top)
############################################################################################
data = read(sys.argv[3]+'@:')
gs4 = gridspec.GridSpec(1, 1)
ax4 = fig.add_subplot(gs4[0])
n_images=len(data)
atoms=data[0]
atoms1=atoms
#atoms.rotate('y',np.radians(-90),rotate_cell=True)
#atoms.rotate('z',np.radians(-90),rotate_cell=True)
# Add the atoms to the plot as circles.
for atom in atoms:
    color = jmol_colors[atom.number]
    radius = radii[atom.number]
    circle = Circle((atom.x, atom.y), radius*0.90, facecolor=color,
                            edgecolor='k',linewidth=0.5,alpha=1.0,zorder=2)
    ax4.add_patch(circle)
ax4.axis('equal')
ax4.set_xticks([])
ax4.set_yticks([])
ax4.axis('off')

xlim = ax4.get_xlim()
ylim = ax4.get_ylim()
box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
ax4.add_patch(
     patches.Rectangle(
        (box_x[0],box_y[0]),
        xlim[1]-xlim[0],
        ylim[1]-ylim[0],
        fill=True,facecolor='white', clip_on=False, zorder=1      # remove background
     ) )
ax4.plot(box_x, box_y, color='black',linewidth=0.5)
gs4.tight_layout(fig, rect=[0.81, 0.1, 0.6, 0.6],h_pad=0.2,w_pad=0.3) # (left, bottom, right, top)
#############################################################################################
#pyplot.subplots_adjust(top=0.96, bottom=0.00, left=0.18, right=0.95, hspace=0.07,wspace=0.35)
pyplot.show()
fig.savefig('methodology_surface_Al2O3.pdf')
