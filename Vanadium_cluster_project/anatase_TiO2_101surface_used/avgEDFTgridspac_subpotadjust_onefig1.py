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

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)

def plot_atoms(ax, atoms, xyz, acols, alp, z):

    ecols = [[0, 0, 0] for col in atoms]

    indices = range(len(atoms))
    for ia in indices:
        acol = acols[ia]
        ecol = ecols[ia]
        arad = aradii[atoms[ia].number]
        apos = atoms[ia].position
        eps = arad

        circ = Circle([apos[xyz[0]], apos[xyz[1]]],
                      fc = acol,
                      ec = ecol,
                      radius = arad,
                      lw = 0.5,
                      alpha = alp[ia],
                      zorder = 1 - apos[1]/1000
                      )
        ax.add_patch(circ)



def plot_conf(ax, atoms, colorlenth,rot=False):
    colors = np.array([jmol_colors[atom.number] for atom in atoms])
    positions =atoms.get_positions()
    for i, atom in enumerate(atoms):
        if (atom.number ==23):
           colors[i] =[76/255, 153/255, 0/255]
        if (atom.number ==8 and positions[i,2] >=16.0):
           colors[i] =[153/255, 0/255, 0/255]
        if (atom.number ==1):
           colors[i] =[255/255, 255/255, 255/255]

    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Ti' or a.symbol == 'O':
            if a.position[2] < 10.50:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)


plt.rcParams['savefig.facecolor'] = "1.0"
#######################################
fig = pyplot.figure(figsize=(8, 8))
data = read(sys.argv[1]+'@:')
gs1 = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs1[0])
#-------------------------------------------------------------#
global_ax = ax
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
for j in range(0,len(data)):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-8 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+11.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-11.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+11.8]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-11.10]

    colorlenth = len(atoms)

    cell = atoms.get_cell()

    # 0 0
  #  dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
  #  xy = transform((j+0.8, energydif[j]))
  #  print(dy, xy)
  #  ax = plt.axes([xy[0], xy[1]-(dy)/6.50, 0.170, 0.170])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-11.50, centreofmass[0]+11.50])
    ax.set_ylim([10.7, 20.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    #----------------- drawing box -------------------------------#
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
    box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
    ax.add_patch(
        patches.Rectangle(
           (box_x[0],box_y[0]),
           xlim[1]-xlim[0],
           ylim[1]-ylim[0],
           fill=True,facecolor='white', clip_on=False,zorder =0.8) )
    ax.plot(box_x, box_y, color='blue',linewidth=5.0)
    plt.axes(ax)

    # 0 1                                                                                                                          
   # ax = plt.subplot()
   # dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    #xy = transform((j+0.8, energydif[j]))
   # print(dy, xy)
   # ax = plt.axes([xy[0], xy[1]-(dy)/2.6, 0.170, 0.170])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-11.5, centreofmass[0]+11.5])
    ax.set_ylim([centreofmass[1]-10.5, centreofmass[1]+11.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
    box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
    ax.add_patch(
        patches.Rectangle(
           (box_x[0],box_y[0]),
           xlim[1]-xlim[0],
           ylim[1]-ylim[0],
           fill=True,facecolor='white', clip_on=False,zorder =0.8) )
    ax.plot(box_x, box_y, color='blue',linewidth=5.0)
    plt.axes(ax)

#-------------------------------------------------------------#
pyplot.show()
fig.savefig('TiO2_anatase_101surface.png')
