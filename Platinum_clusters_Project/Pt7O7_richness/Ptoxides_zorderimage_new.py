from __future__ import division
import matplotlib
#matplotlib.use('Agg')  # Can also use 'tkagg' or 'webagg'
#from plot_neb_tio2 import *
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import matplotlib.patches as patches
from math import ceil, floor
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.visualize import view
import matplotlib.patches as mpatches
from ase.data.colors import jmol_colors
from decimal import Decimal 
from pylab import *
from ase.data import covalent_radii as aradii
from matplotlib.patches import Circle
from math import atan2,pi
import matplotlib.gridspec as gridspec

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
        if (atom.number ==78):
           colors[i] =[0.1, 0.6, 0.6]
        if (atom.number ==6):
           colors[i] =[0.0, 0.0, 0.0]
        if (atom.number ==8 and positions[i,2]>12.2):
           colors[i] =[128/255, 0/255, 128/255]
           #colors[i] =[0.0, 128/255,0.0]
     #   if (atom.number ==8 and i >=colorlenth*5-8):
     #      colors[i] =[102/255, 0/255, 0/255]
       # if (atom.number ==8 and i >= 135+colorlenth*2 and i <colorlenth*3 ):
       #    colors[i] =[102/255, 0/255, 0/255]
       # if (atom.number ==8 and i >= 135+colorlenth*3 and i <colorlenth*4 ):
       #    colors[i] =[102/255, 0/255, 0/255]
      #  if (atom.number ==8 and i >= 135+colorlenth*4 and i <colorlenth*5 ):
      #     colors[i] =[102/255, 0/255, 0/255]
      #  if (atom.number ==8 and i >= 135+colorlenth*5 and i <colorlenth*6 ):
      #     colors[i] =[102/255, 0/255, 0/255]

    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Al' or a.symbol == 'O':
            if a.position[2] < 9.7:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)

#-----------------------------------------------------------#
fig = plt.figure(figsize=(13.0,10.5))
outer = gridspec.GridSpec(4, 9, wspace=0.04, hspace=0.2)

color_lib = ['#00FF00','#377eb8','#4daf4a','#00FFFF','#a65628','#FF0000','#0000FF', '#FF00FF','#FFFF00','#000000']
#---------------------- Pt7 clusters -------------------------------------#
data=read(sys.argv[1]+'@:')
energydif =np.zeros(len(data))
for j in range(len(data)):
    GM_energy = data[0].get_potential_energy()
    energydif[j] = (data[j].get_potential_energy() - GM_energy)

for j in range(0,len(data)):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=outer[j], wspace=0.00, hspace=0.0, height_ratios=[6.86,9.9])
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    print(colorlenth)
   # write('newimage.traj',atoms)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-15 or atom.index >=colorlenth*5]]
    #view(atoms)
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.8]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]

    colorlenth = len(atoms)
    #view(atoms)
    cell = atoms.get_cell()
    # 0 0
    ax = plt.Subplot(fig, inner[0])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    ax.set_ylim([10.7, 20.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    fig.add_subplot(ax)
    if (j ==0):
       name2 ='Pt$_7$O$_7$ Lowest Isomers'
       ax.set_title(name2)

    #----------------- drawing box -------------------------------#
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    #print(xlim)
    #print(ylim)
    box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
    box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
    ax.add_patch(
        patches.Rectangle(
           (box_x[0],box_y[0]),
           xlim[1]-xlim[0],
           ylim[1]-ylim[0],
           fill=True,facecolor='white', clip_on=False,zorder =0.8) )
    ax.plot(box_x, box_y, color='blue',linewidth=5.0)

    # 0 1                                                                                                                          
    ax = plt.Subplot(fig, inner[1])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)

    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.50])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
    name ='$\Delta E = {:3.3f}$ eV'.format(energydif[j])
    ax.text(0.05, -0.14, name, transform=ax.transAxes,fontsize=10)
    name1 = "S$_{"+ str(j+1) + "}$"
    ax.text(0.05, 1.6, name1, transform=ax.transAxes,fontsize=10)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    #----------------- drawing box -------------------------------#
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    #print(xlim)
    #print(ylim)
    box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
    box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
    ax.add_patch(
        patches.Rectangle(
           (box_x[0],box_y[0]),
           xlim[1]-xlim[0],
           ylim[1]-ylim[0],
           fill=True,facecolor='white', clip_on=False,zorder =0.8) )
    ax.plot(box_x, box_y, color='blue',linewidth=5.0)

    fig.add_subplot(ax)
name = sys.argv[2]
name =name
savefig(name,bbox_inches='tight')
show()
exit()
