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
                     # alpha = alp[ia],
                      #zorder = 1 - apos[1]/1000
                      zorder =1.0
                      )
        ax.add_patch(circ)



def plot_conf(ax, atoms, colorlenth,rot=False):
    colors = np.array([jmol_colors[atom.number] for atom in atoms])
    positions =atoms.get_positions()
    for i, atom in enumerate(atoms):
        if (positions[i,0]<=12.5 or positions[i,0]>=24.75):
           colors[i] =[255/255, 255/255, 255/255]

    alp = [None] * colors.shape[0]
    #for i,a in enumerate(atoms):
    #    if a.symbol == 'Al' or a.symbol == 'O':
    #        if a.position[2] < 9.7:
    #            alp[i] = 1.0

    if rot:
       # atoms.rotate('x',pi/2)
       atoms.rotate('y',np.radians(-90),rotate_cell=True)
       atoms.rotate('z',np.radians(-90),rotate_cell=True)

    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)

#-----------------------------------------------------------#
fig = plt.figure(figsize=(6.0,5))
outer = gridspec.GridSpec(1, 1, wspace=0.04, hspace=0.2)

color_lib = ['#00FF00','#377eb8','#4daf4a','#00FFFF','#a65628','#FF0000','#0000FF', '#FF00FF','#FFFF00','#000000']
#---------------------- Pt7 clusters -------------------------------------#
data=read(sys.argv[1]+'@:')
for j in range(0,len(data)):
    #inner = gridspec.GridSpecFromSubplotSpec(2, 1,subplot_spec=outer[j], wspace=0.00, hspace=0.0, height_ratios=[6.86,9.9])
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(1,3,1)
    print(colorlenth)
   # write('newimage.traj',atoms)
    a=atoms
    #view(atoms)
    del atoms[atoms.positions[:,1] <=6.0]
    del atoms[atoms.positions[:,1] >=28.0]
    #view(atoms)
    colorlenth = len(atoms)
    cell = atoms.get_cell()
    # 0 0
    ax = plt.Subplot(fig, outer[0])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth,rot=False)

    #ax.set_xlim([0.20, 0.60])
    #ax.set_ylim([0.0, 20.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    fig.add_subplot(ax)
    #----------------- drawing box -------------------------------#
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    print(xlim)
    print(ylim)
    box_x = [xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]]
    box_y =[ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]]
    ax.add_patch(
        patches.Rectangle(
           (box_x[0],box_y[0]),
           xlim[1]-xlim[0],
           ylim[1]-ylim[0],
           fill=True,facecolor='white', clip_on=False,zorder =0.8) )
    ax.plot(box_x, box_y, color='blue',linewidth=5.0)

name = sys.argv[2]
name =name
savefig(name,bbox_inches='tight')
show()
exit()
