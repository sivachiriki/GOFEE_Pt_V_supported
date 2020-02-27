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
import matplotlib.cbook as cbook
import matplotlib.dates as dates
import matplotlib.ticker as ticker

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
           colors[i] =[0.1, 0.2, 0.9]
        if (atom.number ==8 and positions[i,2]>12.2):
           colors[i] =[128/255, 0/255, 128/255]

    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Al' or a.symbol == 'O':
            if a.position[2] < 9.7:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)



#fig = plt.figure(figsize=(10.0,10.50))
fig, [ax1,ax2] = plt.subplots(nrows=2, ncols=1,figsize=(10.0,10.50))
for i in range(0,8):
    ax1.plot(i,0.0)
    ax2.plot(i,-2.0)
#---------------------------------------------------------------------------------------#
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
#---------------------------------------------------------------------------------------#
ax1.set_xlim(0.4, 7.5)
ax2.set_xlim(0.4, 7.5)
ax1.yaxis.set_ticks(np.arange(-0.30, 0.9, 0.1))
ax2.yaxis.set_ticks(np.arange(-3.0, -1.90, 0.1))
#----------------------------------------------------------------------------------------#
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
fig.text(0.5, 0.04, 'Lowest Isomers found for Pt$_7$ and Pt$_7$CO with GOFEE', ha='center',fontsize=14)
fig.text(0.04, 0.5, 'Stability of Isomers (eV)', va='center', rotation='vertical',fontsize=14)
plt.subplots_adjust(wspace=0.02,hspace=0.07)
#plt.show()
#-----------------------------------------------------------------------------------------#
data=read(sys.argv[1]+'@:')
energydif =np.zeros(len(data))
for j in range(len(data)):
    GM_energy = data[0].get_potential_energy()
    energydif[j] = (data[j].get_potential_energy() - GM_energy)
    print(j,energydif[j])
global_ax = ax1
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
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.8]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]
  
    colorlenth = len(atoms) 
    
    cell = atoms.get_cell()
   
    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j+0.5, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/30.50, 0.10, 0.10])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-6.50, centreofmass[0]+7.50])
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
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j+0.5, energydif[j]))
    print(dy, xy) 
    ax = plt.axes([xy[0], xy[1]-(dy)/3.40, 0.10, 0.10])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.5])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
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

#-------------------- set 1 Pt7+CO on surface ---------------#
CO =read('CO_molecule_DFTrelaxed.traj@:')
E_CO = CO[0].get_potential_energy()
data=read(sys.argv[2]+'@:')
energydif =np.zeros(len(data))
for j in range(len(data)):
    energydif[j] = (data[j].get_potential_energy()- GM_energy -E_CO)
    print(j,energydif[j])

global_ax = ax2
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
for j in range(0,7):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-10 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.8]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]

    colorlenth = len(atoms)

    cell = atoms.get_cell()
    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j+0.5, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/15.50, 0.10, 0.10])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    ax.set_ylim([10.7, 20.0])
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
    # 0 1                                                                                                                          
   # ax = plt.subplot()
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j+0.5, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/3.15, 0.10, 0.10])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.5])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
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

ax1.set_xticks([])

name = sys.argv[3]
savefig(name,bbox_inches='tight')
show()
