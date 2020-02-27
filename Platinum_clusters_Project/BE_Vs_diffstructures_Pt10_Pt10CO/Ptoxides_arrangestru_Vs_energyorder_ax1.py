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



#fig  = plt.figure(figsize=(14.0,10.50))
fig, [ax1,ax2] = plt.subplots(nrows=1, ncols=2,figsize=(14.0,10.50),gridspec_kw={'width_ratios': [0.9, 1.1]})

for i in range(0,5):
    ax1.plot(i,0.0)
for i in range(4,10):
    ax2.plot(i,-2.0)

d = .01 # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d,1+d), (-d,+d), **kwargs)
ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d,+d), (1-d,1+d), **kwargs)
ax2.plot((-d,+d), (-d,+d), **kwargs)
#-----------------------------------------------------------------------------------------#
data=read(sys.argv[1]+'@:')
energydif =np.zeros(len(data))
for j in range(0,4):
    GM_energy = data[5].get_potential_energy()
    energydif[j] = (data[j].get_potential_energy() - GM_energy)
#    ax1.plot(j,energydif[j])
    print(j,energydif[j])
#ax1.set_xticks(np.arange(1, 5))
#ax1.set_xlim(0.5,4.5)
ax1.set_yticks(np.arange(-4.0, 3.50, 0.5))
ax1.set_xlabel('Structures with CO$_3$ formation',color='#FF0000')
#----------------------------------------------------------------------------------------#
#for j in range(4,10):
#    GM_energy = data[5].get_potential_energy()
#    energydif[j] = (data[j].get_potential_energy() - GM_energy)
#    ax2.plot(j,energydif[j])
#    print(j,energydif[j])
#ax2.set_xticks(np.arange(4.0,10.0))
#ax2.set_xlim(4.5,9.5)
#ax2.set_yticks([])
ax2.set_yticks(np.arange(-4.0, 3.50, 0.5))
ax2.set_xlabel('Structures with CO$_2$ formation',color='#FF00FF')
#----------------------------------------------------------------------------------------#
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
fig.text(0.5, 0.04, 'Lowest Isomer found for Pt$_7$O$_{10}$ and Pt$_7$O$_{10}$CO with GOFEE', ha='center',fontsize=18)
fig.text(0.04, 0.5, 'Stability of Isomers (eV)', va='center', rotation='vertical',fontsize=18)
plt.subplots_adjust(wspace=0.05)
#plt.show()
#exit()
#-----------------------------------------------------------------------------------------#
global_ax = ax1
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
for j in range(0,4):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-18 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.0]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.0]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.0]
  
    colorlenth = len(atoms) 
    
    cell = atoms.get_cell()
   
    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.007, xy[1]-(dy)/300.50, 0.066, 0.08])
    #ax = plt.axes([xy[0], xy[1]-(dy)/300.50, 0.08, 0.08])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
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
    xy = transform((j, energydif[j]))
    print(dy, xy) 
    ax = plt.axes([xy[0], xy[1]-(dy)/13.3, 0.08, 0.08])
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
global_ax = ax1
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
for j in range(0,4):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-20 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]
    if (j ==6):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+9.0]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-7.10]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]

    colorlenth = len(atoms)

    cell = atoms.get_cell()
    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.007, xy[1]-(dy)/300.0, 0.066, 0.08])
    #ax = plt.axes([xy[0], xy[1]-(dy)/300.0, 0.08, 0.08])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    if (j==6):
       ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0])
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
    xy = transform((j, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/13.3, 0.08, 0.08])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.5])
    if (j ==6):
       ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0])
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

#--------------------------------------------------------------------------------------#
data=read(sys.argv[3]+'@:')
for j in range(0,len(data)):
    GM_energy = data[5].get_potential_energy()
    energydif[j] = (data[j].get_potential_energy() - GM_energy)
    print(j,energydif[j])
#-----------------------------------------------------------------------------------------#
global_ax = ax2
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
for j in range(4,9):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-18 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.0]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.0]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.0]

    colorlenth = len(atoms)

    cell = atoms.get_cell()

    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j-0.2, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.032, xy[1]-(dy)/300.50, 0.066, 0.08])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
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
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j-0.2, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.025, xy[1]-(dy)/13.3, 0.08, 0.08])
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
data=read(sys.argv[4]+'@:')
energydif =np.zeros(len(data))
for j in range(len(data)):
    energydif[j] = (data[j].get_potential_energy()- GM_energy -E_CO)
    print(j,energydif[j])
for j in range(4,9):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-20 or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.10]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.10]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]
    if (j ==6):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+9.0]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-7.10]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]

    colorlenth = len(atoms)

    cell = atoms.get_cell()
    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j-0.2, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.032, xy[1]-(dy)/300.50, 0.066, 0.08])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    if (j==6):
       ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0])
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
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((j-0.2, energydif[j]))
    print(dy, xy)
    ax = plt.axes([xy[0]+0.025, xy[1]-(dy)/13.3, 0.08, 0.08])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.5])
    if (j ==6):
       ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0])
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
#--------------------------------------------------------------------------------------#
#xint = range(math.ceil(0.5), math.ceil(4.5))
ax1.set_xlim(0.5,4.5)
ax1.set_xticks([1,2,3,4])
#ax1.set_xticks(xint)
#ax1.set_xticks(np.arange(0.5,6,1))
ax2.set_xlim(4,9.8)
ax2.set_xticks([5,6,7,8,9])
ax2.set_yticks([])
#--------------------------------------#
name = sys.argv[5]
#savefig(name,bbox_inches='tight',format='png')
savefig(name,format='pdf')
show()
