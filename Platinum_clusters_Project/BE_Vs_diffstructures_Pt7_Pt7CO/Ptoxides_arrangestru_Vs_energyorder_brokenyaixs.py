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

def plot_conf1(ax, atoms, colorlenth,rot=False):
    colors = np.array([jmol_colors[atom.number] for atom in atoms])
    positions =atoms.get_positions()
    for i, atom in enumerate(atoms):
        if (atom.number ==78):
           colors[i] =[0.1, 0.6, 0.6]
        if (atom.number ==6):
           colors[i] =[0.1, 0.2, 0.9]
        if (atom.number ==8 and positions[i,2]>12.2):
           colors[i] =[128/255, 0/255, 128/255]
        if (positions[i,2]<12.7 ):
           colors[i] =[255/255, 255/255, 255/255]

    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Al' or a.symbol == 'O':
            if a.position[2] < 9.7:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)



#fig  = plt.figure(figsize=(14.0,10.50))
fig, [ax1,ax2,ax3,ax4] = plt.subplots(nrows=4, ncols=1,figsize=(13.0,12.50))
for i in range(2,15):
    ax1.plot(i,1)
    ax2.plot(i,1)
    ax3.plot(i,1)
    ax4.plot(i,1)
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
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax3.transAxes)  # switch to the bottom axes
ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax3.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
#---------------------------------------------------------------------------------------#
d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax3.transAxes, color='k', clip_on=False)
ax3.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

kwargs.update(transform=ax4.transAxes)  # switch to the bottom axes
ax4.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
ax4.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
#---------------------------------------------------------------------------------------#
ax1.set_xlim(2, 16)
ax2.set_xlim(2, 16)
ax3.set_xlim(2, 16)
ax4.set_xlim(2, 16)
ax1.set_yticks([])
ax2.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
#----------------------------------------------------------------------------------------#
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax4.spines['top'].set_visible(False)
ax1.set_ylabel('CO form', color='#0f0f0f')
ax2.set_ylabel('CO$_2$(attched) form', color='#8EBA42')
ax3.set_ylabel('CO$_2$ (detached) form', color='#FF00FF')
ax4.set_ylabel('CO$_3$ form', color='#FF0000')
fig.text(0.5, 0.04, 'No. of Oxygens(y) in Pt$_7$O$_{y}$ + CO', ha='center',fontsize=14)
fig.text(0.04, 0.5, 'Type of Structures', va='center', rotation='vertical',fontsize=14)
plt.subplots_adjust(wspace=0.05,hspace=0.02)
#plt.show()
#-----------------------------------------------------------------------------------------#
data=read(sys.argv[1]+'@:')
#-----------------------------------------------------------------------------------------#
global_ax = ax1
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
decrement =10
NO_of_oxygens =[2,4,6,8,10,12,14]
for j in range(0,len(data)):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    decrement = (decrement + 2)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-decrement or atom.index >=colorlenth*5]]
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
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    #if (j !=1):
    ax = plt.axes([xy[0]+0.0055, xy[1]-(dy)/15.50, 0.0892, 0.1])
    #if (j==1):
    #   ax =plt.axes([xy[0]+0.0066, xy[1]-(dy)/300.50, 0.067, 0.08])
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
    ax.plot(box_x, box_y, color='#0f0f0f',linewidth=5.0)
    plt.axes(ax)
 
    # 0 1                                                                                                                          
   # ax = plt.subplot()
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy) 
    ax = plt.axes([xy[0], xy[1]-(dy)/0.64, 0.1, 0.1])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.5])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
    name1 = "A$_{"+ str(j+1) + "}$"
    ax.text(0.05, 1.53, name1, transform=ax.transAxes,fontsize=10)
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
    ax.plot(box_x, box_y, color='#0f0f0f',linewidth=5.0)
    plt.axes(ax)
#---------------------------------- CO2 attached ---------------------------------------------#
data=read(sys.argv[2]+'@:')
#-----------------------------------------------------------------------------------------#
global_ax = ax2
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
decrement =10
NO_of_oxygens =[2,4,6,8,10,12,14]
for j in range(0,len(data)-1):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    decrement = (decrement + 2)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-decrement or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+8.0]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-8.0]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-7.0]
    #if (j ==4):
    #   del atoms[atoms.positions[:,0] >=centreofmass[0]+9.0]
    #   del atoms[atoms.positions[:,0] <= centreofmass[0]-7.10]
    #   del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
    #   del atoms[atoms.positions[:,1] <= centreofmass[1]-7.10]

    colorlenth = len(atoms)
    #view(atoms)
    cell = atoms.get_cell()

    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
   # if (j !=4):
    ax = plt.axes([xy[0]+0.0055, xy[1]-(dy)/15.50, 0.0892, 0.1])
    #if (j==4):
    #   ax = plt.axes([xy[0]+0.0056, xy[1]-(dy)/15.50, 0.089, 0.1])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)
    #if (j!=4):
    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    #if (j==4):
    #   ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0]) 
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
    ax.plot(box_x, box_y, color='#8EBA42',linewidth=5.0)
    plt.axes(ax)

    # 0 1                                                                                                                          
   # ax = plt.subplot()
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/0.64, 0.1, 0.1])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.5])
    #if (j==4):
    #   ax.set_xlim([centreofmass[0]-7.0, centreofmass[0]+8.0])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
    name1 = "B$_{"+ str(j+1) + "}$"
    ax.text(0.05, 1.53, name1, transform=ax.transAxes,fontsize=10)
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
    ax.plot(box_x, box_y, color='#8EBA42',linewidth=5.0)
    plt.axes(ax)
#--------------------------- CO2 detached    ----------------------------------------------#
data=read(sys.argv[3]+'@:')
#-----------------------------------------------------------------------------------------#
global_ax = ax3
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
decrement =10
NO_of_oxygens =[2,4,6,8,10,12,14]
for j in range(0,len(data)-1):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    decrement = (decrement + 2)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-decrement or atom.index >=colorlenth*5]]
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    if (j !=0 and j !=1 and j !=2):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+8.0]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-8.0]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-7.0]
    if (j==0):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+10.0]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-7.0]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+7.3]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-7.0]
    if (j ==1):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+8.0]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-8.0]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+8.3]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-8.0]
    if (j ==2):
       del atoms[atoms.positions[:,0] >=centreofmass[0]+7.50]
       del atoms[atoms.positions[:,0] <= centreofmass[0]-8.50]
       del atoms[atoms.positions[:,1] >= centreofmass[1]+5.80]
       del atoms[atoms.positions[:,1] <= centreofmass[1]-7.50]

    colorlenth = len(atoms)
   # if (j==2):
   #    view(atoms)
    cell = atoms.get_cell()

    # 0 0
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    #if (j !=1):
    ax = plt.axes([xy[0]+0.0055, xy[1]-(dy)/15.50, 0.0892, 0.1])
    #if (j==1):
    #   ax =plt.axes([xy[0]+0.0066, xy[1]-(dy)/300.50, 0.067, 0.08])
    img = atoms.copy()
    if (j!=6):
       plot_conf(ax, img,colorlenth)
    if (j==6): 
       plot_conf1(ax, img,colorlenth)

    ax.set_xlim([centreofmass[0]-7.50, centreofmass[0]+7.50])
    if (j==0):
       ax.set_xlim([centreofmass[0]-6.50, centreofmass[0]+8.50])
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
    ax.plot(box_x, box_y, color='#FF00FF',linewidth=5.0)
    plt.axes(ax)
    # 0 1                                                                                                                          
   # ax = plt.subplot()
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/0.64, 0.1, 0.1])
    cell = atoms.get_cell()
    img = atoms.copy()
    if (j!=6):
       plot_conf(ax, img,colorlenth, rot=True)
    if (j==6):
       plot_conf1(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.5])
    if (j==0):
       ax.set_xlim([centreofmass[0]-6.50, centreofmass[0]+8.50])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
    if (j==2):
       ax.set_ylim([centreofmass[1]-7.0, centreofmass[1]+6.50])
    if (j==1):
       ax.set_ylim([centreofmass[1]-5.5, centreofmass[1]+8.0])
    name1 = "C$_{"+ str(j+1) + "}$"
    ax.text(0.05, 1.53, name1, transform=ax.transAxes,fontsize=10)
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
    ax.plot(box_x, box_y, color='#FF00FF',linewidth=5.0)
    plt.axes(ax)
#------------------------------- CO3 form -------------------------------------------------#
data=read(sys.argv[4]+'@:')
#-----------------------------------------------------------------------------------------#
global_ax = ax4
transform = lambda x: fig.transFigure.inverted().transform(global_ax.transData.transform(x))
inverse = lambda x: global_ax.transData.inverted().transform(fig.transFigure.transform(x))
decrement =10
NO_of_oxygens =[2,4,6,8,10,12,14]
for j in range(0,len(data)-1):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    #print(colorlenth)
    decrement = (decrement + 2)
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-decrement or atom.index >=colorlenth*5]]
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
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    #if (j !=1):
    ax = plt.axes([xy[0]+0.0055, xy[1]-(dy)/15.50, 0.089, 0.1])
    #if (j==1):
    #   ax =plt.axes([xy[0]+0.0066, xy[1]-(dy)/300.50, 0.067, 0.08])
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
    ax.plot(box_x, box_y, color='#ff796c',linewidth=5.0)
    plt.axes(ax)
    # 0 1                                                                                                                          
    dy = (inverse((1, 0)) - inverse((1, -0.1)))[1]
    xy = transform((NO_of_oxygens[j], 1))
    print(dy, xy)
    ax = plt.axes([xy[0], xy[1]-(dy)/0.64, 0.1, 0.1])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    ax.set_xlim([centreofmass[0]-7.5, centreofmass[0]+7.5])
    ax.set_ylim([centreofmass[1]-6.5, centreofmass[1]+7.0])
    name1 = "D$_{"+ str(j+1) + "}$"
    ax.text(0.05, 1.53, name1, transform=ax.transAxes,fontsize=10)
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
    ax.plot(box_x, box_y, color='#ff796c',linewidth=5.0)
    plt.axes(ax)
#------------------------------------------------------------------------------------------#
ax1.set_xticks([])
ax2.set_xticks([])
ax3.set_xticks([])
name = sys.argv[5]
#savefig(name,format='png')
savefig(name,format='pdf')
show()
