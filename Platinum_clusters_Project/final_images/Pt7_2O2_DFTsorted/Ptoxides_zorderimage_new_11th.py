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


#matplotlib.rcParams['xtick.direction'] = 'out'
#matplotlib.rcParams['ytick.direction'] = 'out'
#matplotlib.rc('font',**{'family':'sans-serif',
#                        'sans-serif':['Helvetica'],
#                        'size':14})
#matplotlib.rc('text',usetex=True)

#matplotlib.rcParams['text.latex.unicode']=False
#matplotlib.rcParams['text.latex.preamble']=['\usepackage{bm}']
#matplotlib.rcParams['text.latex.preamble']=['\usepackage{xfrac}']
#matplotlib.rcParams['mathtext.default'] = 'regular'
#matplotlib.rcParams['ps.usedistiller'] = 'xpdf'

matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)



def plot_atoms(ax, atoms, xyz, acols, alp, z):

    ecols = [[0, 0, 0] for col in atoms]

    indices = range(len(atoms))
    for ia in indices:
        acol = acols[ia]
        ecol = ecols[ia]
       # if atoms[ia].symbol == 'Ti':
       #     arad = aradii[atoms[ia].number] #* 0.9 * 0.5
       # else:
        arad = aradii[atoms[ia].number] #* 0.9 
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
           colors[i] =[102/255, 0/255, 0/255]
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



data=read(sys.argv[1]+'@:')
energydif =np.zeros(len(data))
for j in range(len(data)):
    GM_energy = data[0].get_potential_energy()
    energydif[j] = (data[j].get_potential_energy() - GM_energy)
    #print('{:3.3f}'.format(energydif[j]))
for j in range(46,len(data)):
    atoms = data[j]
    colorlenth = len(atoms)
    atoms =atoms*(3,3,1)
    print(colorlenth)
   # write('newimage.traj',atoms)
    #exit()
    a=atoms
    del atoms[[atom.index for atom in atoms if atom.index <=colorlenth*5-8 or atom.index >=colorlenth*5]]
    #view(atoms)
    centreofmass = a.get_center_of_mass()
    atoms = data[j]*(3,3,1)
    a=atoms
    del atoms[atoms.positions[:,0] >=centreofmass[0]+9.0]
    del atoms[atoms.positions[:,0] <= centreofmass[0]-9.0]
    del atoms[atoms.positions[:,1] >= centreofmass[1]+9.0]
    del atoms[atoms.positions[:,1] <= centreofmass[1]-8.30]
  
    colorlenth = len(atoms) 
    #view(atoms)
    #exit() 
    plt.figure(figsize=(4.0,5.0))
    gs = gridspec.GridSpec(2, 1,
                           height_ratios=[6.86,11.80])
    
    cell = atoms.get_cell()
    # 0 0
    ax = plt.subplot(gs[0, 0])
    img = atoms.copy()
    plot_conf(ax, img,colorlenth)
    
    ax.set_xlim([centreofmass[0]-8.0, centreofmass[0]+8.0])
    ax.set_ylim([10.7, 20.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    
    # 0 1                                                                                                                          
    ax = plt.subplot(gs[1, 0])
    cell = atoms.get_cell()
    img = atoms.copy()
    plot_conf(ax, img,colorlenth, rot=True)
    
    ax.set_xlim([centreofmass[0]-8.0, centreofmass[0]+8.0])
    ax.set_ylim([centreofmass[1]-7.5, centreofmass[1]+8.5])
    #name ='$\Delta E = {}$ eV'.format(math.ceil(energydif[j],3))
    name ='$\Delta E = {:3.3f}$ eV'.format(energydif[j])
    ax.text(0.1, -0.1, name, transform=ax.transAxes,fontsize=20)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    
    gs.update(wspace=0.00,hspace=0.00)
    plt.tight_layout()
    name = sys.argv[2]
    name =name+'_{}'.format(j)
    savefig(name,bbox_inches='tight')
    show()
    exit()