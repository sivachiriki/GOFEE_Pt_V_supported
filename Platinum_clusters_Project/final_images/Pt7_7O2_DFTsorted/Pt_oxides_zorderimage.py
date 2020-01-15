from __future__ import division
import matplotlib
#matplotlib.use('Agg')  # Can also use 'tkagg' or 'webagg'
#from plot_neb_tio2 import *
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import matplotlib.patches as patches

import matplotlib.pyplot as plt
from ase.io import read, write
from ase.visualize import view
import matplotlib.patches as mpatches
from ase.data.colors import jmol_colors

from pylab import *
from ase.data import covalent_radii as aradii
from matplotlib.patches import Circle
from math import atan2,pi
import matplotlib.gridspec as gridspec


matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rc('font',**{'family':'sans-serif',
                        'sans-serif':['Helvetica'],
                        'size':14})
matplotlib.rc('text',usetex=True)

matplotlib.rcParams['text.latex.unicode']=True
#matplotlib.rcParams['text.latex.preamble']=['\usepackage{bm}']
#matplotlib.rcParams['text.latex.preamble']=['\usepackage{xfrac}']
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'

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



def plot_conf(ax, atoms, rot=False):
    colors = np.array([jmol_colors[atom.number] for atom in atoms])
    for i, atom in enumerate(atoms):
        if (atom.number ==78):
           colors[i] =[0.1, 0.6, 0.6]
        if (atom.number ==6):
           colors[i] =[0.1, 0.2, 0.9]
        if (atom.number ==8 and i >= 135 and i <=156 ):
           colors[i] =[102/255, 0/255, 0/255]
        if (atom.number ==8 and i >= 291 and i <=312 ):
           colors[i] =[102/255, 0/255, 0/255]
        if (atom.number ==8 and i >= 449 and i <=468 ):
           colors[i] =[102/255, 0/255, 0/255]
        if (atom.number ==8 and i >= 603 and i <=624 ):
           colors[i] =[102/255, 0/255, 0/255]

    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Al' or a.symbol == 'O':
            if a.position[2] < 9.7:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)



data=read(sys.argv[1]+'@:')

for j in range(len(data)):
    image = data[j] #* (2,2,1) 
    for i,a in enumerate(image):
        if i ==48 and image.positions[i,1]>8.0 :
            image.positions[i,1] =image.positions[i,1]-12.429
            image.positions[i,0] =image.positions[i,0]+7.176
        if i ==93 and image.positions[i,1]>8.0 :
            image.positions[i,1] =image.positions[i,1]-12.429
            image.positions[i,0] =image.positions[i,0]+7.176
        if i ==3  and image.positions[i,1]>8.0 :
            image.positions[i,1] =image.positions[i,1]-12.429
            image.positions[i,0] =image.positions[i,0]+7.176
        #if i ==143:
        #    image.positions[i,1] =image.positions[i,1]+12.429
        #    image.positions[i,0] =image.positions[i,0]-7.176
        #if i ==142:
        #    image.positions[i,1] =image.positions[i,1]+12.429
        #    image.positions[i,0] =image.positions[i,0]-7.176

        if a.position[0] >7.100 and i<=151:
           image.positions[i,0] =image.positions[i,0]-14.352
    
    #write('newimage.traj',image)
    
    plt.figure(figsize=(6.0,7.0))
    
    gs = gridspec.GridSpec(2, 1,
                           height_ratios=[6.86,10.32])
    
    cell = image.get_cell()
    
    # 0 0
    ax = plt.subplot(gs[0, 0])
    image = image * (2,1,1)
    img = image.copy()
    plot_conf(ax, img)
    
    ax.set_xlim([-6.2, 10.5])
    ax.set_ylim([10.7, 20.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    
    # 0 1                                                                                                                          
    ax = plt.subplot(gs[1, 0])
    image = image * (1,2,1)
    write('newimage.traj',image)
    
    cell = image.get_cell()
    img = image.copy()
    plot_conf(ax, img, rot=True)
    
    ax.set_xlim([-6.2, 10.50])
    ax.set_ylim([0.0, 14.0])
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set(aspect=1)
    
    gs.update(wspace=0.00,hspace=0.00)
    plt.tight_layout()
    name = sys.argv[2]
    name =name+'_{}'.format(j)
    savefig(name,bbox_inches='tight')
    #show()
    #exit()
