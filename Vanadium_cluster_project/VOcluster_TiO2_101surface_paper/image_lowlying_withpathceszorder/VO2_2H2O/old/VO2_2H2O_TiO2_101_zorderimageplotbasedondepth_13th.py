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
        if atoms[ia].symbol == 'Ti':
            arad = aradii[atoms[ia].number] #* 0.9 * 0.5
        else:
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
    alp = [None] * colors.shape[0]
    for i,a in enumerate(atoms):
        if a.symbol == 'Ti' or a.symbol == 'O':
            if a.position[2] < 11.95:
                alp[i] = 0.3

    if rot:
        atoms.rotate('x',pi/2)
    plot_atoms(ax, atoms, [0,2,1], colors, alp, z=-1)



data=read('VO2_2H2O_TiO2_101_DFTlowlyingsorted.traj@:')
image = data[12] * (2,2,1)
#image = read('anataseTi24O48_101surface_optPBEesben.traj') * (2,2,1)
write('221image.traj',image)
# Make array of indices for atoms that should be repeated in x and y directions                                       
plt.figure(figsize=(6.0,10.0))

gs = gridspec.GridSpec(2, 1,
                       height_ratios=[7.0,9.5])

cell = image.get_cell()

#image.translate(-cell[:,0]/3.)
#image.translate(-cell[:,1]/3.)
#view(image)
# 0 0
ax = plt.subplot(gs[0, 0])
img = image.copy()
plot_conf(ax, img)

print(cell[0,0])
print(cell[1,1])
ax.set_xlim([1.5, cell[0,0]/1.8])
ax.set_ylim([13.5, 22.0])
ax.set_yticks([])
ax.set_xticks([])
ax.set(aspect=1)

# 0 1                                                                                                                          
ax = plt.subplot(gs[1, 0])
img = image.copy()
plot_conf(ax, img, rot=True)

ax.set_xlim([1.5, cell[0,0]/1.8])
ax.set_ylim([0, 11.40])
ax.set_yticks([])
ax.set_xticks([])
ax.set(aspect=1)

gs.update(wspace=0.00,hspace=0.00)
plt.tight_layout()
savefig('VO2_2H2O_TiO2_101_DFTopt_13th.pdf',bbox_inches='tight')
show()
