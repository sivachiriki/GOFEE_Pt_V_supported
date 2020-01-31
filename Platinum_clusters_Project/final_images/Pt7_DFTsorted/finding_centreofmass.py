from __future__ import division
import matplotlib
#matplotlib.use('Agg')  # Can also use 'tkagg' or 'webagg'
#from plot_neb_tio2 import *
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
import matplotlib.patches as patches
from ase import Atoms
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

a = read('Pt7_Al2O3_Planar_GMDFTrelaxed.traj')
a = a * (3,3,1)
atoms=a
#view(a)
#a.center()
#del atoms[atoms.numbers <=702]
#del atoms[atoms.numbers >=710]
del atoms[[atom.index for atom in atoms if atom.index <=702 or atom.index >=710]]
#del atoms[[atom.index for atom in atoms if atom.index >=710]]
#del atoms[[atom.index for atom in atoms if atom.symbol == 'O']]
#del atoms[[atom.index for atom in atoms if atom.symbol == 'Al']]
centreofmass = a.get_center_of_mass()
print(centreofmass)
a=read('Pt7_Al2O3_Planar_GMDFTrelaxed.traj')
a = a * (3,3,1)
atoms=a
#positions = atoms.get_positions()
#print(positions[0,0])
#for i, atom in enumerate(atoms):
del atoms[atoms.positions[:,0] >=centreofmass[0]+7.0]
del atoms[atoms.positions[:,0] <= centreofmass[0]-7.0]
del atoms[atoms.positions[:,1] >= centreofmass[1]+9.0]
del atoms[atoms.positions[:,1] <= centreofmass[1]-9.0]
#del atoms[[atom.index for atom in atoms if positions[atom.number,0] <=702]]
view(a)
