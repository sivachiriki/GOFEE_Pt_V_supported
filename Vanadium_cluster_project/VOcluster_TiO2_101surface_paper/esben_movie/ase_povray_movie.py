# -*- coding: utf-8 -*-
from ase.io import read, write
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.neb import NEB
from ase_povray import povray_parameter
import os
import numpy as np
from ase.visualize import view


def make_pic(atoms, xll, yll, xul, yul, jvej):
    if jvej < 150 or jvej > 170:
        return
    atoms = atoms * (2,2,1)
    del atoms[range(312,320)]

    PovRay = povray_parameter(atoms,
                              atoms_radii={'Ti':1.3,
                                           'O':0.74,
                                           'H':0.3,
                                           'N':0.8})
    PovRay.set_textures({'Ti':'pale',
                         'O':'pale',
                         'H':'pale',
                         'N':'pale'})

    PovRay.set_colors({'Ti':(0.749, 0.761, 0.78),
                       'O':(214./255, 81./255, 83./255),
                       'N':(0.188, 0.314, 0.973),
                       'H':(1.,1.,1.)})

    PovRay.kwargs['rotation']='0x,0y,0z'
    PovRay.kwargs['bbox']=(xll, yll, xul, yul)
    PovRay.kwargs['canvas_width']=800
    PovRay.kwargs['transparent']=False
    write('figure'+str(jvej).zfill(4)+'.pov',atoms,**PovRay.kwargs)


# load in rough path
short = read('../../left_full.traj',index=':')
trajs = []

# increase resolution by X times as many images using the idpp interpolation method
n_between = 12


for i in range(len(short) - 1):
    mstart = short[i]
    mfinal = short[i + 1]

    tot = [mstart]
    for i in range(n_between):
        tot.append(mstart.copy())
    tot.append(mfinal)

    n = NEB(tot)
    n.interpolate(method='idpp', mic=True)
    trajs += tot[:-1]

trajs.append(short[-1])


# Bounding box for images
xll = 8.489
yll = 0.089
xul = 18.403
yul = 12.403

for i, atoms in enumerate(trajs):
    make_pic(atoms, xll, yll, xul, yul, i)

os.system('mencoder mf://*.png -mf w=800:fps=20:type=png -ovc lavc -oac copy -lavcopts vcodec=mpeg4:vbitrate=5000 -o example.avi')

#os.system('rm fig*.ini fig*.pov fig*.png')
