from gpaw import *
from matplotlib.pyplot import *
from gpaw.utilities.dos import print_projectors

#params={'legend.fontsize':13,
#        'xtick.labelsize':13,
#        'ytick.labelsize':13,
#        'axes.labelsize':13,
#        'axes.titlesize':14,
#        'lines.linewidth':1.5,
#        'axes.color_cycle':['#f781bf','#4daf4a','#a65628','#ff7f00','#e41a1c','#ffff33','#999999','#377eb8','#984ea3']}

#rcParams.update(params)

# Density of States
dirs = ['V2O5_2H2O_GM_Ov','V2O5_2H2O_LM1_Ov']
naming =['GM','LM1']
for m, d in enumerate(dirs):
    slab, calc = restart('../'+d+'/top.gpw')
    e, dos = calc.get_dos(spin=0, npts=2001, width=0.2)
    e_f = calc.get_fermi_level()
  #  plot(e-e_f, dos, lw=2, label='DOS_'+d)
    
  #  xlim([-10,5])
  #  ylim([0,120])
    
    #savefig('DOS.eps')
    
    #Ti = [i for i,a in enumerate(slab) if a.number == 22 and a.position[2] > 9.5]
    #V = [i for i,a in enumerate(slab) if a.number == 23 and a.position[2] > 9.5]
    #O = [i for i,a in enumerate(slab) if a.number == 8 and a.position[2] > 9.5]
    
    Ti = [i for i,a in enumerate(slab) if a.number == 22]
    V = [i for i,a in enumerate(slab) if a.number == 23]
    O = [i for i,a in enumerate(slab) if a.number == 8]
    
    #combs = [[Sn,{'5s':['s'],'5p':['p'],'4d':['d']}],[O,{'2s':['s'],'2p':['p']}]]
    combs = [[Ti,{'total':['d']}],[V,{'total':['d']}],[O,{'total':['p']}]]
    #combs = [[Sn,{'4d':['d'],'5px':[3],'5py':[1],'5pz':[2],'5s':[0]}],[O,{'2px':[3],'2py':[1],'2pz':[2],'2s':[0]}]]
    
    pdos_t = 0.
    for comb in combs:
        print(slab[comb[0][0]].symbol) #type of atoms
        for key in sorted(comb[1]):
            print(key) # summation over different orbitals for one type of atom
            #pdos = 0.
            for orb in comb[1][key]:
                print('    {}'.format(orb)) # different type of orbitals plotting 
                pdos = 0.
                for i in comb[0]: # summation over all the atoms position above 9.5 in z-axis.
                    e, pdosi = calc.get_orbital_ldos(a=i, angular=orb, npts=2001, width=0.5)
                    pdos += pdosi
                plot(e-e_f, pdos, label='{} {} ({})'.format(naming[m],slab[comb[0][0]].symbol,orb))
    
axis([-10, 5, None, 120])
xlabel('Energy (eV)')
ylabel('PDOS')

legend(loc=2)

savefig('PDOS.eps')

show()
