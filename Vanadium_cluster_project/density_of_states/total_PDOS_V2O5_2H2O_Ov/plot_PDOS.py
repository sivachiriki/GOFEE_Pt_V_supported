import numpy as np
from matplotlib.pyplot import *
import glob
from gpaw import *

params={'legend.fontsize':15,
        'xtick.labelsize':16,
        'ytick.labelsize':16,
        'axes.labelsize':16,
        'axes.titlesize':16,
        'lines.linewidth':1.5}
rcParams.update(params)

dirs = ['V2O5_2H2O_GM_Ov','V2O5_2H2O_LM1_Ov']

for m, d in enumerate(dirs):
    print(d)
    slab, calc = restart('../'+d+'/top.gpw',txt=None)
    e, dos = calc.get_dos(spin=0, npts=2001, width=0.2)
    #plot(e-e_f,dos,'k-')
    e_f = calc.get_fermi_level()
    print('DOS done')

    plot(e-e_f,dos,label=d)

axis([-10, 5,0,120])
xlabel('Energy (eV)')
ylabel('DOS')
legend(loc=2)
savefig('DOS.pdf')
show()
