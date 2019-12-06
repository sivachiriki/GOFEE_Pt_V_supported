from ase.infrared import InfraRed
from ase.thermochemistry import IdealGasThermo
from ase.io import read,write
from matplotlib.pyplot import *
import numpy as np
import os

#rcParams['mathtext.default'] = 'default'

fontProperties = {'family':'sans-serif','sans-serif':['Helvetica']}
rc('font',**fontProperties)
rcParams['text.latex.unicode'] = True
rc('text', usetex=True)

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Helvetica'
matplotlib.rcParams['mathtext.it'] = 'Helvetica:italic'
matplotlib.rcParams['mathtext.bf'] = 'Helvetica:bold'
matplotlib.rcParams['mathtext.sf'] = 'Helvetica'

params={'lines.linewidth':1.5,
        'legend.fontsize':14,
        'xtick.labelsize':16,
        'ytick.labelsize':16,
        'axes.labelsize':16,
        'axes.linewidth':1,
        'legend.linewidth':1}

rcParams.update(params)

#rcParams['mathtext.default'] = 'regular'
e_H2O = read('/home/matjoerg/phd/ga/SnO2/ref_structures/H2O/optimisation.traj').get_potential_energy()
e_H2 = read('/home/matjoerg/phd/ga/SnO2/ref_structures/H2/optimisation.traj').get_potential_energy()

kJ_to_eV = 6.242e21
Avo = 6.022e23
H_H2O = -241.826 * kJ_to_eV / Avo
#e_O = e_H2O - e_H2 - H_H2O
e_O = 0.5 * -9.158
e_Sn = 0.25 * read('/home/matjoerg/phd/ga/SnO2/ref_structures/bulk/Sn/tetra/a1.02c0.985/optimisation.traj').get_potential_energy()
e_SnO2 = 0.5 * read('/home/matjoerg/phd/ga/SnO2/ref_structures/bulk/k8x8x6/a1.0c0.99/optimisation.traj').get_potential_energy()
e_SnO2_formation = e_SnO2 - (e_Sn + e_O*2) # Lower boundary
A = 89.69
kB=8.61733*10**-5
T = 1100

resolution = 300
alpha_back = 0.35
alpha_front = 1.0

#gaO4 = read('/home/matjoerg/phd/ga/SnO2/red/gaO4/optimisation.traj')
gaO6 = read('/home/matjoerg/phd/ga/SnO2/red/gaO6/bbottom/optimisation.traj')
batzill = read('/home/matjoerg/phd/ga/SnO2/red/batzill/bbottom/optimisation.traj')
gaO8 = read('/home/matjoerg/phd/ga/SnO2/red/gaO8c/optimisation.traj')
gaO10 = read('/home/matjoerg/phd/ga/SnO2/red/gaO10/optimisation.traj')
gaO12 = read('/home/matjoerg/phd/ga/SnO2/red/gaO12c/optimisation.traj')
atrei = read('/home/matjoerg/phd/ga/SnO2/red/atrei/reduced/optimisation.traj')
Obri = read('/home/matjoerg/phd/ga/SnO2/ref_structures/SnO2_correct/bulkbottom/optimisation.traj')
gaOsto = read('/home/matjoerg/phd/ga/SnO2/ref_structures/SnO2_correct/not_reduced/bulkbottom/optimisation.traj')

ga = [gaO6, gaO8, batzill, gaO10, gaO12, atrei, Obri, gaOsto]

G_slab = gaOsto.get_potential_energy()
N_Sn = len([a for a in gaOsto if a.number == 50])
N_O = len([a for a in gaOsto if a.number == 8])
e_surfb = 1/(2*A) * (G_slab - 0.5*N_O*e_SnO2 - (N_Sn - 0.5*N_O) * (e_Sn))

#color_lib = ['#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#8c2d04','b']
#color_lib = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494']
#color_lib = ['#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#91bfdb','#4575b4']
color_lib = ['#377eb8','#4daf4a','#984ea3','#a65628','#ffff33','#ff7f00','#e41a1c','#999999']
#color_lib = ['#ffffb2','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#b10026']
#color_lib = ['#8c510a','#bf812d','#dfc27d','#c7eae5','#80cdc1','#35978f','#01665e']


Gs_slab = {}
colors = {}
for i,atoms in enumerate(ga):
    N_Sn = len([a for a in atoms if a.number == 50])
    N_O = len([a for a in atoms if a.number == 8])
    e = atoms.get_potential_energy()
    key = (i,N_Sn,N_O)
    Gs_slab[key] = e
    colors[key] = color_lib[i]

cwd = os.getcwd()
os.chdir('/home/matjoerg/phd/ga/SnO2/ref_structures/O2/vib/')
O2 = read('optimisation.traj')
O2vib = InfraRed(O2)
e_O2vib = O2vib.get_energies()
thermo = IdealGasThermo(vib_energies=e_O2vib,
                        atoms=O2,
                        geometry='linear',
                        symmetrynumber=2, spin=1)
os.chdir(cwd)

#mu0s = {100: -0.08, 200: -0.17, 300: -0.27, 400: -0.38, 500: -0.50, 600: -0.61, 700: -0.73, 800: -0.85, 900: -0.98, 1000: -1.10}

def free_e(mu,i,N_Sn,N_O):
    G_slab = Gs_slab[(i,N_Sn,N_O)]

    gamma_Opoor = 1/A * (G_slab - 0.5*N_O*e_SnO2 - (N_Sn - 0.5*N_O) * (e_Sn)) - e_surfb# - e_SnO2))
    gamma_Orich = gamma_Opoor - 1/A * (N_Sn - 0.5*N_O) * e_SnO2_formation
    slope = (gamma_Orich - gamma_Opoor) / (e_O - e_SnO2_formation)

    e = gamma_Orich + slope * mu
    return e

def get_lowest_e(mu):
    lowest_e = None
    lowest_N = None
    for N in Gs_slab.keys():
        if free_e(mu,N[0],N[1],N[2]) < lowest_e or lowest_e == None:
            lowest_e = free_e(mu,N[0],N[1],N[2])
            lowest_N = N
    return [lowest_N,lowest_e]

def get_highest_e(mu):
    highest_e = None
    highest_N = None
    for N in Gs_slab.keys():
        if free_e(mu,N[0],N[1],N[2]) > highest_e or highest_e == None:
            highest_e = free_e(mu,N[0],N[1],N[2])
            highest_N = N
    return [highest_N,highest_e]

def get_p(mu,T):
    G = thermo.get_gibbs_energy(temperature=T, pressure=thermo.referencepressure, verbose=False)
    mu0 = 0.5*G
#    mu0 = mu0s[T]

    p = np.exp(2*(mu-mu0)/(kB*T)) * 1.01325 * 1e3
    return p

def get_mu(T,p):
    G = thermo.get_gibbs_energy(temperature=T, pressure=thermo.referencepressure, verbose=False)
    mu0 = 0.5*G
#    mu0 = mu0s[T]

    mu = mu0 + 0.5*kB*T*np.log(p / (1.01325 * 1e3))
    return mu

fig = figure(figsize=(14,6))
ax1 = fig.add_subplot(121)
#ax2 = ax1.twiny()

mus = np.array([0.5 * e_SnO2_formation, 0.])

for N in sorted(Gs_slab.keys()):
    label = r'Sn$_{}$O$_{}$'.format(str(N[1]-24),'{'+str(N[2]-48)+'}')
    if N[0] == 2:
        label += ' (Batzill)'
    if N[0] == 5:
        label += ' (Atrei)'
    ax1.plot(mus,free_e(mus,N[0],N[1],N[2]),'-',color=colors[N],label=label,zorder=2)

ax1.set_xlim(mus)
ylims = [get_lowest_e(mus[0])[1]*1.1,get_highest_e(mus[-1])[1]*1.1]
ax1.set_ylim(ylims)

#ax2.set_ylim(ylims)
#ax2.set_xlim(ax1.get_xlim())

ps = np.power(10,[-10.,-5.,0.,5.,10.])
mu_p = [get_mu(T,p) for p in ps]
#ax2.set_xticks(mu_p)
#ax2labels = [r'{:.0e}'.format(p) for p in [get_p(mu,T) for mu in mu_p]]
#ax2labels = ps
#ax2.set_xticklabels(ax2labels)
#ax2.set_xlabel(r'O$_2$ pressure at {} K (mbar)'.format(T))

ax1.set_xlabel(r'$\mu_\mathrm{O}(T,p)$ $-$ $\frac{1}{2}E_{\mathrm{O}_2}$ (eV)')
ax1.set_ylabel(r'Surface energy (eV/Å$^2$)')
xticks([-2.5,-2.0,-1.5,-1.0,-0.5,0.0],[r'-2.5','-2.0','-1.5','-1.0','-0.5','0.0'])
yticks([-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3],['-0.4','-0.3','-0.2','-0.1','0.0','0.1','0.2','0.3'])

# Lowest line and fill
mu_first = mus[0]
for mu_last in np.linspace(mus[0],mus[-1]+1e-10, 200):
    off = 0.
    if get_lowest_e(mu_first)[0] != get_lowest_e(mu_last)[0] or mu_last > mus[-1]:
        if mu_first > -1.:
            off = 0.01
        print mu_first
        ax1.plot([mu_first+off, mu_last], [get_lowest_e(mu_first)[1], get_lowest_e(mu_last)[1]], 
             '-', lw=4, color=colors[get_lowest_e(mu_first)[0]], alpha=alpha_front, zorder=3)
        ax1.fill_between([mu_first, mu_last],[ylims[-1], ylims[-1]], [ylims[0],ylims[0]], lw=0.,
                         color=colors[get_lowest_e(mu_first)[0]], alpha=alpha_back, zorder=1)
        mu_first = mu_last

leg = ax1.legend(loc=2)

# p-T diagram
ax3 = fig.add_subplot(122)

T_range = np.linspace(200, 1200, resolution)
p_range = np.power(10,np.linspace(-15, 10, resolution))

z = np.zeros((len(T_range), len(p_range)))
mu_z = np.zeros((len(T_range), len(p_range)))
levels = []
for i,T in enumerate(T_range):
    for j,p in enumerate(p_range):
        mu_z[i][j] = get_mu(T,p)
        z[i][j] = get_lowest_e(mu_z[i][j])[0][2]
        if get_lowest_e(mu_z[i][j])[0] not in levels:
            levels.append(get_lowest_e(mu_z[i][j])[0])

mu_z = np.transpose(mu_z)
z = np.transpose(z)

color_levels = []
for level in levels:
    color_levels.append(colors[level])
levels = [l[2] for l in levels]
levels.append(0)
levels = levels[::-1]
color_levels = color_levels[::-1]

ax3.contourf(T_range, p_range, z, levels=levels, colors=color_levels, alpha=alpha_back, antialiased=True)

ax3.set_yscale('log')
ax3.set_ylim([p_range[0],p_range[-1]])
ax3.set_xlim([T_range[0],T_range[-1]])
yticks(np.power(10,np.linspace(-13,9,12)),range(-13,11,2))
xticks(range(200,1201,200),range(200,1201,200))

ax3.set_ylabel(r'log$_{10}$($\frac{p}{p_0}$)')
ax3.set_xlabel(r'$T$ (K)')

text(650,5e-3,'Sn$_6$O$_6$',fontsize=16)
text(515,5e0,'Sn$_8$O$_{16}$',fontsize=16)
text(750,2e-5,'Sn$_6$O$_4$',fontsize=16)
text(850,1e-7,'Sn$_6$O$_2$',fontsize=16)

tight_layout()
savefig('phase_diagram_pT.pdf')

show()
