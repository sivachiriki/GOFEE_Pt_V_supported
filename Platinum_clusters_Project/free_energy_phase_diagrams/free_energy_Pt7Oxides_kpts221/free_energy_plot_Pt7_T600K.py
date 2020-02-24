from ase.vibrations import Infrared
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.io import read,write
import matplotlib.pyplot as plt
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from ase.calculators.vasp import Vasp
from matplotlib.pyplot import *
import numpy as np
import os

#constants
kB=8.61733*10**-5 # eV/K
T = 600 # K
p0 =1.01325 * 1e3 # bar
#all the GO structures
O2 = read('O2_molecule_DFTrelaxed_kpts221.traj')
Pt7 = read('Pt7_Al2O3_DFTrelaxed_GMplanar.traj')
Pt7_O2 = read('Pt7O2_Al2O3_GMDFTrelaxed.traj')
Pt7_2O2 = read('Pt7O4_Al2O3_GMDFTrelaxed.traj')
Pt7_3O2 = read('Pt7O6_Al2O3_GMDFTrelaxed.traj')
Pt7_4O2 = read('Pt7O8_Al2O3_GMDFTrelaxed.traj')
Pt7_5O2 = read('Pt7O10_Al2O3_GMDFTrelaxed.traj')
Pt7_6O2 = read('Pt7O12_Al2O3_GMDFTrelaxed.traj')
Pt7_7O2 = read('Pt7O14_Al2O3_GMDFTrelaxed.traj')
#Pt7_8O2 = read('Pt7O16_Al2O3_GMDFTrelaxed.traj')

ga_stru =[Pt7_O2,Pt7_2O2,Pt7_3O2,Pt7_4O2,Pt7_5O2,Pt7_6O2,Pt7_7O2]

N_O = [2,4,6,8,10,12,14]

# energies of structures
E_O2 = O2.get_potential_energy()
E_Pt7 = Pt7.get_potential_energy()
E_Pt7O2 = Pt7_O2.get_potential_energy()
E_Pt72O2 = Pt7_2O2.get_potential_energy()
E_Pt73O2 = Pt7_3O2.get_potential_energy()
E_Pt74O2 = Pt7_4O2.get_potential_energy()
E_Pt75O2 = Pt7_5O2.get_potential_energy()
E_Pt76O2 = Pt7_6O2.get_potential_energy()
E_Pt77O2 = Pt7_7O2.get_potential_energy()
#E_Pt78O2 = Pt7_8O2.get_potential_energy()
###########################################
Gs_slab ={}
colors = {}
#color_lib = ['#377eb8','#4daf4a','#984ea3','#a65628','#ffff33']
color_lib = ['#377eb8','#4daf4a','#984ea3','#a65628','#999999','#fdbf6f', '#ff7f00','#ffff33']
print(color_lib)
resolution = 300
alpha_back =0.70
alpha_front =1.0

for i,atoms in enumerate(ga_stru):
    e = atoms.get_potential_energy()
    Gs_slab[i] = e
    colors[i] = color_lib[i]
# vibrational frequency calculations
calc = GPAW(mode=PW(500),xc='PBE',
            spinpol=True,
            kpts=(2,2,1),
            symmetry={'point_group': False},
            basis='dzp')
O2.set_calculator(calc)
O2vib = Infrared(O2)
#print(O2vib)
#O2vib.run()
#O2vib.summary()
e_O2vib = O2vib.get_energies()
#print(e_O2vib)

thermo = IdealGasThermo(vib_energies=e_O2vib,
                        atoms=O2,
                        geometry='linear',
                        symmetrynumber=2, spin=1)

#print(e_O2vib)
def get_mu(T,p):
    G = thermo.get_gibbs_energy(temperature=T, pressure=p0, verbose=False)
    #print(p0,T)
    mu0 = 0.5*G
    mu = 0.5*E_O2 + mu0 + 0.5*kB*T*np.log(p / (1.01325 * 1e3))
   # print(T,p,np.log(p / (1.01325 * 1e3)),mu)
    return mu


def get_required_deltaG_atmup(i,mu):
    G = delta_E[i]-N_O[i]*mu
    #print(i,G,delta_E[i])
    return G


def get_lowest_e(mu):
    lowest_e = 1000000000000000.000
    lowest_N = None
    for i in range(0,len(N_O)):
        #print(get_required_deltaG_atmup(i,mu))
        if get_required_deltaG_atmup(i,mu) <= lowest_e:
            lowest_e = get_required_deltaG_atmup(i,mu)
            lowest_N = i
    #print(lowest_N,lowest_e)
    return [lowest_N,lowest_e]

def get_lowest_e_pt(mu):
    lowest_e = 100000000000000000.000
    lowest_N = None
    for i in range(0,len(N_O)):
        if get_required_deltaG_atmup(i,mu) <= lowest_e:
            lowest_e = get_required_deltaG_atmup(i,mu)
            N =N_O[i]
            lowest_N = N
    return [lowest_N,lowest_e]


def get_highest_e(mu):
    highest_e = -100000000.00
    highest_N = None
    for i in range(0,len(N_O)):
        if get_required_deltaG_atmup(i,mu) >= highest_e :
            highest_e = get_required_deltaG_atmup(i,mu)
            lowest_N = i
    return [highest_N,highest_e]

def get_highest_e_pt(mu):
    highest_e = -100000000.00
    highest_N = None
    for i in range(0,len(N_O)):
        if get_required_deltaG_atmup(i,mu) >= highest_e :
            highest_e = get_required_deltaG_atmup(i,mu)
            N =N_O[i]
            lowest_N = N
    return [highest_N,highest_e]
 

ps = np.power(10,[-24.,-22.,-20.,-18.,-16.,-14.,-12.,-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.,12.,14.,16., 18.,20.,22.,24.,26.,28.,30.])
#print(ps)
#exit()
mu_p = [get_mu(T,p) for p in ps]
print(mu_p)
delta_E =np.zeros(len(ga_stru))
for i,atoms in enumerate(ga_stru):
    delta_E[i] = (ga_stru[i].get_potential_energy() - E_Pt7) 
   # print(delta_E[i])

#fig, ax = plt.subplots()
fig = figure(figsize=(10,5))
ax1 = fig.add_subplot(121)

delta_G =np.zeros([len(ga_stru),len(mu_p)])
for i in range(len(ga_stru)):
    for j in range(len(mu_p)):
        delta_G[i,j] = delta_E[i]-N_O[i]*mu_p[j]
        print(delta_G[i,j])
    ax1.plot(mu_p, delta_G[i,:],'-',color=colors[i],label='Pt$_{7}$O$_{'+str(N_O[i])+'}$',zorder=2)
   # print(delta_G[i,:])

#xlims =[min(mu_p),max(mu_p)]
xlims =[-7.0,-4.0]
print(xlims)
ax1.set_xlim(xlims)
# Lowest line and fill
#ylims = [get_lowest_e(mu_p[0])[1]*1.1,get_highest_e(mu_p[-1])[1]*1.1]
ylims =[-60.0,20.0]
ax1.set_ylim(ylims)
#print(ylims)
xticks([-7.0,-6.5,-6.0,-5.5,-5.0,-4.5,-4.0],[r'-7.0','-6.5','-6.0','-5.5','-5.0','-4.5','-4.0'])
yticks([-60.0,-50.0,-40.0,-30.0,-20.0,-10.0,0.0,10.0,20.],['-60.0','-50.0','-40.0','-30.0','-20.0','-10.0','0.0','10.0','20.0'])

mu_first = mu_p[0]
#print(mu_first)
#print(np.linspace(mu_p[0],mu_p[-1]+1e-10, 200))
for mu_last in np.linspace(mu_p[0],mu_p[-1]+1e-10, 200):
    off = 0.
    if get_lowest_e(mu_first)[0] != get_lowest_e(mu_last)[0] or mu_last > mu_p[-1]:
        if mu_first > -1.:
            off = 0.01
        #print(mu_first)
        ax1.plot([mu_first+off, mu_last], [get_lowest_e(mu_first)[1], get_lowest_e(mu_last)[1]],
             '-', lw=4, color=colors[get_lowest_e(mu_first)[0]], alpha=alpha_front, zorder=3)
        ax1.fill_between([mu_first, mu_last],[ylims[-1], ylims[-1]], [ylims[0],ylims[0]], lw=0.,
                         color=colors[get_lowest_e(mu_first)[0]], alpha=alpha_back, zorder=1)
        mu_first = mu_last

ax1.legend(loc='upper right', frameon=False)
ax1.set_xlabel(r'$\mu_\mathrm{O}(T,p)$ (eV)')
ax1.set_ylabel(r'$\Delta G$ (eV)')
ax1.set_title('Free energy diagram at 600 K Temp.')
# p-T diagram
ax2 = fig.add_subplot(122)
T_range = np.linspace(100, 1200, resolution)
p_range = np.power(10,np.linspace(-24, 30, resolution))
z = np.zeros((len(T_range), len(p_range)))
mu_z = np.zeros((len(T_range), len(p_range)))

levels = []
levels1 = []
for i,T in enumerate(T_range):
    for j,p in enumerate(p_range):
        mu_z[i][j] = get_mu(T,p)
        #print(get_lowest_e(mu_z[i][j])[0])
        z[i][j] = get_lowest_e_pt(mu_z[i][j])[0]
        if get_lowest_e(mu_z[i][j])[0] not in levels:
            levels.append(get_lowest_e(mu_z[i][j])[0])
            levels1.append(get_lowest_e_pt(mu_z[i][j])[0])
    #if (T >100):
    #   exit() 


mu_z = np.transpose(mu_z)
z = np.transpose(z)
#print(z[:][:])

color_levels = []
for level in levels:
    color_levels.append(colors[level])

print(color_levels)
levels = [l for l in levels1]
print(levels)

levels.append(0)
print(levels)
# levels index
levels = np.sort(levels)
print(levels)
changedvalues =[4,3,2,1,0]
color_levels1=[]
for i,a in enumerate(changedvalues):
    color_levels1.append(color_levels[a])
print(color_levels1)
ax2.contourf(T_range, p_range, z,levels=levels,colors=color_levels1, alpha=alpha_back, antialiased=True)
ax2.set_yscale('log')
ax2.set_ylim([p_range[0],p_range[-1]])
ax2.set_xlim([T_range[0],T_range[-1]])
yticks(np.power(10,np.linspace(-24,30,22)),range(-24,30,2))
xticks(range(100,1201,200),range(100,1201,200))
ax2.set_title('Phase diagram for Pt$_7$O$_y$')
ax2.set_ylabel(r'log$_{10}$($\frac{p}{p_0}$)')
ax2.set_xlabel(r'$T$ (K)')
#ax2.text(0.1, 0.70,'Pt$_7$O$_{16}$',color=color_levels1[5],alpha=1.0,fontsize=16, transform=ax2.transAxes)
ax2.text(0.1, 0.60,'Pt$_7$O$_{14}$',color=color_levels1[4],alpha=1.0,fontsize=16, transform=ax2.transAxes)
#ax2.text(0.2, 0.5,'Pt$_7$O$_{12}$',fontsize=16,color=color_levels1[4],alpha=1.0,transform=ax2.transAxes)
ax2.text(0.5, 0.6,'Pt$_7$O$_{10}$',fontsize=16, color=color_levels1[3],alpha=1.0,transform=ax2.transAxes)
#ax2.text(0.3, 0.4,'Pt$_7$O$_{8}$',fontsize=16, color=color_levels1[3],alpha=1.0,transform=ax2.transAxes)
ax2.text(0.4, 0.3,'Pt$_7$O$_{6}$',fontsize=16, color=color_levels1[2],alpha=1.0,transform=ax2.transAxes)
ax2.text(0.7, 0.45,'Pt$_7$O$_{4}$',fontsize=16,color=color_levels1[1],alpha=1.0,transform=ax2.transAxes)
ax2.text(0.7, 0.2,'Pt$_7$O$_{2}$',fontsize=16,color=color_levels1[0],alpha=1.0,transform=ax2.transAxes)

savefig('free_energy_diagram_Pt7Oxides_600K.png')
plt.show()
