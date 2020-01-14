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
################################# E_form of O2 ####################################################
# read stru of GM Pt7 oxides
O2 = read('O2_kpts221.traj')
Pt7O0 =  read('Pt7_Al2O3_DFTrelaxed_GMplanar.traj')
Pt7O2 =  read('Pt7O2_Al2O3_GMDFTrelaxed.traj')
Pt7O4 =  read('Pt7O4_Al2O3_GMDFTrelaxed.traj')
Pt7O6 =  read('Pt7O6_Al2O3_GMDFTrelaxed.traj')
#Pt7O7 =  read('Pt7O7_Al2O3_GMDFTrelaxed.traj')
Pt7O8 =  read('Pt7O8_Al2O3_GMDFTrelaxed.traj')
Pt7O10 = read('Pt7O10_Al2O3_GMDFTrelaxed.traj')
Pt7O12 = read('Pt7O12_Al2O3_GMDFTrelaxed.traj')
Pt7O14 = read('Pt7O14_Al2O3_GMDFTrelaxed.traj')
Pt7O16 = read('Pt7O16_Al2O3_GMDFTrelaxed.traj')
Pt7O18 = read('Pt7O18_GM.traj')
# energies of GM stru of Pt7 oxides
E_O2 = O2.get_potential_energy()
E_Pt7O0 = Pt7O0.get_potential_energy()
E_Pt7O2 = Pt7O2.get_potential_energy()
E_Pt7O4 = Pt7O4.get_potential_energy()
E_Pt7O6 = Pt7O6.get_potential_energy()
#E_Pt7O7 = Pt7O7.get_potential_energy()
E_Pt7O8 = Pt7O8.get_potential_energy()
E_Pt7O10 = Pt7O10.get_potential_energy()
E_Pt7O12 = Pt7O12.get_potential_energy()
E_Pt7O14 = Pt7O14.get_potential_energy()
E_Pt7O16 = Pt7O16.get_potential_energy()
E_Pt7O18 = Pt7O18.get_potential_energy()

gm_tru =[Pt7O0,Pt7O2,Pt7O4,Pt7O6,Pt7O8,Pt7O10,Pt7O12,Pt7O14,Pt7O16]
N_O = [0,2,4,6,8,10,12,14,16]
#N_O2 =[0,1,2,3,4,5,6,7,8,9]
#print(N_O2)
# color labels
colors = {}
color_lib = ['#bf40bf','#377eb8','#4daf4a','#984ea3','#a65628','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru):
    colors[i] = color_lib[i]
# cal formation energy
Oxygen_correct = -0.5
E_f =np.zeros(len(gm_tru))
E_f_correct =np.zeros(len(gm_tru))
for i,atoms in enumerate(gm_tru):
    E_f_correct[i] = (gm_tru[i].get_potential_energy() - gm_tru[0].get_potential_energy() - (N_O[i])*(0.50000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_f_correct[i])
E_f_correct[0] =0.0
# plot formation E vs No. of Oxygens
fig, ax = plt.subplots(1,1,figsize=(7,7))
ax.plot(N_O,E_f_correct,marker="o",color='#1E90FF',label='E$_{formation}$ of O$_2$ on Pt$_7$/Al$_2$O$_3$(0001)')
for i, nu in enumerate(N_O):
    ax.plot(N_O[i],E_f_correct[i],marker="o",color=colors[i])

################################# E_ads of CO ####################################################
# read stru of GM Pt7 oxides
CO = read('CO_molecule_DFTrelaxed.traj')
# energies of GM stru of Pt7 oxides
E_CO = CO.get_potential_energy()
######### Pt7-Oxides_CO ###########
Pt7O0CO = read('Pt7_CO_Al2O3_Eads_planarDFTrelaxed.traj')
Pt7O2CO = read('Pt7O2_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O4CO = read('Pt7O4_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O8CO = read('Pt7O8_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O10CO = read('Pt7O10_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O12CO = read('Pt7O12_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O14CO = read('Pt7O14_CO_Al2O3_Eads_DFTrelaxed.traj')
Pt7O16CO = read('Pt7O16_CO_Al2O3_Eads_DFTrelaxed.traj')

E_Pt7O0CO = Pt7O0CO.get_potential_energy()
E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()
E_Pt7O16CO = Pt7O14CO.get_potential_energy()

gm_tru_CO =[Pt7O0CO,Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO,Pt7O16CO]
print(gm_tru_CO[-1].get_potential_energy())
N_O = [0,2,4,6,8,10,12,14,16]
# color labels
colors = {}
color_lib = ['#bf40bf','#377eb8','#4daf4a','#984ea3','#a65628','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
ax.plot(N_O,E_ads,marker="o",color='#0f0f0f',label='E$_{ads}$ of CO on Pt$_7$O$_y$/Al$_2$O$_3$(0001)')
for i, nu in enumerate(N_O):
    ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])

####################   CO chemisorption energy for Pt7 #############################################
Pt7O0CO = read('Pt7_CO_Al2O3_Eads_planarDFTrelaxed.traj')
Pt7O2CO = read('Pt7O2_CO_Al2O3_GMDFTrelaxed.traj')
Pt7O4CO = read('Pt7O4CO_GM_Al2O3_0001surface_DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6_CO_Al2O3_GMDFTrelaxed.traj')
#Pt7O7CO = read('Pt7O7CO_Al2O3_GMDFTrelaxed.traj')
Pt7O8CO = read('Pt7O8CO_GM_Al2O3_KRRfund9l_DFTrelaxed.traj')
Pt7O10CO = read('Pt7O10CO_GM_Al2O3_0001surface_DFTrelaxed.traj')
Pt7O12CO = read('Pt7O12_CO_Al2O3_GMDFTrelaxed.traj')
Pt7O14CO = read('Pt7O14CO_Al2O3_GMDFTrelaxed.traj')
Pt7O16CO = read('Pt7O16CO_Al2O3_GMDFTrelaxedsorted.traj')

E_Pt7O0CO = Pt7O0CO.get_potential_energy()
E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
#E_Pt7O7CO = Pt7O7CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()
E_Pt7O16CO = Pt7O16CO.get_potential_energy()

gm_tru_CO =[Pt7O0CO,Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO,Pt7O16CO]
#print(gm_tru_CO[-1].get_potential_energy())
N_O = [0,2,4,6,8,10,12,14,16]
# color labels
colors = {}
#color_lib = ['#bf40bf','#4daf4a','#a65628','#999999', '#ff7f00','#eeefff']
#color_lib = ['#bf40bf','#377eb8','#4daf4a','#984ea3','#a65628','#999999','#fdbf6f', '#ff7f00','#ff7f00','#eeefff','#ffff33']
color_lib = ['#bf40bf','#377eb8','#4daf4a','#984ea3','#a65628','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
ax.plot(N_O,E_ads,marker="o",color='#8EBA42',label='E$_{chem}$ of CO on Pt$_7$O$_y$/Al$_2$O$_3$(0001)')
for i, nu in enumerate(N_O):
    ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])
#########################################################################################
ax.set_ylabel(r'E (eV)')
ax.set_xlabel(r'No. of O')
ax.set_xticks(np.arange(0, 20, step=1))
plt.legend()
savefig('COEads_CoEchem_O2Ef_Pt7Oxides.png')
plt.show()
