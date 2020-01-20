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

#  formation energy of oxygen
Pt7Ox=[[0, 0.0],
[2, -3.228179134242799],
[4, -6.222445977342442],
[6, -8.691274797583972],
[8, -9.618859668972185],
[10, -10.713886005601097],
[12, -10.943386683601517],
[14, -11.218929335924585],
[16, -9.40003696161061]]
print(np.shape(Pt7Ox))
print(len(Pt7Ox))
Pt7Ox = np.array(Pt7Ox)
print((Pt7Ox[:,0]))
#chemisorption of CO molecule 
Pt7OxCO=[[0, -2.6255645259353013],
[2, -5.9666031289775585],
[4, -8.460640759613906],
[6, -10.20056754566815],
[8, -10.702740792354184],
[10, -11.37352059698732],
[12, -11.58456857759127],
[14, -11.843813938516433]]
Pt7OxCO=np.array(Pt7OxCO)
#chemisorption of CO2 molecule attached to the cluster terminal
Pt7OxCO2=[[2, -5.861131701983359],
[4, -8.488919619970769],
[6, -10.543490464432352],
[8, -11.252457493847437],
[10, -12.39804226918558],
[12, -13.352578716735778],
[14, -14.377973328460513],
[16, -12.45595536962648]]
Pt7OxCO2 = np.array(Pt7OxCO2)
# chemisorption of CO3 molecule form on within cluster
Pt7OxCO3=[[2, -5.043113192502336],
[4, -6.898934154012689],
[6, -9.387412271531389],
[8, -11.208205733651717],
[10, -12.52795522340297],
[12, -14.14148595635477],
[14, -14.119295751828105],
[16, -13.057021308935617]]
Pt7OxCO3=np.array(Pt7OxCO3)
# chemisorption of CO2 dettached from cluster 
Pt7OxCO2detached=[[4, -8.091406898979184],
[6, -10.461499640039165],
[8, -11.783901333699482],
[10, -12.165900048601216],
[12, -13.622882709995864],
[14, -12.96284899020182],
[16, -10.095627075372448]]
Pt7OxCO2detached =np.array(Pt7OxCO2detached)
# color labels
colors = {}
color_lib = ['#bf40bf','#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(Pt7OxCO):
    colors[i] = color_lib[i]
# cal formation energy
######################### CO attached ###############################################################
delta_E_f_CO =np.zeros(len(Pt7OxCO))
for i in range(len(Pt7OxCO)):
    delta_E_f_CO[i] =(Pt7OxCO[i,1]-Pt7Ox[i,1])
    print(Pt7OxCO[i,0],delta_E_f_CO[i])
# plot formation E vs No. of Oxygens
fig, ax = plt.subplots(1,1,figsize=(7,7))
ax.plot(Pt7OxCO[:,0],delta_E_f_CO,marker="o",color='#0f0f0f',label='Pt$_7$O$_y$CO')
for i in range(len(Pt7OxCO)):
    ax.plot(Pt7OxCO[i,0],delta_E_f_CO[i],marker="o",color=colors[i])
######################### Co2 attached ##########################################################
colors = {}
color_lib = ['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(Pt7OxCO2):
    colors[i] = color_lib[i]
delta_E_f_CO2 =np.zeros(len(Pt7OxCO2))
for i in range(len(Pt7OxCO2)):
    delta_E_f_CO2[i] =(Pt7OxCO2[i,1]-Pt7Ox[i+1,1])
    print(Pt7OxCO2[i,0],delta_E_f_CO2[i])
# plot formation E vs No. of Oxygens
ax.plot(Pt7OxCO2[:,0],delta_E_f_CO2,marker="o",color='#8EBA42',label='Pt$_7$O$_{(y-1)}$CO$_2$')
for i in range(len(Pt7OxCO2)):
    ax.plot(Pt7OxCO2[i,0],delta_E_f_CO2[i],marker="o",color=colors[i])
######################## CO3 attached ############################################################
colors = {}
color_lib = ['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(Pt7OxCO3):
    colors[i] = color_lib[i]
delta_E_f_CO3 =np.zeros(len(Pt7OxCO3))
for i in range(len(Pt7OxCO3)):
    delta_E_f_CO3[i] =(Pt7OxCO3[i,1]-Pt7Ox[i+1,1])
    print(Pt7OxCO3[i,0],delta_E_f_CO3[i])
# plot formation E vs No. of Oxygens
ax.plot(Pt7OxCO3[:,0],delta_E_f_CO3,marker="o",color='#FF0000',label='Pt$_7$O$_{(y-2)}$CO$_3$')
for i in range(len(Pt7OxCO3)):
    ax.plot(Pt7OxCO3[i,0],delta_E_f_CO3[i],marker="o",color=colors[i])
############################### CO2 detached #####################################################
colors = {}
color_lib = ['#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(Pt7OxCO2detached):
    colors[i] = color_lib[i]
delta_E_f_CO3 =np.zeros(len(Pt7OxCO2detached))
for i in range(len(Pt7OxCO2detached)):
    delta_E_f_CO3[i] =(Pt7OxCO2detached[i,1]-Pt7Ox[i+2,1])
    print(Pt7OxCO2detached[i,0],delta_E_f_CO3[i])
# plot formation E vs No. of Oxygens
ax.plot(Pt7OxCO2detached[:,0],delta_E_f_CO3,marker="o",color='#FF00FF',label='Pt$_7$O$_{(y-1)}--$CO$_2$')
for i in range(len(Pt7OxCO2detached)):
    ax.plot(Pt7OxCO2detached[i,0],delta_E_f_CO3[i],marker="o",color=colors[i])
##################################################################################################
ax.set_ylabel(r'Extra Stability from Pt_Oxides [$\Delta$E$_f$ (eV)]')
ax.set_xlabel(r'No. of O')
ax.set_xticks(np.arange(0, 20, step=1))
plt.legend()
savefig('delta_COEf_CO2Ef_CO3Ef_O2Ef_Pt7Oxides.png')
plt.show()
