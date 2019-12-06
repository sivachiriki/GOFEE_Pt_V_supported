# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import sys
from ase.io import read
from ase.infrared import InfraRed
from ase.thermochemistry import IdealGasThermo
import os
import matplotlib as mpl
from ase import units
from matplotlib.patches import Rectangle


#energy inputs
energy=np.loadtxt(sys.argv[1])
number_zn=int(sys.argv[1][2])
en_o = 0.5* -10.161432
en_zn = -1.36
en_surf=-703.581017
en_h2o=-14.349087
en_h2=-6.729133
en_o_zpe=0.0
kB=8.61733*10**-5

#boundaries
zno_decomp=-7.91
surf_alloy=-6.79
bulk_cu2o=-6.14

#plot specific stuff
resolution=10 #increase to 300 to have round lines in the p-T-plot
alpha_back=0.5 
alpha_front=0.7


def create_color_dictionary(numbers):
    colors=['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a']        
    color_dict={}
    for i,number in enumerate(numbers):
        color_dict[int(number)] = colors[i]
    return color_dict

def en(number_oxygen,mu):
    # calculate Free energy for a given nmber of oxygen. 
    return (energy[number_oxygen][1]+0.09*number_oxygen - number_zn*(en_zn) - en_surf - (number_oxygen)*(mu))

def get_lowest(mu):
     return lowest free energy value and the number of oxygens  at given mu
    lowest = None
    lowest_number=None
    for number, value in energy:
        if en(number,mu) < lowest or lowest == None:
            lowest=en(number,mu)
            lowest_number=number
    return [lowest_number,lowest]


def get_mu(p, T):
    # calculates mu for a given pressure and temperature
    os.chdir('vibO2')
    atoms = read('relax.traj')
    vib = InfraRed(atoms)
    vib_energies = vib.get_energies()
    thermo = IdealGasThermo(vib_energies=vib_energies,
                            atoms=atoms,
                            geometry='linear',
                            symmetrynumber=2, spin=1)
    G = thermo.get_gibbs_energy(temperature=T, pressure=p, verbose=False)
    delta_mu = 0.5*G#-thermo.get_ZPE_correction())
    os.chdir('..')
    return delta_mu+en_o

def get_p(mu, T):
    # calculates p for a given mu and temperature
    os.chdir('vibO2')
    atoms = read('relax.traj')
    vib = InfraRed(atoms)
    vib_energies = vib.get_energies()
    thermo = IdealGasThermo(vib_energies=vib_energies,
                            atoms=atoms,
                            geometry='linear',
                            symmetrynumber=2, spin=1)
    G = thermo.get_gibbs_energy(temperature=T, pressure=thermo.referencepressure, verbose=False)
    mu_0 = 0.5*G+en_o
    p=thermo.referencepressure*np.exp(2*(mu-mu_0)/(T*kB))
    os.chdir('..')
    return p


colors = create_color_dictionary(energy[:,0])
mu = np.linspace(-8.5,-5,2)
ylimits=[-20,10]


ax=plt.subplot(121)
# plot free energy diagram
for number, value in energy:
    plt.plot(mu-(en_o+en_o_zpe), en(number,mu), '-', label ='Zn$_%d$O$_%d$' % (number_zn, number), color=colors[int(number)], lw=2)

#plot lowest line thick and fill space with color
start=mu[0]
for number in np.linspace(mu[0],mu[-1], 200):
    if get_lowest(number)[0] != get_lowest(start)[0]:
        plt.plot([start-(en_o+en_o_zpe), number-(en_o+en_o_zpe)], [get_lowest(start)[1], get_lowest(number)[1]],
                 '-', lw=4, color=colors[int(get_lowest(start)[0])], alpha=alpha_front)
        plt.fill_between([start-(en_o+en_o_zpe), number-(en_o+en_o_zpe)],[ylimits[-1], ylimits[-1]], [ylimits[0],ylimits[0]],
                         color=colors[int(get_lowest(start)[0])], alpha=alpha_back)
        start=number

plt.plot([start-(en_o+en_o_zpe), 0.2-(en_o+en_o_zpe)],[get_lowest(start)[1], get_lowest(0.2)[1]], '-',
         lw=4, color=colors[int(get_lowest(start)[0])], alpha=alpha_front)
plt.fill_between([start-(en_o+en_o_zpe), 0.2-(en_o+en_o_zpe)],[ylimits[-1], ylimits[-1]], [ylimits[0],ylimits[0]],
                 color=colors[int(get_lowest(start)[0])], alpha=alpha_back)


#plot boundaries between oxides etc
for temp_energy in [zno_decomp, surf_alloy, bulk_cu2o]:
    plt.axvline(temp_energy-en_o-en_o_zpe, color='w', linestyle='-', linewidth=3)

#plot experimental area
plt.axvline(-2.56, color='k', linestyle='--')
plt.axvline(-2.82, color='k', linestyle='--')
plt.fill_between([-2.56, -2.82],[ylimits[-1], ylimits[-1]], [ylimits[0],ylimits[0]], color="none", hatch="\\", edgecolor="k", linewidth=1.0)


plt.xlabel('$\mu_O$(T,p)-$\\frac{1}{2}E_{\\mathrm{O}_2}$', fontsize=16)
plt.ylabel('$\Delta$ G [eV]', fontsize=16)


#plot legend
handles, labels = ax.get_legend_handles_labels()

p = Rectangle((0, 0), 1, 1, hatch="\\\\", edgecolor="k", fill=False)
handles.append(p)
labels.append("CH$_3$OH \nSynthesis")

legend=plt.legend(handles, labels, bbox_to_anchor=(0., 1.02, 2.3, .102), loc=3, ncol=number_zn+2, mode="expand", borderaxespad=0., handletextpad=0.2, fontsize=16)
for j,legobj in enumerate(legend.legendHandles):
    if j==len(legend.legendHandles)-1:
        continue
    legobj.set_linewidth(3.0)
    legobj.set_alpha(alpha_front)

plt.ylim(ylimits)
plt.xlim([mu[0]-(en_o+en_o_zpe), mu[1]-(en_o+en_o_zpe)])

#plot additional x axis with const temp
plt.subplots_adjust(bottom=0.3)
plt.subplots_adjust(top=0.88)
fig=plt.gcf()
fig.subplots_adjust(wspace=.3)
temps=[400,700,900]
for i in range(3):
    # second x-axis
    ax2 = ax.twiny()
    temp=temps[i]
    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.15-i*0.07))
    ticklabelpad = mpl.rcParams['xtick.major.pad']
    ax2.annotate('T='+str(temp)+'K, log$_{10}$($\\frac{p}{1 \mathrm{Pa}}$)', xy=(1.03,-0.15-i*0.07), xytext=(5, -ticklabelpad), ha='left', va='top',
                 xycoords='axes fraction', textcoords='offset points', fontsize=16)
    new_tick_locations = np.linspace(mu[0]-(en_o+en_o_zpe), mu[1]-(en_o+en_o_zpe),10)
    pressures=np.zeros(len(new_tick_locations))
    for i,tick in enumerate(new_tick_locations):
        pressures[i]=get_p(tick+en_o+en_o_zpe,temp)
    # Turn on the frame for the twin axis, but then hide all 
    # but the bottom spine
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    for sp in ax2.spines.itervalues():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xticks(new_tick_locations)
    pressures=[int(round(np.log10(pressure),0)) for pressure in pressures]
    ax2.set_xticklabels(pressures)
    ax2.set_xlim([mu[0]-(en_o+en_o_zpe), mu[1]-(en_o+en_o_zpe)])


# add xaxis constant pressure
ax2 = ax.twiny()
temps=[100,300,500,700,900,1100,1300, 1500]
# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")
# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.15-3*0.07))
ticklabelpad = mpl.rcParams['xtick.major.pad']
ax2.annotate('T, p=1O$^{-6}$ Pa', xy=(1.03,-0.15-3*0.07), xytext=(5, -ticklabelpad), ha='left', va='top',
                 xycoords='axes fraction', textcoords='offset points', fontsize=16)

new_tick_locations = [get_mu(10**-6,temp)-en_o-en_o_zpe for temp in temps]
# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
for sp in ax2.spines.itervalues():
    sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(temps)
ax2.set_xlim([mu[0]-(en_o+en_o_zpe), mu[1]-(en_o+en_o_zpe)])



#add text
ax.text(0.01, 0.01, 'Bulk Zn',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes, fontsize=16)

ax.text(0.92, 0.01, 'Bulk Cu$_2$O',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=16)

ax.text(0.57, 0.01, 'Surface \n Oxide Cu',
        verticalalignment='bottom', horizontalalignment='center',
        transform=ax.transAxes, fontsize=16)

ax.text(0.32, 0.01, 'Zn \n Oxides',
        verticalalignment='bottom', horizontalalignment='center',
        transform=ax.transAxes, fontsize=16)



# plot-p-T diagram
ax=plt.subplot(122)
temp_range=np.linspace(400, 1000, resolution)

pressure_range=np.power(10,np.linspace(-35, 1, resolution))*10**5

z=np.zeros((len(temp_range), len(pressure_range)))
mu_z=np.zeros((len(temp_range), len(pressure_range)))
levels=[int(energy[0][0]-1)]
for i,T in enumerate(temp_range):
    for j,p in enumerate(pressure_range):
        mu_z[i][j]=get_mu(p,T)
        z[i][j]=int(get_lowest(mu_z[i][j])[0])
        if int(get_lowest(mu_z[i][j])[0]) not in levels:
            levels.append(int(get_lowest(mu_z[i][j])[0]))

mu_z=np.transpose(mu_z)
z=np.transpose(z)
levels.sort()

color_levels=['' for level in levels[1:]]
for i,level in enumerate(levels[1:]):
    color_levels[i]=colors[level]
    
plt.contourf(temp_range,pressure_range,z, levels=levels, colors=color_levels, alpha=alpha_back)

plt.yscale('log')
plt.ylim([pressure_range[0],pressure_range[-1]])
plt.xlim([temp_range[0],temp_range[-1]])


levels=[zno_decomp, surf_alloy, bulk_cu2o]
mpl.rcParams['contour.negative_linestyle'] = 'solid'
plt.contour(temp_range, pressure_range, mu_z, colors='w', levels=levels, linestyle='-', linewidths=3)

#add text

ax.text(0.01, 0.9, 'Bulk Cu$_2$O',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes, fontsize=16)


ax.text(0.92, 0.01, 'Bulk Zn',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=16)


ax.text(0.92, 0.5, 'Zn Oxides',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=16)

ax.text(0.92, 0.8, 'Surface \n Oxide Cu',
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=16)


ax.set_yticklabels([int(tick) for tick in np.log10(ax.get_yticks())])
plt.ylabel('log$_{10}$($\\frac{p}{1 \mathrm{Pa}}$)', fontsize=16)
plt.xlabel('Temperature [K]', fontsize=16)
plt.show()
