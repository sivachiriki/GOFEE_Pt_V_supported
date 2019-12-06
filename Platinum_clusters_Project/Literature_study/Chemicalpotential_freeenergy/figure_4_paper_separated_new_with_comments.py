import numpy as np
from ase import * 
from ase.io import read,write
import matplotlib.pyplot as plt
import operator
from matplotlib import lines
from kmos import evaluate_rate_expression
from matplotlib.ticker import AutoMinorLocator, AutoLocator, MultipleLocator

l = [1,2,3,4,5,6]
#l = [1,2,3,4,5]


slab_sr = []
h2o_sr = []
hoh_sr = []
h_sr = []
vac_sr = []

slab_sp = []
h2o_sp = []
hoh_sp = []
h_sp = []
vac_sp = []

slab_spu = []
h2o_spu = []
hoh_spu = []
h_spu = []
vac_spu = []

hohb_spu = []

for i in range(len(l)):

    slab_sr.append(read('adm_{}x4_source_fd_done_sr.traj'.format(l[i])).get_potential_energy())
    h2o_sr.append(read('adm_{}x4_H2O_fd_done_sr.traj'.format(l[i])).get_potential_energy()) 
    hoh_sr.append(read('adm_{}x4_HOH_fd_done_sr.traj'.format(l[i])).get_potential_energy())
    
    slab_sp.append(read('adm_{}x4_source_fd_sp_done.traj'.format(l[i])).get_potential_energy())
    h2o_sp.append(read('adm_{}x4_H2O_fd_sp_done.traj'.format(l[i])).get_potential_energy()) 
    hoh_sp.append(read('adm_{}x4_HOH_fd_sp_done.traj'.format(l[i])).get_potential_energy())
    
    slab_spu.append(read('adm_{}x4_sourcepbeu_sp_done.traj'.format(l[i])).get_potential_energy())
    h2o_spu.append(read('adm_{}x4_H2Opbeu_sp_done.traj'.format(l[i])).get_potential_energy()) 
    hoh_spu.append(read('adm_{}x4_HOHpbeu_sp_done.traj'.format(l[i])).get_potential_energy()) 
    hohb_spu.append(read('adm_{}x4_HOHbpbeu_sp_done.traj'.format(l[i])).get_potential_energy()) 

e_o2_pbe = read('o2fd_sp_done.traj').get_potential_energy()
e_h2o_pbe = read('h2ofd_sp_done.traj').get_potential_energy()
e_h2_pbe = read('h2fd_done.traj').get_potential_energy()

wi, le = 3.5,3.5*0.75
fig1,((ax1)) = plt.subplots(1,1,figsize=(wi,le))
fig2,((ax2)) = plt.subplots(1,1,figsize=(wi,le))
fig3,((ax3)) = plt.subplots(1,1,figsize=(wi,le))
fig3,((ax4)) = plt.subplots(1,1,figsize=(wi,le))

hohb_2x3 = read('adm_2x3x4_HOHbpbeu_sp_done.traj').get_potential_energy()
hohb_2x2 = read('adm_2x2x4_HOHbpbeu_sp_done.traj').get_potential_energy()

ads_aux1 = (np.array(hohb_2x3) - np.array(slab_spu[5]) - 2*e_h2o_pbe)/2.
ads_aux2 = (np.array(hohb_2x2) - np.array(slab_spu[3]) - 2*e_h2o_pbe)/2.

# this calculates the changes in surface energy upon h2o adsorption
# which is the adsorption energy normalized by the surface area

n_ads_hohb = (np.array(hohb_spu) - np.array(slab_spu) - e_h2o_pbe)/np.array(l)



ads_hohb = (np.array(hohb_spu) - np.array(slab_spu) - e_h2o_pbe)

eds = []
eds.append(ads_hohb[5])
eds.append((ads_hohb[4]*12.-ads_hohb[5]*10)/2.)
eds.append((ads_hohb[3]*15.-ads_hohb[4]*12.)/3.)
eds.append((ads_hohb[2]*20.-ads_hohb[3]*15.)/5.)
eds.append((ads_hohb[1]*30.-ads_hohb[2]*20.)/10.)
eds.append((ads_hohb[0]*60.-ads_hohb[1]*30.)/30.)


mumin = -1.5
mumax = 0.2

mu = np.linspace(mumin, mumax,200)
ratios = [0., 1./6, 1./5., 1/4., 1./3, 1./2, 1.]
ratios = np.array(ratios)



#  we define the array states, which includes the surface energy changes for each state
# we also add a state for the bare surface, for which the surface energy change is zero
states = np.insert(n_ads_hohb[::-1], 0, 0.)



# This defines  the array that will contain the delta gamma variations as a function of mu
deltas=np.zeros((len(states),len(mu)))

colors = [
'#a6cee3',
'#1f78b4',
'#b2df8a',
'#33a02c',
'#fb9a99',
'#e31a1c',
'#fdbf6f',
]



# this calculates and plots the evolution of the delta gamma for each state
for i in range(len(states)):
    for j in range(len(mu)):
       deltas[i,j] = -ratios[i]*mu[j]/4. + states[i]/4.
    ax4.plot(mu,deltas[i],colors[i])


emin,emax = -0.3, 0.1
ax4.axis(xmin=mumin, xmax=mumax, ymin=emin, ymax=emax)


low_g_mo = np.zeros(len(mu))
low_g_mo_ind = np.zeros(len(mu))
mu_limits_mo = []
alfa_l = 0.5

print ''

# This identifies the most stable state for each value of my and colors the plot accordingly

fin = None
for j in range(len(mu)):
    if mu[j] >= mumin and len(mu_limits_mo)==0:
        mu_limits_mo.append(j)
    low_g_mo_ind[j], low_g_mo[j] = min(enumerate(deltas[i][j] for i in range(len(states))),key=operator.itemgetter(1))
    if low_g_mo_ind[j] != low_g_mo_ind[j-1]:
        mu_limits_mo.append(j)
    if mu[j] >= mumax and fin==None:
        mu_limits_mo.append(j)
        fin = 1
for i in range(len(mu_limits_mo)-1):
    lim0,lim1 = mu_limits_mo[i], mu_limits_mo[i+1]
    ind1 = int(low_g_mo_ind[(lim0+lim1)/2])
    print lim0, lim1, ind1
    ax4.plot(mu[lim0:lim1],low_g_mo[lim0:lim1],color=colors[ind1],lw=3.0)
    ax4.fill_between([mu[lim0],mu[lim1]],[emin,emin],[emax,emax],color = colors[ind1],alpha=alfa_l)




#hohb ads
ms = 5
ax1.plot(l,ads_hohb[::-1], 'b-o', label='OH+H abs')
print l
print ads_hohb[::-1]
#ax1.plot(3,ads_aux1,'^r',markersize=ms)
#ax1.plot(2,ads_aux2,'^r',markersize=ms)
ax2.plot(l,n_ads_hohb[::-1], 'r-o', label='OH+H norm')
#ax2.plot(3,ads_aux1/3.,'^b',markersize=ms)
#ax2.plot(2,ads_aux2/2.,'^b',markersize=ms)
print l
l5 = range(1,6)
print l5
ax3.plot(l,eds, 'k-o', label='OH+H abs')
labs_seq = [
r'0 $\rightarrow$ 1/6',
r'1/6 $\rightarrow$ 1/5',
r'1/5 $\rightarrow$ 1/4',
r'1/4 $\rightarrow$ 1/3',
r'1/3 $\rightarrow$ 1/2',
r'1/2 $\rightarrow$ 1/1',
]


print labs_seq

y_ticks = np.arange(0,-3,1.0)

plt.sca(ax3)
plt.xticks(l, labs_seq,rotation=30, fontsize = 7)
ax3.yaxis.set_major_locator(MultipleLocator(1.0))
ax3.yaxis.set_minor_locator(MultipleLocator(0.5))


labs_seq = [
r'0 $\rightarrow$ 1/6',
r'0 $\rightarrow$ 1/5',
r'0 $\rightarrow$ 1/4',
r'0 $\rightarrow$ 1/3',
r'0 $\rightarrow$ 1/2',
r'0 $\rightarrow$ 1/1',
]


#ax1.set_xticklabels(l, labs_seq)
#for tick in ax1.get_xticklabels():
#    tick.set_rotation(45)
#
plt.sca(ax1)
plt.xticks(l, labs_seq,rotation=30,fontsize=7)
ax1.yaxis.set_major_locator(MultipleLocator(1.0))
ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
plt.sca(ax2)
plt.xticks(l, labs_seq,rotation=30, fontsize=7)
ax2.yaxis.set_major_locator(MultipleLocator(.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))



ax1.set_ylabel(r'E$_{ads}$ (H$_{2}$O) [eV]')
ax3.set_ylabel(r'$\Delta$E$_{ads}$ (H$_{2}$O) [eV]')
ax2.set_ylabel(r'E$_{ads}$(H$_{2}$O)$\cdot$$\Theta$ [eV]')
ax4.set_ylabel(r'$\Delta \gamma$ [eV / unit cell]')
ax4.set_xlabel(r'$\Delta \mu_{\mathrm{H_{2}O}}$ [eV]')

cmin,cmax = 0.5, 6.5
amin,amax = -3.7, 0.5
amin2, amax2 = -1.0, -0.5

ax1.axis(xmin=cmin, xmax=cmax, ymin=amin, ymax= amax)
ax2.axis(xmin=cmin, xmax=cmax, ymin=amin2, ymax= amax2)
ax3.axis(xmin=cmin, xmax=cmax, ymin=amin, ymax= amax)

ax4.text(-1.3,-0.25,r'$\Theta$=1/4')
ax4.text(-0.50,-0.25,r'$\Theta$=1/3')
ax4.text(-0.15,-0.05,r'$\Theta$=1')



ax4.plot([-0.580 , -0.580], [-0.5 , 0.5], 'k--', alpha=0.4)
ax4.yaxis.set_major_locator(MultipleLocator(.1))
ax4.yaxis.set_minor_locator(MultipleLocator(0.05))


scale1 = ax4.twiny()
#scale2 = ax4.twiny()

scale1.set_xlim(mumin, mumax)
plt.sca(scale1)

labs_seq = []
labs_seq =  [1e-18,1e-14,1e-10,1e-6,1e-2,1e2,1e6,1e10]
labs_seq_strings = []

marks = []
for T in [298.15]:
#for T in [200]:
   for p in labs_seq:
       labs_seq_strings.append(r'$10^{%d}$'%(np.log10(p)))
       print  evaluate_rate_expression('mu_H2Ogas',
              {'T':{'value':T},
               'p_H2Ogas': {'value':p}} )

       value= evaluate_rate_expression('mu_H2Ogas',
              {'T':{'value':T},
               'p_H2Ogas': {'value':p}} )

       marks.append(value)



plt.xticks(marks, labs_seq_strings)
scale1.tick_params(labelsize =7)
scale1.set_xlabel('$p\mathrm{_{H_{2}O}}$ $\mathrm{[bar], RT}$', fontsize = 8)
plt.tight_layout()


print ads_hohb
#ax1.legend()
#ax2.legend()

n=1
for ax in [ax1,ax2,ax3,ax4]:
    plt.sca(ax)
    plt.tight_layout()
    plt.savefig('panel{}.pdf'.format(n),transparent=True)
    #plt.savefig('h2o_ridge_new.eps',transparent=True)
    n +=1 

plt.show()
