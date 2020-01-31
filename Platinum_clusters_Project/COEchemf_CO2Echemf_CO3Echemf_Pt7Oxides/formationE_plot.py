from ase.vibrations import Infrared
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.io import read,write
import matplotlib.pyplot as plt
import matplotlib.lines
from matplotlib.transforms import Bbox, TransformedBbox
from matplotlib.legend_handler import HandlerBase
from matplotlib.image import BboxImage
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from ase.calculators.vasp import Vasp
from matplotlib.pyplot import *
import numpy as np
import os

class HandlerLineImage(HandlerBase):

    def __init__(self, path, space=15, offset = 10 ):
        self.space=space
        self.offset=offset
        self.image_data = plt.imread(path)
        super(HandlerLineImage, self).__init__()

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):

        l = matplotlib.lines.Line2D([xdescent+self.offset,xdescent+(width-self.space)/3.+self.offset],
                                     [ydescent+height/2., ydescent+height/2.])
        l.update_from(orig_handle)
        l.set_clip_on(False)
        l.set_transform(trans)

        bb = Bbox.from_bounds(xdescent +(width+self.space)/3.+self.offset,
                              ydescent,
                              height*self.image_data.shape[1]/self.image_data.shape[0],
                              height)

        tbb = TransformedBbox(bb, trans)
        image = BboxImage(tbb)
        image.set_data(self.image_data)

        self.update_prop(image, orig_handle, legend)
        return [l,image]
################################# E_form of O2 ####################################################
# read stru of GM Pt7 oxides
O2 = read('O2_kpts221.traj')
Pt7O0 =  read('Pt7_Al2O3_Planar_GMDFTrelaxed.traj')
Pt7O2 =  read('Pt7O2_Al2O3_chem_GMDFTrelaxed.traj')
Pt7O4 =  read('Pt7O4_Al2O3_chem_GMDFTrelaxed.traj')
Pt7O6 =  read('Pt7O6_Al2O3_chem_DFTrelaxed.traj')
Pt7O7 =  read('Pt7O7_Al2O3_Chem_GMDFTrelaxed.traj')
Pt7O8 =  read('Pt7O8_Al2O3_Chem_GMDFTrelaxed.traj')
Pt7O10 = read('Pt7O10_Al2O3_Chem_GMDFTrelaxed.traj')
Pt7O12 = read('Pt7O12_Al2O3_Chem_GMDFTrelaxed.traj')
Pt7O14 = read('Pt7O14_Al2O3_Chem_GMDFTrelaxed.traj')
Pt7O16 = read('Pt7O16_Al2O3_chem_GMDFTrelaxed.traj')
# energies of GM stru of Pt7 oxides
E_O2 = O2.get_potential_energy()
E_Pt7O0 = Pt7O0.get_potential_energy()
E_Pt7O2 = Pt7O2.get_potential_energy()
E_Pt7O4 = Pt7O4.get_potential_energy()
E_Pt7O6 = Pt7O6.get_potential_energy()
E_Pt7O7 = Pt7O7.get_potential_energy()
E_Pt7O8 = Pt7O8.get_potential_energy()
E_Pt7O10 = Pt7O10.get_potential_energy()
E_Pt7O12 = Pt7O12.get_potential_energy()
E_Pt7O14 = Pt7O14.get_potential_energy()
E_Pt7O16 = Pt7O16.get_potential_energy()

gm_tru =[Pt7O0,Pt7O2,Pt7O4,Pt7O6,Pt7O7,Pt7O8,Pt7O10,Pt7O12,Pt7O14,Pt7O16]
N_O = [0,2,4,6,7,8,10,12,14,16]
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
PtOxides, = ax.plot(N_O,E_f_correct,marker="o",color='#1E90FF')
for i, nu in enumerate(N_O):
    #ax.plot(N_O[i],E_f_correct[i],marker="o",color=colors[i])
    ax.plot(N_O[i],E_f_correct[i],marker="o",color='#1E90FF')

################################# E_ads of CO ####################################################
CO = read('CO_molecule_DFTrelaxed.traj')
E_CO = CO.get_potential_energy()
Pt7O0CO = read('Pt7CO_Chem_COform_Al2O3_GMDFTrelaxed.traj')
Pt7O2CO = read('Pt7O2CO_Al2O3_chem_COform_GMDFTrelaxed.traj')
Pt7O4CO = read('Pt7O4CO_Al2O3_chem_COform_LM2DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6CO_Al2O3_Chem_COform_LM2DFTrelaxed.traj')
Pt7O8CO = read('Pt7O8CO_Al2O3_Chem_COform_LM1DFTrelaxed.traj')
Pt7O10CO = read('Pt7O10CO_Al2O3_Chem_COfrom_LM2DFTrelaxed.traj')
Pt7O12CO = read('Pt7O12CO_Al2O3_Chem_COform_LM4DFTrelaxed.traj')
Pt7O14CO = read('Pt7O14CO_Al2O3_Chem_COform_LM3DFTrelaxed.traj')
#Pt7O16CO = read('Pt7O16CO_Al2O3_Chem_COform_LM3DFTrelaxed.traj')

E_Pt7O0CO = Pt7O0CO.get_potential_energy()
E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()
#E_Pt7O16CO = Pt7O16CO.get_potential_energy()

gm_tru_CO =[Pt7O0CO,Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO]
print(gm_tru_CO[-1].get_potential_energy())
N_O = [0,2,4,6,8,10,12,14]
# color labels
colors = {}
color_lib =['#bf40bf','#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
coattached, =ax.plot(N_O,E_ads,marker="o",color='#0f0f0f')
for i, nu in enumerate(N_O):
    #ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])
    ax.plot(N_O[i],E_ads[i],marker="o",color='#0f0f0f')
#################################### E_chem CO2 attached to edge of cluster #######################
CO = read('CO_molecule_DFTrelaxed.traj')
E_CO = CO.get_potential_energy()
Pt7O2CO = read('Pt7O2CO_Al2O3_chem_CO2form_LMDFTrelaxed.traj')
Pt7O4CO = read('Pt7O4CO_Al2O3_chem_CO2form_LM1DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6CO_Al2O3_Chem_CO2form_GMDFTrelaxed.traj')
Pt7O8CO = read('Pt7O8CO_Al2O3_Chem_CO2form_LM2DFTrelaxed.traj')
Pt7O10CO = read('Pt7O10CO_Al2O3_Chem_CO2from_LM1DFTrelaxed.traj')
Pt7O12CO = read('Pt7O12CO_Al2O3_Chem_CO2form_LM1DFTrelaxed.traj')
Pt7O14CO = read('Pt7O14CO_Al2O3_Chem_CO2form_LMDFTrelaxed.traj')

E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()

gm_tru_CO =[Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO]
print(gm_tru_CO[-1].get_potential_energy())
N_O = [2,4,6,8,10,12,14]
# color labels
colors = {}
color_lib =['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
co2attached, = ax.plot(N_O,E_ads,marker="o",color='#8EBA42')
for i, nu in enumerate(N_O):
    #ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])
    ax.plot(N_O[i],E_ads[i],marker="o",color='#8EBA42')
################################## E_chem of CO2 detached from cluster ###########################
CO = read('CO_molecule_DFTrelaxed.traj')
E_CO = CO.get_potential_energy()
Pt7O2CO = read('Pt7O2CO_Al2O3_chem_CO2formonsurface_LMDFTrelaxed.traj')
Pt7O4CO = read('Pt7O4CO_Al2O3_chem_CO2formonsurface_LM3DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6CO_Al2O3_Chem_CO2formonsurface_LM1DFTrelaxed.traj')
Pt7O8CO = read('Pt7O8CO_Al2O3_Chem_CO2formonsurface_LMDFTrelaxed.traj')
Pt7O10CO = read('Pt7O10CO_Al2O3_Chem_CO2fromonsurface_LMnewDFTrelaxed.traj')
Pt7O12CO = read('Pt7O12CO_Al2O3_Chem_CO2formonsurface_LM2DFTrelaxed.traj')
Pt7O14CO = read('Pt7O14CO_Al2O3_Chem_CO2formonsurface_LM2DFTrelaxed.traj')
Pt7O16CO = read('Pt7O16CO_Al2O3_Chem_CO2formonsurface_LMDFTrelaxed.traj')

E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()
E_Pt7O16CO = Pt7O16CO.get_potential_energy()

gm_tru_CO =[Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO,Pt7O16CO]
print(gm_tru_CO[-1].get_potential_energy())
N_O = [2,4,6,8,10,12,14,16]
# color labels
colors = {}
color_lib =['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f','#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
co2detached, = ax.plot(N_O,E_ads,marker="o",color='#FF00FF')
for i, nu in enumerate(N_O):
    #ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])
    ax.plot(N_O[i],E_ads[i],marker="o",color='#FF00FF')
##################################################################################################
########################## E_chem of CO3 (carbene) ###############################################
CO = read('CO_molecule_DFTrelaxed.traj')
E_CO = CO.get_potential_energy()
Pt7O2CO = read('Pt7O2CO_Al2O3_chem_CO3form_LM2DFTrelaxed.traj')
Pt7O4CO = read('Pt7O4CO_Al2O3_chem_CO3form_LM4DFTrelaxed.traj')
Pt7O6CO = read('Pt7O6CO_Al2O3_Chem_CO3form_LMDFTrelaxed.traj')
#Pt7O7CO = read('Pt7O7CO_Al2O3_chem_CO3form_GMDFTrelaxed.traj')
Pt7O8CO = read('Pt7O8CO_Al2O3_Chem_CO3form_GMDFTrelaxed.traj')
Pt7O10CO = read('Pt7O10CO_Al2O3_Chem_CO3from_GMnewDFTrelaxed.traj')
Pt7O12CO = read('Pt7O12CO_Al2O3_Chem_CO3form_GMDFTrelaxed.traj')
Pt7O14CO = read('Pt7O14CO_Al2O3_Chem_CO3form_GMDFTrelaxed.traj')
Pt7O16CO = read('Pt7O16CO_Al2O3_Chem_CO3form_GMDFTrelaxed.traj')

E_Pt7O2CO = Pt7O2CO.get_potential_energy()
E_Pt7O4CO = Pt7O4CO.get_potential_energy()
E_Pt7O6CO = Pt7O6CO.get_potential_energy()
#E_Pt7O7CO = Pt7O7CO.get_potential_energy()
E_Pt7O8CO = Pt7O8CO.get_potential_energy()
E_Pt7O10CO = Pt7O10CO.get_potential_energy()
E_Pt7O12CO = Pt7O12CO.get_potential_energy()
E_Pt7O14CO = Pt7O14CO.get_potential_energy()
E_Pt7O16CO = Pt7O16CO.get_potential_energy()

gm_tru_CO =[Pt7O2CO,Pt7O4CO,Pt7O6CO,Pt7O8CO,Pt7O10CO,Pt7O12CO,Pt7O14CO,Pt7O16CO]
print(gm_tru_CO[-1].get_potential_energy())
N_O = [2,4,6,8,10,12,14,16]
# color labels
colors = {}
color_lib =['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(gm_tru_CO):
    colors[i] = color_lib[i]
# cal formation energy
E_ads =np.zeros(len(gm_tru_CO))
for i,atoms in enumerate(gm_tru_CO):
    E_ads[i] = (gm_tru_CO[i].get_potential_energy() - gm_tru[0].get_potential_energy() - CO.get_potential_energy() - (N_O[i])*(0.5000*O2.get_potential_energy()+Oxygen_correct))
    print(N_O[i],E_ads[i])
# plot formation E vs No. of Oxygens
co3attached, = ax.plot(N_O,E_ads,marker="o",color='#FF0000')
for i, nu in enumerate(N_O):
    #ax.plot(N_O[i],E_ads[i],marker="o",color=colors[i])
    ax.plot(N_O[i],E_ads[i],marker="o",color='#FF0000') 
##################################################################################################
plt.legend([PtOxides,coattached,co2attached,co2detached,co3attached], ["","","","", ""],
   handler_map={PtOxides:HandlerLineImage("Pt7Oxides_template.png"),coattached:HandlerLineImage("Pt7Oxides_template_COform.png"),co2attached:HandlerLineImage("Pt7Oxides_template_CO2form1.png"),co2detached: HandlerLineImage("Pt7Oxides_template_CO2form4.png"), co3attached: HandlerLineImage("Pt7Oxides_template_CO3form.png")},
   handlelength=3.50, labelspacing=0.0, fontsize=40, borderpad=0.2, loc=1,
    handletextpad=0.0, borderaxespad=0.1)
ax.set_ylabel(r'E$_f$ (eV)')
ax.set_xlabel(r'No. of O')
ax.set_xticks(np.arange(0, 22, step=1))
#plt.legend()
savefig('COEf_CO2Ef_CO3Ef_O2Ef_Pt7Oxides.png')
plt.show()
