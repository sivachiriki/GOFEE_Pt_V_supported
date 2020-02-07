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


#  formation energy of oxygen
Pt7Ox=[[0, 0.0],
[2, -3.228179134242799],
[4, -6.222445977342442],
[6, -8.691274797583972],
[8, -9.618859668972185],
[10, -10.713886005601097],
[12, -10.943386683601517],
[14, -11.218929335924585],
[16, -9.866595826652471]]
print(np.shape(Pt7Ox))
print(len(Pt7Ox))
Pt7Ox = np.array(Pt7Ox)
print((Pt7Ox[:,0]))
#chemisorption of CO molecule 
Pt7OxCO=[[0, -2.6255645259353013],
[2, -5.9666031289775585],
[4, -8.460640759613906],
[6, -10.23591680017853],
[8, -10.702740792354184],
[10, -11.37352059698732],
[12, -11.582457655820079],
[14, -11.852233422852834]]
Pt7OxCO=np.array(Pt7OxCO)
#chemisorption of CO2 molecule attached to the cluster terminal
Pt7OxCO2=[[2,-5.988002673455972],
[4, -8.621413739587236],
[6, -10.543490464432352],
[8, -11.508377394322402],
[10, -12.39804226918558],
[12, -13.352578716735778],
[14, -14.377973328460513],
[16, -11.946205594281054]]
Pt7OxCO2 = np.array(Pt7OxCO2)
# chemisorption of CO3 molecule form on within cluster
Pt7OxCO3=[[2, -5.043113192502336],
[4,-7.856083072271016],
[6,  -10.031012823815558],
[8, -12.233746351882424],
[10, -13.887649384911711],
[12, -14.14148595635477],
[14, -14.84690225602074],
[16, -13.057021308935617]]
Pt7OxCO3=np.array(Pt7OxCO3)
# chemisorption of CO2 dettached from cluster
Pt7OxCO2detached=[[2,  -5.092979606472898],
[4, -8.091406898979184],
[6, -10.461499640039165],
[8, -11.880547430604516],
[10, -13.376545594653308],
[12, -14.277087877371812],
[14, -15.018589767618806],
[16, -12.45595536962648]]
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
fig, ax = plt.subplots(1,1,figsize=(11,9))
#ax.plot(Pt7OxCO[:,0],delta_E_f_CO,marker="o",color='#0f0f0f',label='Pt$_7$O$_y$CO')
coattached, =ax.plot(Pt7OxCO[:,0],delta_E_f_CO,marker="o",color='#0f0f0f')
for i in range(len(Pt7OxCO)):
    ax.plot(Pt7OxCO[i,0],delta_E_f_CO[i],marker="o",color='#0f0f0f')
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
#ax.plot(Pt7OxCO2[:,0],delta_E_f_CO2,marker="o",color='#8EBA42',label='Pt$_7$O$_{(y-1)}$CO$_2$')
co2attached, =ax.plot(Pt7OxCO2[:,0],delta_E_f_CO2,marker="o",color='#8EBA42')
for i in range(len(Pt7OxCO2)):
    ax.plot(Pt7OxCO2[i,0],delta_E_f_CO2[i],marker="o",color='#8EBA42')
############################### CO2 detached #####################################################
colors = {}
color_lib = ['#377eb8','#4daf4a','#984ea3','#999999','#fdbf6f', '#ff7f00','#eeefff','#ffff33']
for i,atoms in enumerate(Pt7OxCO2detached):
    colors[i] = color_lib[i]
delta_E_f_CO2detached =np.zeros(len(Pt7OxCO2detached))
for i in range(len(Pt7OxCO2detached)):
    delta_E_f_CO2detached[i] =(Pt7OxCO2detached[i,1]-Pt7Ox[i+1,1])
    print(Pt7OxCO2detached[i,0],delta_E_f_CO2detached[i])
# plot formation E vs No. of Oxygens
co2detached, = plt.plot(Pt7OxCO2detached[:,0],delta_E_f_CO2detached,marker="o",color='#FF00FF')
for i in range(len(Pt7OxCO2detached)):
    ax.plot(Pt7OxCO2detached[i,0],delta_E_f_CO2detached[i],marker="o",color='#FF00FF')
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
#ax.plot(Pt7OxCO3[:,0],delta_E_f_CO3,marker="o",color='#FF0000',label='Pt$_7$O$_{(y-2)}$CO$_3$')
co3attached, = plt.plot(Pt7OxCO3[:,0],delta_E_f_CO3,marker="o",color='#FF0000')
for i in range(len(Pt7OxCO3)):
    ax.plot(Pt7OxCO3[i,0],delta_E_f_CO3[i],marker="o",color='#FF0000')
##################################################################################################
plt.legend([coattached,co2attached,co2detached,co3attached], ["","","", ""],
   handler_map={coattached:HandlerLineImage("Pt7Oxides_template_COform.png"),co2attached:HandlerLineImage("Pt7Oxides_template_CO2form1.png"),co2detached: HandlerLineImage("Pt7Oxides_template_CO2form4.png"), co3attached: HandlerLineImage("Pt7Oxides_template_CO3form.png")},
   handlelength=3.50, labelspacing=0.0, fontsize=40, borderpad=0.2, loc=1,
    handletextpad=0.0, borderaxespad=0.1)
################################################################################################
ax.set_ylabel(r'Extra Stability of Pt$_7$-Oxides with CO [$\Delta$E$_f$ (eV)]')
ax.set_xlabel(r'No. of O')
ax.set_xticks(np.arange(0, 22, step=1))
#plt.legend()
savefig('delta_COEf_CO2Ef_CO3Ef_O2Ef_Pt7Oxides.png')
plt.show()
