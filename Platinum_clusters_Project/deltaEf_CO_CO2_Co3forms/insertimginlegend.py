import os
from ase.vibrations import Infrared
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo
from ase.io import read,write
import matplotlib.pyplot as plt
from gpaw import GPAW, FermiDirac, PoissonSolver, Mixer,PW
from ase.calculators.vasp import Vasp
from matplotlib.pyplot import *
import numpy as np
from matplotlib.transforms import TransformedBbox
from matplotlib.image import BboxImage
from matplotlib.legend_handler import HandlerBase
from matplotlib._png import read_png

class ImageHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):

        # enlarge the image by these margins
        sx, sy = self.image_stretch 

        # create a bounding box to house the image
        bb = Bbox.from_bounds(xdescent - sx,
                              ydescent - sy,
                              width + sx,
                              height + sy)

        tbb = TransformedBbox(bb, trans)
        image = BboxImage(tbb)
        image.set_data(self.image_data)

        self.update_prop(image, orig_handle, legend)

        return [image]

    def set_image(self, image_path, image_stretch=(0, 0)):
        if not os.path.exists(image_path):
            sample = get_sample_data("Pt7Oxides_template_CO2form1.png", asfileobj=False)
            self.image_data = read_png(sample)
        else:
            self.image_data = read_png(image_path)

        self.image_stretch = image_stretch

# random data
x = np.random.randn(100)
y = np.random.randn(100)
y2 = np.random.randn(100)

# plot two series of scatter data
s = plt.scatter(x, y, c='b')
s2 = plt.scatter(x, y2, c='r')

# setup the handler instance for the scattered data
custom_handler = ImageHandler()
custom_handler.set_image("[/home/siva/JCP_PtOxide_Al2O3_0001surface/deltaEf_CO_CO2_Co3forms/Pt7Oxides_template_CO2form1.png]",
                         image_stretch=(0, 20)) # this is for grace hopper

# add the legend for the scattered data, mapping the
# scattered points to the custom handler
plt.legend([s, s2],
           ['Scatters 1', 'Scatters 2'],
           handler_map={s: custom_handler, s2: custom_handler},
           labelspacing=2,
           frameon=False)

plt.show()
