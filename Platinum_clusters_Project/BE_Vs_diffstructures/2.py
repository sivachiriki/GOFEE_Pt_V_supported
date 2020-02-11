import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes


x = np.arange(-10, 10, 0.01)
sinx = np.sin(x)
tanx = np.tan(x)

fig, ax = plt.subplots( 1, 2, sharey='row', figsize=(9, 3) )

for i, f in enumerate([sinx, tanx]):
    ax[i].plot( x, f, color='red' )
    ax[i].set_ylim([-2, 2])

    # create an inset axe in the current axe:
    inset_ax = inset_axes(ax[i],
                          height="30%", # set height
                          width="30%", # and width
                          loc=10) # center, you can check the different codes in plt.legend?
    inset_ax.plot(x, f, color='green')
    inset_ax.set_xlim([0, 5])
    inset_ax.set_ylim([0.75, 1.25])
plt.show()
