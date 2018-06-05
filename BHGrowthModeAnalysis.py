from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.logspace(-4,-1)
Y = np.logspace(-4,-1)
X, Y = np.meshgrid(X, Y)
Z = np.log10((2-57*Y)*(30/57)*X/Y)

# Plot the surface.
ax.plot_surface(X,Y,Z,linewidth=0, antialiased=False)

# Customize the z axis.
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.view_init(90,0)
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
