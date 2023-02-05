import math
import panelGeometry
import plotting
import panelMethods
import matplotlib.pyplot as plt
import os

# Source/Vortex Panel Method For Varying Alphas

# Inputs
freestreamVelocity = 1
alphaMinDeg = -20
alphaMaxDeg = 20
fileName = "NACA-2412_Geom.txt"  # The data file must be in the same folder as this file.
#   fileName = "NACA_0012_b.txt"  # The data file must be in the same folder as this file.
#   fileName = "Cyl_Geom.txt"  # The data file must be in the same folder as this file.
#   experimentalFileName = "NACA_0012_cl_a.txt" #https://ntrs.nasa.gov/api/citations/19880019495/downloads/19880019495.pdf
experimentalFileName = "NACA-2412_Book-fig4-10.txt"
seperator = " "  # The seperator ie comma, space etc.

# Convert alpha to radians
alphaDegs = list(range(alphaMinDeg, alphaMaxDeg))
#   alphaDegs = list(0.01*deg for deg in alphaDegs)

# Import the data from the specified file and create points/panels.
path = os.path.join(os.getcwd(), fileName)
points = panelGeometry.importPoints(path, seperator)

# Compute the cl for various angle of attacks
cls = []
for alphaDeg in alphaDegs:
    alpha = alphaDeg * math.pi / 180

    panels = panelGeometry.createPanelsFromPoints(points, alpha)

    # Compute coefficients
    cps, cl, cd, cm = panelMethods.findSourceVortexPanelCoefficients(
        panels, freestreamVelocity, alpha
    )
    cls.append(cl)

# Plot the cls vs alpha
plotting.plotAlphaAndCls(alphaDegs, cls)

# Plot the experimental data
if experimentalFileName:
    path = os.path.join(os.getcwd(), experimentalFileName)
    alphs, cls = panelGeometry.importTuples(path, ",")
    plt.plot(alphs, cls, "*", c="k", label='Experimental')

plt.suptitle(os.path.splitext(fileName)[0])

plt.legend(loc="upper left")
plt.show()
