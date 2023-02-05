import math
import panelGeometry
import plotting
import panelMethods
import matplotlib.pyplot as plt
import os

# Source/Vortex Panel Method of an Imported Body

# Inputs
freestreamVelocity = 1
alphaDeg = 0
fileName = "NACA-2412_Geom.txt"  # The data file must be in the same folder as this file.
#   fileName = "NACA_0012_b.txt"  # The data file must be in the same folder as this file.
#   fileName = "Cyl_Geom.txt"  # The data file must be in the same folder as this file.
experimentalFileName = "NACA_0012_experimental.txt"  # https://ntrs.nasa.gov/api/citations/19880009181/downloads/19880009181.pdf
seperator = " "  # The seperator ie comma, space etc.

# Convert alpha to radians
alpha = alphaDeg * math.pi / 180

# Import the data from the specified file and create points/panels.
path = os.path.join(os.getcwd(), fileName)
points = panelGeometry.importPoints(path, seperator)
panels = panelGeometry.createPanelsFromPoints(points, alpha)

# Compute coefficients
cps, cl, cd, cm = panelMethods.findSourceVortexPanelCoefficients(
    panels, freestreamVelocity, alpha
)

# Plot the imported points and panels
plt.subplot(2, 1, 1)
plotting.plotPanelsAndCps(panels, cps)
plt.xlabel("x")
plt.ylabel("y")
plt.title(r'$\alpha = $' + str(alphaDeg) + u'\N{DEGREE SIGN}')

# Plot the cps
plt.subplot(2, 1, 2)
plotting.plotCps(points, cps, label=(r'$\alpha = $' + str(alphaDeg) + u'\N{DEGREE SIGN}'))
plt.title(str.format("$C_L = {:.3}$", cl))

# Plot experimental cps for naca 0012
if fileName == "NACA_0012_b.txt" and alpha == 0:
    path = os.path.join(os.getcwd(), experimentalFileName)
    xs, cps = panelGeometry.importTuples(path, seperator)
    plt.plot(xs, cps, "*", c="k")

panels = panelGeometry.createPanelsFromPoints(points, alpha=(0 * math.pi / 180))
cps, cl, cd, cm = panelMethods.findSourceVortexPanelCoefficients(
    panels, freestreamVelocity, alpha=(0 * math.pi / 180)
)

plotting.plotCps(points, cps, color="green", label=(r'$\alpha = $' + str(0) + u'\N{DEGREE SIGN}'))
plt.gca().invert_yaxis()

panels = panelGeometry.createPanelsFromPoints(points, alpha=(10 * math.pi / 180))
cps, cl, cd, cm = panelMethods.findSourceVortexPanelCoefficients(
    panels, freestreamVelocity, alpha=(10 * math.pi / 180)
)

plotting.plotCps(points, cps, color="red", label=(r'$\alpha = $' + str(10) + u'\N{DEGREE SIGN}'))
plt.gca().invert_yaxis()

plt.suptitle(os.path.splitext(fileName)[0])

plt.legend(loc="upper right")

plt.show()
