import math
import panelGeometry
import plotting
import panelMethods
import matplotlib.pyplot as plt
import os

# Source Panel Method of an Imported Body

# Inputs
freestreamVelocity = 1
alpha = 0
fileName = "NACA_0012_b.txt"  # The data file must be in the same folder as this file.
#   fileName = "Cyl_Geom.txt"  # The data file must be in the same folder as this file.
experimentalFileName = "NACA_0012_experimental.txt"  # https://ntrs.nasa.gov/api/citations/19880009181/downloads/19880009181.pdf
seperator = " "  # The seperator ie comma, space etc.

# Convert alpha to radians
alpha = alpha * math.pi / 180

# Import the data from the specified file and create points/panels.
path = os.path.join(os.getcwd(), fileName)
points = panelGeometry.importPoints(path, seperator)
panels = panelGeometry.createPanelsFromPoints(points, alpha)

# Compute coefficients
cps, cl, cd, accuracy = panelMethods.findSourcePanelCoefficients(
    panels, freestreamVelocity, alpha
)

# Plot the imported points and panels
plt.subplot(2, 1, 1)
plotting.plotPanelsAndCps(panels, cps)
plt.xlabel("x")
plt.ylabel("y")

# Plot the cps
plt.subplot(2, 1, 2)
plotting.plotCps(points, cps)
plt.title(str.format("$\sum \lambda _j S_j = {:.2e}$", accuracy))

# Plot experimental cps for naca 0012
if fileName == "NACA_0012_b.txt" and alpha == 0:
    path = os.path.join(os.getcwd(), experimentalFileName)
    xs, cps = panelGeometry.importTuples(path, seperator)
    plt.plot(xs, cps, "*", c="k")

plt.suptitle(os.path.splitext(fileName)[0])

plt.show()
