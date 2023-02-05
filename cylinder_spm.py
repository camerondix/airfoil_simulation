import math
import panelGeometry
import plotting
import panelMethods
import matplotlib.pyplot as plt

# Source Panel Method of a Custom Cylinder

# Inputs
freestreamVelocity = 10
alpha = 0
radius = 1
divisions = 24


# Convert alpha to radians
alpha = alpha * math.pi / 180

# Create and plot panels
panels = panelGeometry.createCirclePanels(radius, divisions, alpha)
panels.reverse()

# Compute coefficients
cps, cd, cl, accuracy = panelMethods.findSourcePanelCoefficients(
    panels, freestreamVelocity, alpha
)

# Plot the points and panels
plt.subplot(1, 2, 1)
plotting.plotPanelsAndCps(panels, cps)
plt.xlabel("x")
plt.ylabel("y")

# Plot the cps
plt.subplot(1, 2, 2)
plotting.plotCircleAnalyticalCps(200)
plotting.plotCpsForCircle(panels, cps)
plt.title(str.format("$\sum \lambda _j S_j = {:.2e}$", accuracy))

plt.suptitle(str(divisions) + " Panel Cylinder")

plt.show()
