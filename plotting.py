import math
import matplotlib.pyplot as plt
from panelGeometry import Point, Panel


# Methods
def plotPoints(points: list):
    """
    Scatter plots a list of Points.
    """
    if isinstance(points, list):
        for point in points:
            plt.scatter(point.x, point.y, color="k")
    else:
        plt.scatter(points.x, points.y)
    setAxes()


def plotPanels(panels: list):
    """
    Plots a list of Panels as lines.
    """
    if isinstance(panels, list):
        for panel in panels:
            plt.plot(
                [panel.startPoint.x, panel.endPoint.x],
                [panel.startPoint.y, panel.endPoint.y],
                color="k",
            )
    else:
        plt.plot(
            [panels.startPoint.x, panels.endPoint.x],
            [panels.startPoint.y, panels.endPoint.y],
        )
    setAxes()


def plotBody(xVals: list, yVals: list):
    """
    Plots a body by connecting xVals and yVals in order.
    """
    for x, y in zip(xVals, yVals):
        plt.plot(x, y)


def plotPanelsAndCps(panels: list, cps: list):
    """
    Plots the panels with cp vectors at each control point.
    """
    plotPanels(panels)
    # plotPoints(list(panel.controlPoint for panel in panels))
    for i, panel in enumerate(panels):
        surfaceX = panel.controlPoint.x + 0.01 * math.cos(panel.delta)
        surfaceY = panel.controlPoint.y + 0.01 * math.sin(panel.delta)
        awayX = surfaceX + abs(cps[i]) * 0.1 * math.cos(panel.delta)
        awayY = surfaceY + abs(cps[i]) * 0.1 * math.sin(panel.delta)
        if cps[i] > 0:
            plt.arrow(awayX, awayY, surfaceX - awayX, surfaceY - awayY, color="red")
        else:
            plt.arrow(
                surfaceX, surfaceY, awayX - surfaceX, awayY - surfaceY, color="blue"
            )


def plotCircle(radius: float, firstPoint: Point):
    """
    Plots a circle given a radius and the first point.
    """
    circle = plt.Circle(
        (firstPoint.x + radius, firstPoint.y + radius), radius, fill=False, color="k"
    )
    axes = plt.gca()
    axes.add_patch(circle)


def plotCps(points: list, cps: list, color="blue", label=""):
    """
    Scatter plots c_p vs x/c given the points and cps in order.
    """
    maxX = max(point.x for point in points)
    minX = min(point.x for point in points)
    chordLength = maxX - minX
    first = True
    for point, cp in zip(points, cps):
        if first:            
            plt.plot([(point.x / chordLength)], [cp], "o", color=color, label=label, markersize=3)
            first = False
        else:
            plt.plot([(point.x / chordLength)], [cp], "o", color=color, markersize=3)
    plt.xlabel("$x/c$")
    plt.ylabel("$c_p$")
    plt.gca().invert_yaxis()


def plotCpsForCircle(panels: list, cps: list):
    """
    Scatter plots c_p vs theta given the panels and cps in order.
    """
    for panel, cp in zip(panels, cps):
        plt.plot([math.degrees(panel.delta)], [cp], "o", color="blue")
    plt.xlabel("$\Theta$")
    plt.ylabel("$c_p$")


def plotCircleAnalyticalCps(resolution: int):
    """
    Plots the curve for the solution cp=1-4sin^2(theta) of a given resolution.
    """
    cps = []
    thetas = []
    xs = range(1, resolution)
    step = resolution / (2 * math.pi)
    xs = [theta / step for theta in xs]
    for theta in xs:
        thetas.append(math.degrees(theta))
        cp = 1 - 4 * (math.sin(theta)) ** 2
        cps.append(cp)
    plt.plot(thetas, cps)


def plotAlphaAndCls(alphas: list, cps: list):
    """
    Scatter plots c_l vs alpha given the angles and cps in order.
    """
    plt.scatter(alphas, cps, label="Calculated")
    plt.xlabel(r"$\alpha$")
    plt.ylabel("$c_l$")


def setAxes():
    """
    Sets the aspect ratio of the axes to 1.
    """
    axes = plt.gca()
    axes.set_aspect(1)
