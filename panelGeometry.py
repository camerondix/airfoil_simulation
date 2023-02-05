import math
import csv


# Classes
class Point:
    """
    A 2D point containing an x and y value.
    """

    def __init__(this, x: float, y: float):
        this.x = float(x)
        this.y = float(y)


class Panel:
    """
    A 2D panel composed of two points.
    """

    def __init__(this, point1: Point, point2: Point, alpha=0):
        if point1.x == point2.x and point1.y == point2.y:
            raise Exception("A line cannot be made from two identical points.")
        this.startPoint = point1
        this.endPoint = point2
        this.controlPoint = Point((point1.x + point2.x) / 2, (point1.y + point2.y) / 2)
        this.dx = point2.x - point1.x
        this.dy = point2.y - point1.y
        this.length = math.sqrt((point2.x - point1.x) ** 2 + (point2.y - point1.y) ** 2)
        phi = math.atan2(this.dy, this.dx)
        if phi < 0.0:
            phi += 2 * math.pi
        this.phi = phi
        delta = this.phi + math.pi / 2
        while delta >= 2 * math.pi:
            delta -= 2 * math.pi
        this.delta = delta
        this.beta = this.delta - alpha


# Methods
def importPoints(fileName: str, seperator: str) -> list:
    """
    Imports a list of points from a file.
    """
    lines = []
    with open(fileName, "r") as file:
        reader = csv.reader(
            file,
            skipinitialspace=True,
            delimiter=seperator,
        )
        for line in reader:
            if len(line) != 2:
                raise Exception("Data must contain only two columns.")
            lines.append(line)
    points = []
    for line in lines:
        points.append(Point(line[0], line[1]))
    if not checkIfPointsAreCW(points):
        points.reverse()
    return points


def createPointsFromArrays(xs: list, ys: list) -> list:
    """
    Creates a list of points from two ordered lists of x and y values. The default direction is clockwise.
    """
    if len(xs) != len(ys):
        raise Exception("X and Y coordinates must be the same length.")
    points = []
    for x, y in zip(xs, ys):
        points.append(Point(x, y))
    if not checkIfPointsAreCW(points):
        points.reverse()
    return points


def createPanelsFromPoints(points: list, alpha=0) -> list:
    """
    Creates a list of panels from an ordered list of points at the angle of attack alpha.
    """
    panels = []
    for index, point in enumerate(points):
        panel = Panel(points[index - 1], point, alpha)
        panels.append(panel)
    # Remove the TE panel if it is vertical
    if panels[0].dx == 0.0:
        panels.remove(panels[0])
    if panels[-1].dx == 0.0:
        panels.remove(panels[-1])
    return panels


def createCirclePanels(radius: float, divisions: int, alpha=0) -> list:
    """
    Creates a list of panels from a circle at the angle of attack alpha.
    """
    panels = []
    step = (2 * math.pi) / divisions
    initialTheta = 0 - step / 2
    theta = initialTheta
    count = 0
    previousPoint = Point(0, 0)
    while count <= divisions:
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        point = Point(x, y)
        if count > 0:
            panel = Panel(point, previousPoint, alpha)
            panels.append(panel)
        previousPoint = point
        theta += step
        count += 1
    return panels


def importTuples(fileName: str, seperator: str) -> tuple:
    lines = []
    with open(fileName, "r", encoding='utf-8-sig') as file:
        reader = csv.reader(
            file,
            skipinitialspace=True,
            delimiter=seperator,
        )
        for line in reader:
            if len(line) != 2:
                raise Exception("Data must contain only two columns.")
            lines.append(line)
    xs = []
    tuples = []
    for line in lines:
        xs.append(float(line[0]))
        tuples.append(float(line[1]))
    return xs, tuples


def checkIfPointsAreCW(points: list) -> bool:
    summation = 0
    for i in range(len(points)):
        start = points[i - 1]
        end = points[i]
        summation += (end.x - start.x) * (end.y + start.y)
    if summation < 0:
        return False
    else:
        return True


def findIij(paneli: Panel, panelj: Panel) -> float:
    """
    Normal velocity geometric integral of panel i relative to panel j.
    """
    x_i = paneli.controlPoint.x
    y_i = paneli.controlPoint.y
    phi_i = paneli.phi
    x_j = panelj.startPoint.x
    y_j = panelj.startPoint.y
    phi_j = panelj.phi
    s_j = panelj.length
    a = -(x_i - x_j) * math.cos(phi_j) - (y_i - y_j) * math.sin(phi_j)
    b = (x_i - x_j) ** 2 + (y_i - y_j) ** 2
    c = math.sin(phi_i - phi_j)
    d = -(x_i - x_j) * math.sin(phi_i) + (y_i - y_j) * math.cos(phi_i)
    e = (b - a ** 2) ** 0.5
    i_ij = (c / 2) * math.log((s_j ** 2 + 2 * a * s_j + b) / b) + ((d - a * c) / e) * (
        math.atan((s_j + a) / e) - math.atan(a / e)
    )
    return i_ij


def findJij(paneli: Panel, panelj: Panel) -> float:
    """
    Tangential velocity geometric integral of panel i relative to panel j.
    """
    x_i = paneli.controlPoint.x
    y_i = paneli.controlPoint.y
    phi_i = paneli.phi
    x_j = panelj.startPoint.x
    y_j = panelj.startPoint.y
    phi_j = panelj.phi
    s_j = panelj.length
    a = -(x_i - x_j) * math.cos(phi_j) - (y_i - y_j) * math.sin(phi_j)
    b = (x_i - x_j) ** 2 + (y_i - y_j) ** 2
    c = -math.cos(phi_i - phi_j)
    d = (x_i - x_j) * math.cos(phi_i) + (y_i - y_j) * math.sin(phi_i)
    e = (b - a ** 2) ** 0.5
    j_ij = (c / 2) * math.log((s_j ** 2 + 2 * a * s_j + b) / b) + ((d - a * c) / e) * (
        math.atan((s_j + a) / e) - math.atan(a / e)
    )
    return j_ij


def findLij(paneli: Panel, panelj: Panel) -> float:
    """
    Tangential velocity geometric integral of panel i relative to panel j.
    """
    x_i = paneli.controlPoint.x
    y_i = paneli.controlPoint.y
    phi_i = paneli.phi
    x_j = panelj.startPoint.x
    y_j = panelj.startPoint.y
    phi_j = panelj.phi
    s_j = panelj.length
    a = -(x_i - x_j) * math.cos(phi_j) - (y_i - y_j) * math.sin(phi_j)
    b = (x_i - x_j) ** 2 + (y_i - y_j) ** 2
    c = math.sin(phi_j - phi_i)
    d = (x_i - x_j) * math.sin(phi_i) - (y_i - y_j) * math.cos(phi_i)
    e = (b - a ** 2) ** 0.5
    l_ij = (c / 2) * math.log((s_j ** 2 + 2 * a * s_j + b) / b) + ((d - a * c) / e) * (
        math.atan((s_j + a) / e) - math.atan(a / e)
    )
    return l_ij
