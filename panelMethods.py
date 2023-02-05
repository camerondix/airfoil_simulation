import math
import panelGeometry as pg
import numpy as np

def findSourcePanelStrengths(panels: list, freestreamVelocity: float) -> list:
    """
    Finds the source panel strengths using the source panel method.
    """
    matrixA = []
    matrixB = []
    for i, paneli in enumerate(panels):
        rowA = []
        for j, panelj in enumerate(panels):
            if i is j:
                rowA.append(math.pi)
            else:
                i_ij = pg.findIij(paneli, panelj)
                rowA.append(i_ij)
        matrixA.append(rowA)
        rowB = -freestreamVelocity * 2 * math.pi * math.cos(paneli.beta)
        matrixB.append(rowB)
    lambdas = np.linalg.solve(matrixA, matrixB)
    return lambdas


def findSourcePanelCoefficients(
    panels: list, freestreamVelocity: float, alpha: float
) -> tuple:
    """
    Finds the pressure coefficient at each panel and the total lift and drag coefficients using the source panel method.
    """
    cps = []
    clFirstTerm = 0
    clSecondTerm = 0
    cdFirstTerm = 0
    cdSecondTerm = 0
    accuracy = 0
    lambdas = findSourcePanelStrengths(panels, freestreamVelocity)
    for i, paneli in enumerate(panels):
        v = 0
        accuracy += paneli.length * lambdas[i]
        for j, panelj in enumerate(panels):
            if i is j:
                v += freestreamVelocity * math.sin(paneli.beta)
            else:
                j_ij = pg.findJij(paneli, panelj)
                v += (lambdas[j] / (2 * math.pi)) * j_ij
        cp = 1 - (v / freestreamVelocity) ** 2
        cps.append(cp)
        cn = -cp * paneli.length * math.sin(paneli.beta)
        ca = -cp * paneli.length * math.cos(paneli.beta)
        clFirstTerm += cn * math.cos(alpha)
        clSecondTerm += ca * math.sin(alpha)
        cdFirstTerm += cn * math.sin(alpha)
        cdSecondTerm += ca * math.cos(alpha)
    cl = clFirstTerm - clSecondTerm
    cd = cdFirstTerm + cdSecondTerm
    return cps, cl, cd, accuracy


def findVortexPanelStrengths(panels: list, freestreamVelocity: float) -> list:
    """
    Finds the vortex panel strengths using the vortex panel method.
    """
    matrixA = []
    matrixB = []
    rowA = []
    for i, paneli in enumerate(panels):
        rowA = []
        for j, panelj in enumerate(panels):
            if i is j:
                rowA.append(0)
            else:
                j_ij = pg.findJij(paneli, panelj)
                rowA.append(-j_ij)
        matrixA.append(rowA)
        rowB = -freestreamVelocity * 2 * math.pi * math.cos(paneli.beta)
        matrixB.append(rowB)
    # Apply the Kutta condition
    kuttaConditionA = [row * 0 for row in rowA]
    kuttaConditionA[0] = 1
    kuttaConditionA[-1] = 1
    kuttaConditionB = 0
    matrixA[-1] = kuttaConditionA
    matrixB[-1] = kuttaConditionB
    gammas = np.linalg.solve(matrixA, matrixB)
    return gammas


def findVortexPanelCoefficients(
    panels: list, freestreamVelocity: float, alpha: float
) -> tuple:
    """
    Finds the pressure coefficient at each panel and the total lift, drag, and moment coefficients using the vortex panel method.
    """
    cps = []
    clFirstTerm = 0
    clSecondTerm = 0
    cdFirstTerm = 0
    cdSecondTerm = 0
    accuracy = 0
    cm = 0
    gammas = findVortexPanelStrengths(panels, freestreamVelocity)
    for i, paneli in enumerate(panels):
        v = 0
        accuracy += paneli.length * gammas[i]
        for j, panelj in enumerate(panels):
            if i is j:
                v += freestreamVelocity * math.sin(paneli.beta) + gammas[i] / 2
            else:
                l_ij = pg.findLij(paneli, panelj)
                v -= (gammas[j] / (2 * math.pi)) * l_ij
        cp = 1 - (v / freestreamVelocity) ** 2
        cps.append(cp)
        cn = -cp * paneli.length * math.sin(paneli.beta)
        ca = -cp * paneli.length * math.cos(paneli.beta)
        clFirstTerm += cn * math.cos(alpha)
        clSecondTerm += ca * math.sin(alpha)
        cdFirstTerm += cn * math.sin(alpha)
        cdSecondTerm += ca * math.cos(alpha)
        cm += cp * (paneli.controlPoint.x - 0.25) * paneli.length * math.cos(paneli.phi)
    cl = clFirstTerm - clSecondTerm
    cd = cdFirstTerm + cdSecondTerm
    return cps, cl, cd, cm, accuracy


def findSourceVortexPanelStrengths(panels: list, freestreamVelocity: float) -> list:
    """
    Finds the source and vortex panel strengths using a source/vortex panel method.
    """
    matrixA = []
    matrixB = []
    for i, paneli in enumerate(panels):
        rowA = []
        sumKs = 0
        for j, panelj in enumerate(panels):
            if i is j:
                rowA.append(math.pi)
            else:
                i_ij = pg.findIij(paneli, panelj)
                rowA.append(i_ij)
                j_ij = pg.findJij(paneli, panelj)
                sumKs += j_ij
        rowA.append(-sumKs)
        matrixA.append(rowA)
        rowB = -freestreamVelocity * 2 * math.pi * math.cos(paneli.beta)
        matrixB.append(rowB)
    # Apply the Kutta condition
    rowA = []
    sumLs = 0
    for j, panelj in enumerate(panels):
        summation = 0
        if j != 0:
            j_1j = pg.findJij(panels[0], panels[j])
            summation += j_1j
            l_1j = pg.findLij(panels[0], panels[j])
            sumLs += l_1j
        if j != len(panels) - 1:
            j_Nj = pg.findJij(panels[-1], panels[j])
            summation += j_Nj
            l_Nj = pg.findLij(panels[-1], panels[j])
            sumLs += l_Nj
        rowA.append(summation)
    rowA.append(-sumLs + 2 * math.pi)
    matrixA.append(rowA)
    rowB = (
        -freestreamVelocity
        * 2
        * math.pi
        * (math.sin(panels[0].beta) + math.sin(panels[-1].beta))
    )
    matrixB.append(rowB)
    lambdasAndGamma = np.linalg.solve(matrixA, matrixB)
    return lambdasAndGamma


def findSourceVortexPanelCoefficients(
    panels: list, freestreamVelocity: float, alpha: float
) -> tuple:
    """
    Finds the pressure coefficient at each panel and the total lift, drag, and moment coefficients using a source/vortex panel method.
    """
    cps = []
    clFirstTerm = 0
    clSecondTerm = 0
    cdFirstTerm = 0
    cdSecondTerm = 0
    cm = 0
    lambdasAndGamma = findSourceVortexPanelStrengths(panels, freestreamVelocity)
    for i, paneli in enumerate(panels):
        sumJ = 0
        sumL = 0
        gamma = lambdasAndGamma[-1]
        for j, panelj in enumerate(panels):
            if i != j:
                lambd = lambdasAndGamma[j]
                J_ij = pg.findJij(paneli, panelj)
                sumJ += lambd * J_ij
                L_ij = pg.findLij(paneli, panelj)
                sumL += L_ij
        v = (
            freestreamVelocity * math.sin(paneli.beta)
            + (1 / (2 * math.pi)) * sumJ
            + gamma / 2
            - (gamma / (2 * math.pi)) * sumL
        )
        cp = 1 - (v / freestreamVelocity) ** 2
        cps.append(cp)
        cn = -cp * paneli.length * math.sin(paneli.beta)
        ca = -cp * paneli.length * math.cos(paneli.beta)
        clFirstTerm += cn * math.cos(alpha)
        clSecondTerm += ca * math.sin(alpha)
        cdFirstTerm += cn * math.sin(alpha)
        cdSecondTerm += ca * math.cos(alpha)
        cm += cp * (paneli.controlPoint.x - 0.25) * paneli.length * math.cos(paneli.phi)
    cl = clFirstTerm - clSecondTerm
    cd = cdFirstTerm + cdSecondTerm
    return cps, cl, cd, cm
