from panelGeometry import Point, Panel
import math
import csv
import matplotlib
import numpy as np


# Methods


def computeLambdas(panels: list, freestreamVelocity: float) -> list:
    for index, panel in enumerate(panels):
        if panel.startPoint.x - panels[index - 1].endPoint.x > 0.00000001:
            raise Exception("Panels must form a closed path.")
        if panel.startPoint.y - panels[index - 1].endPoint.y > 0.00000001:
            raise Exception("Panels must form a closed path.")
    matrixA = []
    matrixB = []
    for paneli in panels:
        lambdas = []
        for panelj in panels:
            if paneli == panelj:
                lambdas.append(1 / 2)
            else:
                a = -(paneli.controlPoint.x - panelj.startPoint.x) * math.cos(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.sin(panelj.phi)
                b = (paneli.controlPoint.x - panelj.startPoint.x) ** 2 + (
                    paneli.controlPoint.y - panelj.startPoint.y
                ) ** 2
                c = math.sin(paneli.phi - panelj.phi)
                d = (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(
                    paneli.phi
                ) - (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(paneli.phi)
                Sj = math.sqrt(
                    (panelj.endPoint.x - panelj.startPoint.x) ** 2
                    + (panelj.endPoint.y - panelj.startPoint.y) ** 2
                )
                e = (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(panelj.phi)
                integral = (c / 2) * math.log((Sj ** 2 + 2 * a * Sj + b) / b) + (
                    (d - a * c) / e
                ) * (math.atan((Sj + a) / e) - math.atan(a / e))
                lambdas.append(integral / (2 * math.pi))
        matrixA.append(lambdas)
        matrixB.append(-freestreamVelocity * math.cos(paneli.beta))
    lambdas = np.linalg.solve(matrixA, matrixB)
    return lambdas


def computeGammas(panels: list, freestreamVelocity: float) -> list:
    for index, panel in enumerate(panels):
        if panel.startPoint.x - panels[index - 1].endPoint.x > 0.00000001:
            raise Exception("Panels must form a closed path.")
        if panel.startPoint.y - panels[index - 1].endPoint.y > 0.00000001:
            raise Exception("Panels must form a closed path.")
    matrixA = []
    matrixB = []
    for paneli in panels:
        gamas = []
        for panelj in panels:
            if paneli == panelj:
                gamas.append(1)
            else:
                a = -(paneli.controlPoint.x - panelj.startPoint.x) * math.cos(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.sin(panelj.phi)
                b = (paneli.controlPoint.x - panelj.startPoint.x) ** 2 + (
                    paneli.controlPoint.y - panelj.startPoint.y
                ) ** 2
                c = math.sin(paneli.phi - panelj.phi)
                d = (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(
                    paneli.phi
                ) - (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(paneli.phi)
                Sj = math.sqrt(
                    (panelj.endPoint.x - panelj.startPoint.x) ** 2
                    + (panelj.endPoint.y - panelj.startPoint.y) ** 2
                )
                e = (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(panelj.phi)
                integral = (c / 2) * math.log((Sj ** 2 + 2 * a * Sj + b) / b) + (
                    (d - a * c) / e
                ) * (math.atan((Sj + a) / e) - math.atan(a / e))
                gamas.append(integral / (2 * math.pi))
        matrixA.append(gamas)
        matrixB.append(freestreamVelocity * math.cos(paneli.beta))
    kuttaCond = list(np.zeros(len(gamas)))
    kuttaCond[0] = 1
    kuttaCond[-1] = 1
    matrixA[-1] = kuttaCond
    matrixB[-1] = 0
    gamas = np.linalg.solve(matrixA, matrixB)
    return gamas


def computeCpsFromLambdas(
    panels: list, lambdas: list, freestreamVelocity: float
) -> tuple:
    accuracy = 0
    cps = []
    for paneli in panels:
        vi = 0
        for panelj, lamj in zip(panels, lambdas):
            if paneli == panelj:
                vi += freestreamVelocity * math.sin(paneli.beta)
                accuracy += lamj * panelj.length
            else:
                a = -(paneli.controlPoint.x - panelj.startPoint.x) * math.cos(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.sin(panelj.phi)
                b = (paneli.controlPoint.x - panelj.startPoint.x) ** 2 + (
                    paneli.controlPoint.y - panelj.startPoint.y
                ) ** 2
                c = math.sin(paneli.phi - panelj.phi)
                d = (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(
                    paneli.phi
                ) - (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(paneli.phi)
                Sj = math.sqrt(
                    (panelj.endPoint.x - panelj.startPoint.x) ** 2
                    + (panelj.endPoint.y - panelj.startPoint.y) ** 2
                )
                e = (paneli.controlPoint.x - panelj.startPoint.x) * math.sin(
                    panelj.phi
                ) - (paneli.controlPoint.y - panelj.startPoint.y) * math.cos(panelj.phi)
                integral = ((d - a * c) / (2 * e)) * math.log(
                    (Sj ** 2 + 2 * a * Sj + b) / b
                ) - c * (math.atan((Sj + a) / e) - math.atan(a / e))
                vi += (lamj / (2 * math.pi)) * integral
        cps.append(1 - (vi / freestreamVelocity) ** 2)
    print(f"Accuracy: {accuracy}")
    return cps, accuracy


def computeCpsFromGammas(panels: list, gammas: list, freestreamVelocity: float) -> list:
    cps = []
    for panelj, gamj in zip(panels, gammas):
        dx = panelj.startPoint.x - panelj.endPoint.x
        dy = panelj.startPoint.y - panelj.endPoint.y
        cps.append(1 - (gamj / freestreamVelocity) ** 2)
    return cps
