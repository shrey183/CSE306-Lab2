import numpy as np
from parallelVornoi import *
from PowerDiagram import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt

"""
The goal is compute the optimum weights
corresponding to the given mass lambda_i

Input:
        Vornoi Sites: y_i
        Mass Value: lambda_i
"""

def g(VornoiCell_list):
    n = len(VornoiCell_list)
    area_term = 0
    for i in range(n):
        Vcell = VornoiCell_list[i]
        area_term += Vcell.poly_integral() - Vcell.W * Vcell.poly_area() \
                     + Vcell.L * Vcell.W
    return area_term

def grad_g(VornoiCell_list):
    grad_vector = []
    n = len(VornoiCell_list)
    for i in range(n):
        Vcell = VornoiCell_list[i]
        grad_i = - Vcell.poly_area() + Vcell.L
        grad_vector.append(grad_i)
    return np.array(grad_vector)

def invHessian(VornoiCell_list):
    n = len(VornoiCell_list)
    H = np.zeros((n,n))

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            Vi, Vj  = VornoiCell_list[i], VornoiCell_list[j]
            Pi, Pj = np.array(Vi.P), np.array(Vj.P)
            denom = np.linalg.norm(Pj - Pi)
            area = Vi.poly_intersection_area(Vj.Polygon)
            H[i, j] = (1/(2*denom))*area

    for i in range(n):
        H[i, i] = - np.sum(H[i])
    invH = np.linalg.pinv(H)
    return invH

def emptyCell(VornoiCell_list):
    num_points = len(VornoiCell_list)
    for i in range(num_points):
        VCell = VornoiCell_list[i].Polygon
        if len(VCell) == 0:
            return True
    return False


def SDOptimalTransport(point_list, lambda_list, weights_init, num_iter=1000, eps = 0.01):
    W = 1
    H = 1
    num_points = len(point_list)
    BoundingBox = box_2d(W, H)

    VornoiCell_list = [VornoiCell(point_list[i], BoundingBox, weights_init[i], lambda_list[i]) for i in range(num_points)]
    g_values = []
    new_weights = []
    for i in range(num_iter):
        print("Iteration: {}\r".format(i+1), end="")
        new_weights.clear()

        # gradient desecent for each weight_i
        g_value = g(VornoiCell_list)
        g_values.append(g_value)

        grad_vector = grad_g(VornoiCell_list)
        invH = invHessian(VornoiCell_list)
        dir = 2*np.matmul(invH, grad_vector)

        if emptyCell(VornoiCell_list): eps = eps/2

        for j in range(num_points):
            VornoiCell_list[j].W += eps*dir[j]
            new_weights.append(VornoiCell_list[j].W)


        # Update the polygons
        helper_Laguerre(VornoiCell_list)

    area_list = [VornoiCell_list[i].poly_area() for i in range(num_points)]

    print("Area of Polygon: ", area_list)
    print("Lambda List: ", lambda_list)
    print("Weight List: ", new_weights)
    print("Max g value: ", max(g_values))

    plt.plot(np.arange(num_iter), np.array(g_values))
    plt.xlabel("Number of Iterations")
    plt.ylabel("g(W)")
    plt.title("g(W) vs Number of Iterations")
    plt.show()

    return new_weights
