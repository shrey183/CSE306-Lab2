import numpy as np
from sutherLand import *
from parallelVornoi import *
from shapely.geometry import Polygon
import math
from tess import Container


def isClockwise(points):
    # points is your list (or array) of 2d points.
    s = 0.0
    for p1, p2 in zip(points, points[1:] + [points[0]]):
        s += (p2[0] - p1[0]) * (p2[1] + p1[1])
    return s > 0.0


def poly_area(Polygon):
    area = 0.0
    n = len(Polygon)
    if n == 0: return 0
    for i in range(n):
        Xi, Yi = tuple(Polygon[i])
        Xi_, Yi_ = tuple(Polygon[i - 1]) if i > 0 else tuple(Polygon[n - 1])
        area += Xi_ * Yi - Xi * Yi_
    return abs(area)/2

class VornoiCell:
    def __init__(self, P, Polygon, weight=0, lambda_=0):
        self.P = P
        self.Polygon = Polygon
        self.W = weight
        self.L = lambda_

    def insideLaguerre(self, X, Pj, wj):
        X, Pi, Pj = np.array(X), np.array(self.P), np.array(Pj)
        M = (Pi + Pj)/2 + (self.W - wj)*(Pi- Pj)/(2 * np.linalg.norm(Pi - Pj)**2)
        return np.dot(X - M, Pj - Pi) < 0

    def intersectionLaguerre(self, A, B, Pj, wj):
        A, B, Pi, Pj = np.array(A), np.array(B), np.array(self.P), np.array(Pj)
        M = (Pi + Pj)/2 + (self.W - wj)*(Pj- Pi)/(2 * np.linalg.norm(Pi - Pj)**2)
        top = np.dot(M - A, Pi - Pj)
        bottom = np.dot(B - A, Pi - Pj)
        t = top/bottom
        P = A + t*(B - A)
        return P

    def clipLaguerre(self, Pj, wj):
        new_points = []
        subjectPolygon = self.Polygon
        poly_size = len(subjectPolygon)
        for i in range(poly_size):
            k = (i+1) % poly_size
            A, B  = subjectPolygon[i], subjectPolygon[k]
            if(self.insideLaguerre(B, Pj, wj)):
                if (not self.insideLaguerre(A, Pj, wj)):
                    P = self.intersectionLaguerre(A, B, Pj, wj)
                    new_points.append(P)
                new_points.append(B)
            elif self.insideLaguerre(A, Pj, wj):
                P = self.intersectionLaguerre(A, B, Pj, wj)
                new_points.append(P)
        self.Polygon = new_points

    def triangle_area(self, tri):
        x1, y1, x2, y2, x3, y3 = tri[0][0], tri[0][1], tri[1][0], tri[1][1], tri[2][0], tri[2][1]
        return abs(0.5 * (((x2-x1)*(y3-y1))-((x3-x1)*(y2-y1))))

    def poly_integral(self):
        n = len(self.Polygon)
        Pix, Piy = self.P[0], self.P[1]
        norm = Pix**2 + Piy**2
        v1 = 0
        tempPoly = self.Polygon
        for k in range(n):
            Xk, Yk = tuple(tempPoly[k])
            Xk_, Yk_ = tuple(tempPoly[k - 1]) if k > 0 else tuple(tempPoly[n - 1])
            v1 += (Xk_ * Yk - Xk * Yk_)*(Xk_**2 + Xk_ * Xk + Xk**2 + Yk_**2 + Yk_ * Yk + \
                    Yk**2 - 4*(Pix*(Xk_ + Xk) + Piy*(Yk_ + Yk))+ 6*norm)
        return abs(v1)/12


    def poly_intersection_area(self, otherPolygon):
        p = Polygon(self.Polygon).buffer(0)
        q = Polygon(otherPolygon).buffer(0)
        v2 = 0
        if p.intersects(q):
            r = p.intersection(q).buffer(0)
            v2 = r.area
        return abs(v2)

    def poly_area(self):
        area = 0.0
        n = len(self.Polygon)
        if n == 0: return 0
        for i in range(n):
            Xi, Yi = tuple(self.Polygon[i])
            Xi_, Yi_ = tuple(self.Polygon[i - 1]) if i > 0 else tuple(self.Polygon[n - 1])
            area += Xi_ * Yi - Xi * Yi_
        return abs(area)/2


def helper_Laguerre(VornoiCell_list):
    num_points = len(VornoiCell_list)

    for i in range(num_points):
        Vcell = VornoiCell_list[i]
        pi = Vcell.P
        wi = Vcell.W
        for j in range(num_points):
            if i == j: continue
            else:
                temp_Vcell = VornoiCell_list[j]
                pj = temp_Vcell.P
                wj = temp_Vcell.W
                Vcell.clipLaguerre(pj, wj)


def Laguerre(point_list, W, H, weights):
    BoundingBox = box_2d(W, H)
    num_points = len(point_list)
    VornoiCell_list = [VornoiCell(point_list[i], BoundingBox, weights[i]) for i in range(num_points)]
    helper_Laguerre(VornoiCell_list)
    outPolygon = {}
    for i in range(num_points):
        Vcell = VornoiCell_list[i]
        outPolygon[tuple(Vcell.P)] = Vcell.Polygon
    return outPolygon
