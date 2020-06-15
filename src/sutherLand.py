# write to svg file
import numpy as np
from copy import deepcopy

def save_svg(polygon_list, filename, fillcol = "none"):
    with open(filename, "w") as f:
        f.write("<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n")
        # for each polygon
        color_list = ["red", "blue", "yellow"]
        for polygon in polygon_list:
            f.write("<g>\n")
            f.write("<polygon points = \"")
            # for each vertex in the polygon
            for j in range(len(polygon)):
                f.write("%3.3f, %3.3f " % (polygon[j][0], 1000-polygon[j][1]))
            f.write("\"\nfill = \"{}\" stroke = \"black\"/>\n".format(color_list.pop(0)))
            f.write("</g>\n")
        f.write("</svg>\n")


def inside(u, v, P):
    x1, y1 = u[0], u[1]
    x2, y2 = v[0], v[1]
    x, y = P[0], P[1]
    return (x2-x1) * (y-y1) - (y2-y1) * (x-x1) < 0

def intersection(u, v, A, B):
    x1, y1 = u[0], u[1]
    x2, y2 = v[0], v[1]
    x3, y3 = A[0], A[1]
    x4, y4 = B[0], B[1]
    x_num = (x1*y2 - y1*x2) * (x3-x4) - (x1-x2) * (x3*y4 - y3*x4)
    xy_den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4)
    y_num = (x1*y2 - y1*x2) * (y3-y4) - (y1-y2) * (x3*y4 - y3*x4)
    return [x_num/xy_den, y_num/xy_den]

def clip(subjectPolygon, u, v):
    x1, y1 = u[0], u[1]
    x2, y2 = v[0], v[1]
    new_points = []
    poly_size = len(subjectPolygon)
    for i in range(poly_size):
        k = (i+1) % poly_size
        A, B  = subjectPolygon[i], subjectPolygon[k]
        if(inside(u, v, B)):
            if (not inside(u, v, A)):
                tx, ty = intersection(u, v, A, B)
                new_points.append([tx, ty])
            new_points.append([B[0], B[1]])
        elif inside(u, v, A):
            tx, ty = intersection(u, v, A, B)
            new_points.append([tx, ty])

    return new_points


def poly_clip(subjectPolygon, clipPolygon):
    outPolygon = None
    clipper_size = len(clipPolygon)
    for i in range(clipper_size):
        k = (i+1) % clipper_size
        u = clipPolygon[i]
        v = clipPolygon[k]
        # clip the subject polygon by the edge (u, v)
        outPolygon = clip(subjectPolygon, u, v)
        subjectPolygon = outPolygon

    return outPolygon
