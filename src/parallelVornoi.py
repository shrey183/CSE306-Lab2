# write to svg file
import numpy as np

color_list = [('darkred', 'rgb(139, 0, 0)'),
        ('darksalmon', 'rgb(233, 150, 122)'),
        ('orange', 'rgb(255, 165, 0)'),
        ('orangered', 'rgb(255, 69, 0)'),
        ('orchid', 'rgb(218, 112, 214)'),
        ('palegoldenrod', 'rgb(238, 232, 170)'),
        ('palegreen', 'rgb(152, 251, 152)'),
        ('paleturquoise', 'rgb(175, 238, 238)'),
        ('palevioletred', 'rgb(219, 112, 147)'),
]


def save_svg_vornoi(polygon_list, filename, fillcol = "none"):
    with open(filename, "w") as f:
        f.write("<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n")
        # for each polygon
        i = 0
        for p, polygon in polygon_list.items():
            f.write("<g>\n")
            f.write("<polygon points = \"")
            # for each vertex in the polygon
            for j in range(len(polygon)):
                f.write("%3.3f, %3.3f " % (polygon[j][0]*1000, 1000- (polygon[j][1]*1000)))
            f.write("\"\nfill = \"{}\" stroke = \"black\"/>\n".format(color_list[i][0]))
            f.write("</g>\n")
            # for each point in the polygon
            f.write("<g>\n")
            f.write("<circle cx = \"%3.3f\" cy = \"%3.3f\" r=\"2\"/>" % (p[0]*1000, 1000 - p[1]*1000))
            f.write("</g>\n")
            i += 1

        f.write("</svg>\n")


def inside(X, Pi, Pj):
    X, Pi, Pj = np.array(X),np.array(Pi), np.array(Pj)
    M = (Pi + Pj)/2
    return np.dot(X - M, Pj - Pi) < 0

def intersection(A, B, Pi, Pj):
    A, B, Pi, Pj = np.array(A), np.array(B), np.array(Pi), np.array(Pj)
    M = (Pi + Pj)/2
    top = np.dot(M - A, Pi - Pj)
    bottom = np.dot(B - A, Pi - Pj)
    t = top/bottom
    P = A + t*(B - A)
    return P

def clip(subjectPolygon, Pi, Pj):
    new_points = []
    poly_size = len(subjectPolygon)
    for i in range(poly_size):
        k = (i+1) % poly_size
        A, B  = subjectPolygon[i], subjectPolygon[k]
        if(inside(B, Pi, Pj)):
            if (not inside(A, Pi, Pj)):
                P = intersection(A, B, Pi, Pj)
                new_points.append(P)
            new_points.append(B)
        elif inside(A, Pi, Pj):
            P = intersection(A, B, Pi, Pj)
            new_points.append(P)

    return new_points

def voronoi(point_list, W, H, d=2):
    outPolygon = {}
    n = len(point_list)
    bbox = box_2d(W, H)
    if d == 3:
        bbox = box_3d(W, H, max(point_list, key=lambda p:p[2])[2])
    for i in range(n):
        pi = point_list[i]
        tempPoly = []
        voronoiCell = bbox
        for j in range(n):
            if i == j: continue
            else:
                pj = point_list[j]
                tempPoly = clip(voronoiCell, pi, pj)
                voronoiCell = tempPoly
        outPolygon[tuple(pi)] = voronoiCell

    return outPolygon

def box_2d(W, H):
    return np.array([[0,0], [0, H], [W, H], [W, 0]]) # clockwise

def box_3d(W, H, D):
    return np.array([[0, 0, 0], [0, H, 0], [W, H, 0], [W, 0, 0],
                     [0, 0, D+2], [0, H, D+2], [W, H, D+2], [W, 0, D+2]])
