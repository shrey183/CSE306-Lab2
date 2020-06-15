from sutherLand import  *
from parallelVornoi import *
from PowerDiagram import *
from SemiDiscreteOptimal import *


def run_SutherLand():
    # Points defined in anti-clockwise manner, Example taken from Rosetta code
    subjectPolygon = [[50, 150], [200, 50], [350, 150], [350, 300], [250, 300],
                      [200, 250], [150, 350], [100, 250], [100, 200]]
    clipPolygon = [[100, 100], [300, 100], [300, 300], [100, 300]]

    # algorithm needs points defined in clockwise manner
    result = poly_clip(subjectPolygon[::-1], clipPolygon[::-1])

    polygon_list = [subjectPolygon, clipPolygon, result]
    save_svg(polygon_list, "../Images/Example_SutherLand.svg")

def run_ParallelVoronoi():
    W = 1
    H = 1
    point_list = [[W/4, H/4],[3*W/4, H/4],[3*W/4, 3*H/4], [W/4, 3*H/4]]
    # algorithm needs points defined in clockwise manner
    polygon_list = voronoi(point_list, W, H)
    save_svg_vornoi(polygon_list, "../Images/Example_Voronoi.svg")

def run_PowerDiagram():
    # Then plot using Vornoi Diagram
    W = 1
    H = 1
    weights = [0.1, 0.2, 0.4, 0.3]
    point_list = [[W/4, H/4],[3*W/4, H/4],[3*W/4, 3*H/4], [W/4, 3*H/4]]
    num_points = len(point_list)
    polygon_list = Laguerre(point_list, W, H, weights)
    VornoiCell_list = []
    for i in range(num_points):
        temp = VornoiCell(point_list[i], polygon_list[tuple(point_list[i])], weights[i])
        VornoiCell_list.append(temp)

    area_list = [VornoiCell_list[i].poly_area() for i in range(num_points)]
    print("Area of Polygon: ", area_list)
    save_svg_vornoi(polygon_list, "../Images/Example_PowerVornoi.svg")


def run_SemiDiscreteOptimalTransport():
    # Simple Example
    W = 1
    H = 1
    A = W*H
    point_list = [[W/4, H/4],[3*W/4, H/4],[3*W/4, 3*H/4], [W/4, 3*H/4]]
    lambda_list = [0.125, 0.5, 0.125, 0.25]
    new_weights = SDOptimalTransport(point_list, lambda_list, np.zeros((4)))
    polygon_list = Laguerre(point_list, W, H, new_weights)
    save_svg_vornoi(polygon_list, "../Images/Example_SDO.svg")


if __name__ == '__main__':
    run_SutherLand()
    run_ParallelVoronoi()
    run_PowerDiagram()
    run_SemiDiscreteOptimalTransport()
