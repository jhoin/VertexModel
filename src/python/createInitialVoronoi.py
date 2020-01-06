import numpy as np
import cv2
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString, MultiPolygon, MultiPoint, Point
import shapely.ops
from math import sqrt
from descartes.patch import PolygonPatch
import seaborn as sns

# Create a vornoi diagram and use it as the initial conditions for a vertex Model

# Replacing the fiugures module
BLUE = '#6699cc'
GRAY = '#999999'
RED = '#ff3333'

GM = (sqrt(5)-1.0)/2.0
W = 8.0
H = W*GM
SIZE = (W, H)

COLOR_ISVALID = {
    True: BLUE,
    False: RED,
}

def plot_coords(ax, ob, color=GRAY, zorder=1, alpha=1):
    x, y = ob.xy
    ax.plot(x, y, '.', color=color, zorder=zorder, alpha=alpha)

def color_isvalid(ob, valid=BLUE, invalid=RED):
    if ob.is_valid:
        return valid
    else:
        return invalid

def set_limits(ax, x0, xN, y0, yN):
    ax.set_xlim(x0, xN)
    ax.set_xticks(range(x0, xN+1))
    ax.set_ylim(y0, yN)
    ax.set_yticks(range(y0, yN+1))
    ax.set_aspect("equal")

def plot_grid(plot_pos, grid):

    # In order to plot a MultiPolygon object, I need to iterate over each oplygon
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(plot_pos)

    for element in grid:
            polygon = Polygon(element)
            plot_coords(ax, polygon.exterior)
            patch = PolygonPatch(polygon, facecolor=color_isvalid(polygon), edgecolor=color_isvalid(polygon, valid=BLUE), alpha=0.5, zorder=2)
            ax.add_patch(patch)


########################################################################################################################
# write the initial conditions file directly from the voronoi diagram
# Arguments: the voronoi diagram object created by the scipy function
# Rreturn: the input files for the vertex model
########################################################################################################################
def createInitialFromVoronoi(voronoiDiagram):
    np.savetxt("icPoints.txt", vorDiagram.vertices, header=str(len(vorDiagram.vertices)), comments='')
    np.savetxt("icEdges.txt", vorDiagram.ridge_vertices, header=str(len(vorDiagram.ridge_vertices)), fmt="%d", comments='')

    # Use the pad function
    with open("icCells.txt","w") as f:
        biggestCell = max([len(region) for region in vorDiagram.regions])
        f.write(str(len(vorDiagram.regions)-2)+" "+str(biggestCell)+"\n")
        for region in vorDiagram.regions:
            if not region:
                continue
            elif biggestCell==len(region):
                fileLine = " ".join(map(str, region))
                f.write(fileLine+"\n")
            else:
                fileLine = " ".join(map(str, np.pad(region,(0,biggestCell-len(region)),'constant',constant_values= -999)))
                f.write(fileLine+"\n")

def createInitialFromPoly(listPolygons):

    # Create a list of points
    listPoints = []
    for poly in listPolygons:
        coordsPoly = poly.exterior.coords
        for coord in coordsPoly:
            if coord not in listPoints:
                listPoints.append(coord)

    # Create the list of cells
    listCells = []
    for poly in listPolygons:
        poly = shapely.geometry.polygon.orient(poly,1)
        coordsPoly = poly.exterior.coords
        coordsIndex = [listPoints.index(coord)+1 for coord in coordsPoly]
        listCells.append(coordsIndex[:-1])

    # Write initial condition to file according to the new system
    with open("ic.txt", "w") as pointFile:
        pointFile.write(str(len(listPoints))+"\n")
        for point in listPoints:
            fileLine = " ".join(map(str, point))
            pointFile.write(fileLine+"\n")
        biggestCell = max([len(cell) for cell in listCells])
        pointFile.write(str(len(listCells))+" "+str(biggestCell)+"\n")
        for cell in listCells:
            if not cell:
                continue
            else:
                fileLine = " ".join(map(str, cell))
                pointFile.write(fileLine+"\n")

########################################################################################################################
# write the initial conditions file after cleaning the voronoi diagram
# Arguments: the voronoi diagram object created by the scipy function
# Rreturn: the input files for the vertex model
########################################################################################################################
def createInitialFromPolygons(listPolygons):

    # Write to the point file
    listPoints = []
    for poly in listPolygons:
        coordsPoly = poly.exterior.coords
        for coord in coordsPoly:
            if coord not in listPoints:
                listPoints.append(coord)

    with open("icPoints.txt", "w") as pointFile:
        pointFile.write(str(len(listPoints))+"\n")
        for point in listPoints:
            fileLine = " ".join(map(str, point))
            pointFile.write(fileLine+"\n")


    with open("icCells.txt", "w") as cellsFile:
        biggestCell = max([len(cell) for cell in listCells])
        cellsFile.write(str(len(listCells))+" "+str(biggestCell)+"\n")
        for cell in listCells:
            if not cell:
                continue
            elif biggestCell==len(cell):
                fileLine = " ".join(map(str, cell))
                cellsFile.write(fileLine+"\n")
            else:
                fileLine = " ".join(map(str, np.pad(cell,(0,biggestCell-len(cell)),'constant',constant_values= -999)))
                cellsFile.write(fileLine+"\n")

########################################################################################################################
# Infer the cell geometries from an input microscopy image
# Arguments: a microscopy image staining the cell menbranes
# Rreturn: a shapely list of polygon
########################################################################################################################
def getPolyFromImage():
    im = cv2.imread('epi3.png')

    # Pre process the image with some filters, in order to minimise the details
    imgray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    ret, thresh = cv2.threshold(imgray, 127, 255, 0)

    listContours, h = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

    # Get the centroids for each contour
    centroid_x = []
    centroid_y = []
    area_list = []
    cellObjs = []
    for eachContour in listContours:

        # Find the moments of the contours, whatever this means
        moment_list = cv2.moments(eachContour)

        # By using the moments, calculate the centroids and the area of each contour
        # Open polygons are open and then return area zero, I exclude these to avoid dividing by 0
        if(moment_list['m00'] != 0.0):
            centroid_x.append(moment_list['m10']/moment_list['m00'])
            centroid_y.append(moment_list['m01']/moment_list['m00'])
            area_list.append(cv2.contourArea(eachContour))

    # Create the voronoi diagram
    centroids = zip(centroid_x,centroid_y)
    vorDiagram = Voronoi(centroids)

    lines = [
        LineString(vorDiagram.vertices[line])
        for line in vorDiagram.ridge_vertices if -1 not in line
    ]
    convex_hull = MultiPoint([Point(i) for i in centroids]).convex_hull.buffer(2)

    result = MultiPolygon(
        [poly.intersection(convex_hull) for poly in shapely.ops.polygonize(lines)])

    plot_grid(111,result)
    return result


########################################################################################################################
# Create the initial conditions
########################################################################################################################

numSites = 50

listPolygons = getPolyFromImage()

# Get numSites points with mean, these will be the diagram sites
#initSites = 5 * np.random.rand(numSites,2) +20


# Create a uniformly distributed grid of points
#initSites =np.random.uniform(0.0, 200.0, size=(35,2))

# Get a regular spaced grid of points
# x = np.linspace(1.0,20.0, num=10)
# y = np.linspace(1.0,20.0,num=10)
# xv, yv = np.meshgrid(x,y)
# initSites= np.concatenate((xv.reshape(100,1),yv.reshape(100,1)), axis=1)
#
# vorDiagram = Voronoi(initSites)
#
# # Initialise figure plot
fig = plt.figure(1, figsize=(5,5), dpi=90)
#
# # This will plot the voronoi diagram
# voronoi_plot_2d(vorDiagram)
# plt.show()
#
# finiteEdges = [
#     shapely.geometry.LineString(vorDiagram.vertices[line])
#     for line in vorDiagram.ridge_vertices
#     if -1 not in line
# ]
#
# # Clean by edge length
# edgeLenList = [line.length for line in finiteEdges]
# print edgeLenList
# edgeLenMean = np.mean(np.array(edgeLenList))
# edgeLenMedian = np.median(np.array(edgeLenList))
# edgesClean = [line for line in finiteEdges if line.length < 2.0*edgeLenMedian]
#
# listPolygons = [poly for poly in shapely.ops.polygonize(edgesClean)]
#
# # clean the polygons, that is, exclude all the polygons that have a greater area than a given deviation
# areaList = [poly.area for poly in listPolygons]
# areaMedian = np.median(np.array(areaList))
# areaMean = np.mean(np.array(areaList))
# polygonClean = [poly for poly in listPolygons if poly.area < 1.2*areaMean]
# areaClean = [poly.area for poly in polygonClean]
#
# # Plot the polygons
ax = fig.add_subplot(111)
for polygon in listPolygons:
        plot_coords(ax, polygon.exterior)
        patch = PolygonPatch(polygon, facecolor=color_isvalid(polygon), edgecolor=color_isvalid(polygon, valid=BLUE), alpha=0.5, zorder=2)
        ax.add_patch(patch)
plt.show()

createInitialFromPoly(listPolygons)
#createInitialFromVoronoi(list)
