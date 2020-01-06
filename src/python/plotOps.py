import numpy as np
import sys
from shapely.ops import polygonize
from shapely.geometry import Polygon, MultiPoint, Point
from descartes.patch import PolygonPatch
import math
from os import listdir
import matplotlib.pyplot as plt

########################################################################################################################
# Plotting operations in general
########################################################################################################################

# Replacing the fiugures module
BLUE = '#6699cc'
GRAY = '#999999'
RED = '#ff3333'

GM = (math.sqrt(5)-1.0)/2.0
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

def plot_labels(plot_pos, points):

    # Create the figure
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(plot_pos)

    # Plot point labels
    for i in range(len(points[:,1])):
        ax.text(points[i,0] * (1 + 0.001), points[i,1] * (1 + 0.001) , i+1, fontsize=5)

    return ax

def plot_grid(ax, grid):
    count = 0
    for element in grid:
        count += 1
        polygon = Polygon(element)
        cellCentroid = list(polygon.centroid.coords)
        #print cellCentroid
        plot_coords(ax, polygon.exterior)
        ax.annotate(str(count), xy=cellCentroid[0], fontsize=8)
        col = BLUE

        col = color_isvalid(polygon)

        patch = PolygonPatch(polygon, facecolor=col, edgecolor=color_isvalid(polygon, valid=BLUE), alpha=0.5, zorder=2)
        ax.add_patch(patch)

def flipEdges(polygonCoords):

        x1 = polygonCoords[5][0]
        y1 = polygonCoords[5][1]
        x2 = polygonCoords[6][0]
        y2 = polygonCoords[6][1]

        cx = (x1 + x2) / 2.0
        cy = (y1 + y2) / 2.0

        x1 = ((x1-cx)*math.cos(1.5) + (y1-cy)*math.sin(1.5))+cx
        y1 = (-(x1-cx)*math.sin(1.5) + (y1-cy)*math.cos(1.5))+cy

        x2 = ((x2-cx)*math.cos(1.5) + (y2-cy)*math.sin(1.5))+cx
        y2 = (-(x2-cx)*math.sin(1.5) + (y2-cy)*math.cos(1.5))+cy

        polygonCoords[5] = [x1,y1]
        polygonCoords[6] = [x2,y2]


def getPolygons(filePref):
    with open(filePref, "r") as pointsFile:

        numPoints = int(pointsFile.readline())

        # This is the list of points
        pointsList = np.zeros((numPoints,2))
        for row in range(numPoints):
            pointsList[row,:] = pointsFile.readline().split()

        numCells, numVerticesCell = [int(each) for each in pointsFile.readline().split()]
        polygonList = []
        for row in range(numCells):
            polygonIndex = [int(indPoint)-1 for indPoint in pointsFile.readline().split() if int(indPoint) > 0 ]
            polygonCoords = [list(pointsList[coord,:]) for coord in polygonIndex]
            polygonList.append(polygonCoords)

def getPolygons1(filePref):
    with open(filePref, "r") as pointsFile:

        numPoints = int(pointsFile.readline())

        # This is the list of points
        pointsList = np.zeros((numPoints,2))
        for row in range(numPoints):
            pointsList[row,:] = pointsFile.readline().split()

        numCells = int(pointsFile.readline())
        polygonList = []
        for row in range(numCells):
            polygonIndex = [int(indPoint)-1 for indPoint in pointsFile.readline().split() if int(indPoint) > 0 ]
            polygonCoords = [list(pointsList[coord,:]) for coord in polygonIndex]
            polygonList.append(polygonCoords)

    # Plot the system
    ax = plot_labels(111, pointsList)
    plot_grid(ax, polygonList)
    plt.savefig(filePref+".png", figsize=(1000,1000), format="png")
    plt.cla()

########################################################################################################################
# Interacting with the user
########################################################################################################################
filePref = sys.argv[1]
if sys.argv[2] == "s":
    getPolygons1(filePref)
else:
    fileNames = [eachFile.split(".")[0] for eachFile in listdir(filePref)]
    for eachFile in fileNames:
        getPolygons(filePref+"/"+eachFile)
