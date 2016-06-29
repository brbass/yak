import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib import cm as cm
from mpl_toolkits.mplot3d import Axes3D

def plot_points_3d(file_path):
    output_xml = et.parse(file_path).getroot()

    dimension = int(output_xml.findtext("./rbf_mesh/dimension"))
    num_points = int(output_xml.findtext("./rbf_mesh/number_of_points"))
    points_data = np.fromstring(output_xml.findtext("./rbf_mesh/point_positions"), sep="\t")
    material = np.fromstring(output_xml.findtext("./rbf_mesh/material"), sep="\t")
    points = np.empty((num_points, dimension))
    colors = cm.rainbow(material)

    if (len(points_data) != dimension * num_points):
        print("error: points data has wrong size")
    
    for i in range(num_points):
        for d in range(dimension):
            points[i, d] = points_data[d + dimension * i]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:,0], points[:,1], points[:,2], c=colors)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.axis('equal')
    
    plt.show()
    
if (len(sys.argv) != 2):
    sys.exit()

file_path = sys.argv[1]

plot_points_3d(file_path)
