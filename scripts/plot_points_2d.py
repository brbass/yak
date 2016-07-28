import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib import cm as cm

def plot_points_2d(file_path):
    output_xml = et.parse(file_path).getroot()

    dimension = int(output_xml.findtext("./rbf_mesh/dimension"))
    num_points = int(output_xml.findtext("./rbf_mesh/number_of_points"))
    num_materials = int(output_xml.findtext("./rbf_mesh/number_of_materials"))
    points_data = np.fromstring(output_xml.findtext("./rbf_mesh/point_positions"), sep="\t")
    material = np.fromstring(output_xml.findtext("./rbf_mesh/material"), sep="\t")
    points = np.empty((num_points, dimension))
    colors = cm.rainbow(1.0 / num_materials * material)
    
    for i in range(num_points):
        for d in range(dimension):
            points[i, d] = points_data[d + dimension * i]

    plt.figure()
    plt.scatter(points[:,0], points[:,1], c=colors)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis('equal')
    plt.show()
    
if (len(sys.argv) != 2):
    sys.exit()

file_path = sys.argv[1]

plot_points_2d(file_path)
