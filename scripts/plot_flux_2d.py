import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib import cm as cm

def plot_flux_2d(file_path):
    output_xml = et.parse(file_path).getroot()
    
    dimension = int(output_xml.findtext("./rbf_mesh/dimension"))
    num_points = int(output_xml.findtext("./rbf_mesh/number_of_points"))
    num_groups = int(output_xml.findtext("./energy_discretization/number_of_groups"))
    points_data = np.fromstring(output_xml.findtext("./rbf_mesh/point_positions"), sep="\t")
    phi_data = np.fromstring(output_xml.findtext("./solution/phi"), sep="\t")
    
    points = np.empty((num_points, dimension))
    phi = np.empty((num_points, num_groups))
    
    for i in range(num_points):
        for d in range(dimension):
            points[i, d] = points_data[d + dimension * i]

    for i in range(num_points):
        for g in range(num_groups):
            phi[i, g] = phi_data[g + num_groups * i]

    for g in range(num_groups):
        plt.figure(g)
        plt.tricontourf(points[:, 0], points[:, 1], phi[:, g])
        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("group " + str(g))
        plt.colorbar()
        plt.axis('equal')
        plt.tight_layout()
        
    plt.show()

if (len(sys.argv) != 2):
    sys.exit()

file_path = sys.argv[1]

plot_flux_2d(file_path)
