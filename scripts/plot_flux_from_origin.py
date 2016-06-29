import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt
from matplotlib import cm as cm

def plot_flux_from_origin(file_path):
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

    distance_from_origin = np.empty(num_points)

    for i in range(num_points):
        distance_from_origin[i] = np.linalg.norm(phi[i, :])

    indices = np.argsort(distance_from_origin)

    distance_from_origin = distance_from_origin[indices]
    for g in range(num_groups):
        phi[:, g] = phi[indices, g]
    
    colors = cm.rainbow(range(num_groups))

    plt.figure()
    
    for g in range(num_groups):
        plt.plot(distance_from_origin, phi[:, g], c=colors[g], label=g)
    plt.xlabel("distance from origin")
    plt.ylabel("scalar flux")
    plt.legend()
    
    plt.show()
        
if (len(sys.argv) != 2):
    sys.exit()

file_path = sys.argv[1]

plot_flux_from_origin(file_path)
