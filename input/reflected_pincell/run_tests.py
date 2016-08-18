import os, sys, subprocess, itertools, multiprocessing, functools
import numpy as np
import xml.etree.ElementTree as et
from matplotlib import pyplot as plt

########################################################
# Run owl in parallel for many values of shape parameter
# Output selected results to file
########################################################

# point to files and commands

command = "~/Documents/research/code/yak/bin/yak"
# command = "echo"

io_folder = "./io/"

template_filename = "template.xml"
template_file = open(template_filename, "r")
template_string = template_file.read()
template_file.close()

results_filename = "output.txt"

# set problem values

shape_values = np.concatenate((np.around(np.arange(0.10, 0.3, 0.01), decimals=2), np.around(np.arange(0.3, 3.1, 0.1), decimals=1)))
distance_values = np.around(np.arange(0.05, 0.15, 0.05), decimals=2)
trans_values = [0, 1]
input_filenames = [io_folder + "pin_{1}_{0}_{2}.xml".format(i, j, k) for i in shape_values for j in distance_values for k in trans_values]
output_filenames = [io_folder + "pin_{1}_{0}_{2}.xml.out".format(i, j, k) for i in shape_values for j in distance_values for k in trans_values]
output_values = [(i, j, k) for i in shape_values for j in distance_values for k in trans_values]
number_of_cases = len(input_filenames)
input_indices = range(number_of_cases)

for i, shape in enumerate(shape_values):
    for j, distance in enumerate(distance_values):
        for k, trans in enumerate(trans_values):
            input_string = template_string.replace("(SHAPE)", str(shape)).replace("(DISTANCE)", str(distance)).replace("(TRANS)", str(trans))
            input_file = open(input_filenames[k + len(trans_values) * (j + len(distance_values) * i)], "w")
            input_file.write(input_string)
            input_file.close()

# storage arrays

k_values = np.empty(number_of_cases)

# functions

def run_problem(index):
    test_description = "\tInd: {2} / {3}\tDist: {0}\tShape: {1}\tTrans: {4}".format(output_values[index][1], output_values[index][0], index + 1, number_of_cases, output_values[index][2])
    print("Begin test" + test_description)
    subprocess.call(["{0} {1}".format(command, input_filenames[index])], shell=True)
    try:
        output_xml = et.parse(output_filenames[index]).getroot()
        total_iterations = int(output_xml.findtext("krylov_iteration/total_iterations"))
        time = float(output_xml.findtext("solution/time"))
        k_eigenvalue = float(output_xml.findtext("./solution/k_eigenvalue"))
        k_values[index] = k_eigenvalue
        
        results_string = str(output_values[index][0]) + "\t" + str(output_values[index][1]) + "\t" + str(output_values[index][2]) + "\t" + str(k_eigenvalue) + "\t" + str(total_iterations) + "\t" + str(time) + "\n"
        
        lock.acquire()
        
        f = open(results_filename, 'a')
        
        f.write(results_string)
        f.close()
        
        lock.release()
        
        print("End test" + test_description + "\tEigenvalue: " + str(k_eigenvalue))
    except:
        results_string = str(output_values[index][0]) + "\t" + str(output_values[index][1]) + "\t" + str(output_values[index][2]) + "\tfailed to execute" + "\n"
        
        lock.acquire()
        
        f = open(results_filename, 'a')
        
        f.write(results_string)
        f.close()
        
        lock.release()
        
        print("Fail test" + test_description)
        
def init_pool(l):
    global lock
    lock = l

def init_results_file():
    result_string = "shape\tdist\ttrans\tk_eff\tkryl_iter\ttime\n"
    
    f = open(results_filename, 'a')
    f.write(result_string)
    f.close()    

# perform solution and postprocessing
    
if __name__ == '__main__':
    init_results_file()
    
    count = 3
    lock_tmp = multiprocessing.Lock()
    pool = multiprocessing.Pool(processes=count, initializer=init_pool, initargs=(lock_tmp,))
    pool.map(run_problem, input_indices)
    pool.close()
    pool.join()
