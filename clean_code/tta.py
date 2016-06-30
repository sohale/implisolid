# author Solene Chauvier, Marc Fraysse
# The purpose of this file is to compare the MC algorithm and the vertex_resampling part
# in python and in cpp
# we want to be sure that both algorithm give data that are compatible

import numpy as np
# import re

python_data = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo.txt", "r")
cpp_data = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo_cpp.txt", "r")
file_test = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_analysis_summary.txt", "w")


old_vertex_p = np.zeros((2400, 3))
new_vertex_p = np.zeros((2400, 3))
faces_p = np.zeros((4796, 3))
centroids_p = np.zeros((4796, 3))

i = 0
for line_p in python_data:
    line_p = line_p.strip()
    line_p = line_p.replace("[", "")
    line_p = line_p.replace("]", "")
    i += 1
    if (i >= 2 and i <= 2401):
        line_p = line_p.replace("      ", " ")
        vertex = line_p.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                old_vertex_p[i-2, k] = float(vertex[j])
                k += 1
    elif(i >= 2404 and i <= 4803):
        line_p = line_p.replace("      ", " ")
        vertex = line_p.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                new_vertex_p[i-2404, k] = float(vertex[j])
                k += 1
    elif(i >= 4805 and i <= 9600):
        line_p = line_p.replace("      ", " ")
        vertex = line_p.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                faces_p[i-4805, k] = int(vertex[j])
                k += 1
    elif(i >= 9602 and i <= 14397):
        line_p = line_p.replace("      ", " ")
        vertex = line_p.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                centroids_p[i-9602, k] = float(vertex[j])
                k += 1

old_vertex_cpp = np.zeros((2406, 3))
new_vertex_cpp = np.zeros((2406, 3))
faces_cpp = np.zeros((4808, 3))
centroids_cpp = np.zeros((4808, 3))

i = 0
for line_cpp in cpp_data:
    line_cpp = line_cpp.strip()
    i += 1
    if (i >= 2 and i <= 2407):
        vertex = line_cpp.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                old_vertex_cpp[i-2, k] = float(vertex[j])
                k += 1
    elif(i >= 2410 and i <= 4815):
        vertex = line_cpp.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                new_vertex_cpp[i-2410, k] = float(vertex[j])
                k += 1
    elif(i >= 4817 and i <= 9624):
        vertex = line_cpp.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                faces_cpp[i-4817, k] = int(vertex[j])
                k += 1
    elif(i >= 9626 and i <= 14433):
        vertex = line_cpp.split(" ")
        k = 0
        for j in range(len(vertex)):
            if vertex[j] != " " and vertex[j] != "":
                centroids_cpp[i-9626, k] = float(vertex[j])
                k += 1


python_data.close()
cpp_data.close()
file_test.close()
