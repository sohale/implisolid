# author Solene Chauvier, Marc Fraysse
# The purpose of this file is to compare the MC algorithm and the vertex_resampling part
# in python and in cpp
# we want to be sure that both algorithm give data that are compatible

import numpy as np
# import re

python_data = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo.txt", "r")
cpp_data = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_algo_cpp.txt", "r")
file_test = open("/home/solene/Desktop/mp5-private/solidmodeler/clean_code/data_analysis_summary.txt", "w")

# ************************
# ******** Python ********
# ************************

old_vertex_p = np.zeros((2400, 3))
new_vertex_p = np.zeros((2400, 3))
faces_p = np.zeros((4796, 3))
centroids_p = np.zeros((4796, 3))
norm_old_vertex_p = np.zeros(2400)
norm_new_vertex_p = np.zeros(2400)
norm_centroids_p = np.zeros(4796)
num_neig_vertex_p = np.zeros(2400)
edge_length_o_p = np.zeros(4796)
edge_length_n_p = np.zeros(4796)
# creating the arrays

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


# *************************
# ********  test1  ********
# *************************

norm_old_vertex_p = np.linalg.norm(old_vertex_p, axis=1)
norm_new_vertex_p = np.linalg.norm(new_vertex_p, axis=1)

mean_norm_old_vertex_p = np.mean(norm_old_vertex_p)
mean_norm_new_vertex_p = np.mean(norm_new_vertex_p)

d_m_ir_o_p = abs(mean_norm_old_vertex_p - 1.)
d_m_ir_n_p = abs(mean_norm_new_vertex_p - 1.)

# *************************
# ********  test2  ********
# *************************

for j in range(faces_p.shape[0]):
    num_neig_vertex_p[faces_p[j, 0]] += 1
    num_neig_vertex_p[faces_p[j, 1]] += 1
    num_neig_vertex_p[faces_p[j, 2]] += 1

mean_num_neig_v_p = np.mean(num_neig_vertex_p)

# *************************
# ********  test3  ********
# *************************

norm_centroids_p = np.linalg.norm(centroids_p, axis=1)

mean_norm_centroids_p = np.mean(norm_centroids_p)

d_m_ir_c_p = abs(mean_norm_centroids_p - 1.)


# *************************
# ********  test4  ********
# *************************

for j in range(faces_p.shape[0]):
    edge_length_o_p[j] += np.linalg.norm((old_vertex_p[faces_p[j][0]][0] - old_vertex_p[faces_p[j][1]][0], old_vertex_p[faces_p[j][0]][1] - old_vertex_p[faces_p[j][1]][1], old_vertex_p[faces_p[j][0]][2] - old_vertex_p[faces_p[j][1]][2]))
    edge_length_o_p[j] += np.linalg.norm((old_vertex_p[faces_p[j][0]][0] - old_vertex_p[faces_p[j][2]][0], old_vertex_p[faces_p[j][0]][1] - old_vertex_p[faces_p[j][2]][1], old_vertex_p[faces_p[j][0]][2] - old_vertex_p[faces_p[j][2]][2]))
    edge_length_o_p[j] += np.linalg.norm((old_vertex_p[faces_p[j][2]][0] - old_vertex_p[faces_p[j][1]][0], old_vertex_p[faces_p[j][2]][1] - old_vertex_p[faces_p[j][1]][1], old_vertex_p[faces_p[j][2]][2] - old_vertex_p[faces_p[j][1]][2]))

    edge_length_n_p[j] += np.linalg.norm((new_vertex_p[faces_p[j][0]][0] - new_vertex_p[faces_p[j][1]][0], new_vertex_p[faces_p[j][0]][1] - new_vertex_p[faces_p[j][1]][1], new_vertex_p[faces_p[j][0]][2] - new_vertex_p[faces_p[j][1]][2]))
    edge_length_n_p[j] += np.linalg.norm((new_vertex_p[faces_p[j][0]][0] - new_vertex_p[faces_p[j][2]][0], new_vertex_p[faces_p[j][0]][1] - new_vertex_p[faces_p[j][2]][1], new_vertex_p[faces_p[j][0]][2] - new_vertex_p[faces_p[j][2]][2]))
    edge_length_n_p[j] += np.linalg.norm((new_vertex_p[faces_p[j][2]][0] - new_vertex_p[faces_p[j][1]][0], new_vertex_p[faces_p[j][2]][1] - new_vertex_p[faces_p[j][1]][1], new_vertex_p[faces_p[j][2]][2] - new_vertex_p[faces_p[j][1]][2]))


edge_mean_o_p = np.mean(edge_length_o_p)/3
edge_mean_n_p = np.mean(edge_length_n_p)/3


# *************************
# *********  C++  *********
# *************************

old_vertex_cpp = np.zeros((2406, 3))
new_vertex_cpp = np.zeros((2406, 3))
faces_cpp = np.zeros((4808, 3))
centroids_cpp = np.zeros((4808, 3))
norm_old_vertex_cpp = np.zeros(2406)
norm_new_vertex_cpp = np.zeros(2406)
num_neig_vertex_cpp = np.zeros(2406)
norm_centroids_cpp = np.zeros(4808)
edge_length_o_cpp = np.zeros(4808)
edge_length_n_cpp = np.zeros(4808)
# creating the arrays

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

# *************************
# ********  test1  ********
# *************************

norm_old_vertex_cpp = np.linalg.norm(old_vertex_cpp, axis=1)
norm_new_vertex_cpp = np.linalg.norm(new_vertex_cpp, axis=1)

mean_norm_old_vertex_cpp = np.mean(norm_old_vertex_cpp)
mean_norm_new_vertex_cpp = np.mean(norm_new_vertex_cpp)

d_m_ir_o_cpp = abs(mean_norm_old_vertex_cpp - 0.8)/0.8
d_m_ir_n_cpp = abs(mean_norm_new_vertex_cpp - 0.8)/0.8


# *************************
# ********  test2  ********
# *************************

for j in range(faces_cpp.shape[0]):
    num_neig_vertex_cpp[faces_cpp[j, 0]] += 1
    num_neig_vertex_cpp[faces_cpp[j, 1]] += 1
    num_neig_vertex_cpp[faces_cpp[j, 2]] += 1

mean_num_neig_v_cpp = np.mean(num_neig_vertex_cpp)


# *************************
# ********  test3  ********
# *************************

norm_centroids_cpp = np.linalg.norm(centroids_cpp, axis=1)

mean_norm_centroids_cpp = np.mean(norm_centroids_cpp)

d_m_ir_c_cpp = abs(mean_norm_centroids_cpp - 0.8)/0.8


# *************************
# ********  test4  ********
# *************************

for j in range(faces_cpp.shape[0]):
    edge_length_o_cpp[j] += np.linalg.norm((old_vertex_cpp[faces_cpp[j][0]][0] - old_vertex_cpp[faces_cpp[j][1]][0], old_vertex_cpp[faces_cpp[j][0]][1] - old_vertex_cpp[faces_cpp[j][1]][1], old_vertex_cpp[faces_cpp[j][0]][2] - old_vertex_cpp[faces_cpp[j][1]][2]))
    edge_length_o_cpp[j] += np.linalg.norm((old_vertex_cpp[faces_cpp[j][0]][0] - old_vertex_cpp[faces_cpp[j][2]][0], old_vertex_cpp[faces_cpp[j][0]][1] - old_vertex_cpp[faces_cpp[j][2]][1], old_vertex_cpp[faces_cpp[j][0]][2] - old_vertex_cpp[faces_cpp[j][2]][2]))
    edge_length_o_cpp[j] += np.linalg.norm((old_vertex_cpp[faces_cpp[j][2]][0] - old_vertex_cpp[faces_cpp[j][1]][0], old_vertex_cpp[faces_cpp[j][2]][1] - old_vertex_cpp[faces_cpp[j][1]][1], old_vertex_cpp[faces_cpp[j][2]][2] - old_vertex_cpp[faces_cpp[j][1]][2]))

    edge_length_n_cpp[j] += np.linalg.norm((new_vertex_cpp[faces_cpp[j][0]][0] - new_vertex_cpp[faces_cpp[j][1]][0], new_vertex_cpp[faces_cpp[j][0]][1] - new_vertex_cpp[faces_cpp[j][1]][1], new_vertex_cpp[faces_cpp[j][0]][2] - new_vertex_cpp[faces_cpp[j][1]][2]))
    edge_length_n_cpp[j] += np.linalg.norm((new_vertex_cpp[faces_cpp[j][0]][0] - new_vertex_cpp[faces_cpp[j][2]][0], new_vertex_cpp[faces_cpp[j][0]][1] - new_vertex_cpp[faces_cpp[j][2]][1], new_vertex_cpp[faces_cpp[j][0]][2] - new_vertex_cpp[faces_cpp[j][2]][2]))
    edge_length_n_cpp[j] += np.linalg.norm((new_vertex_cpp[faces_cpp[j][2]][0] - new_vertex_cpp[faces_cpp[j][1]][0], new_vertex_cpp[faces_cpp[j][2]][1] - new_vertex_cpp[faces_cpp[j][1]][1], new_vertex_cpp[faces_cpp[j][2]][2] - new_vertex_cpp[faces_cpp[j][1]][2]))


edge_mean_o_cpp = np.mean(edge_length_o_cpp)/(3*0.8)
edge_mean_n_cpp = np.mean(edge_length_n_cpp)/(3*0.8)


# writing into the file

file_test.write("Mean distance between the computed points and the ideal point of the sphere for the old vertex in python: " + str(d_m_ir_o_p) + "\n")
file_test.write("Mean distance between the computed points and the ideal point of the sphere for the old vertex in cpp:    " + str(d_m_ir_o_cpp) + "\n")
file_test.write("\n")
file_test.write("Mean distance between the computed points and the ideal point of the sphere for the new vertex in python: " + str(d_m_ir_n_p) + "\n")
file_test.write("Mean distance between the computed points and the ideal point of the sphere for the new vertex in cpp:    " + str(d_m_ir_n_cpp) + "\n")
file_test.write("\n")
file_test.write("Mean number of neighbours of vertex in python: " + str(mean_num_neig_v_p) + "\n")
file_test.write("Mean number of neighbours of vertex in cpp:    " + str(mean_num_neig_v_cpp) + "\n")
file_test.write("\n")
file_test.write("Mean distance between the computed points and the ideal point of the sphere for the centroids in python: " + str(d_m_ir_c_p) + "\n")
file_test.write("Mean distance between the computed points and the ideal point of the sphere for the centroids in cpp:    " + str(d_m_ir_c_cpp) + "\n")
file_test.write("\n")
file_test.write("Mean edge length for the old vertex in python: " + str(edge_mean_o_p) + "\n")
file_test.write("Mean edge length for the old vertex in cpp:    " + str(edge_mean_o_cpp) + "\n")
file_test.write("\n")
file_test.write("Mean edge length for the new vertex in python: " + str(edge_mean_n_p) + "\n")
file_test.write("Mean edge length for the new vertex in cpp:    " + str(edge_mean_n_cpp) + "\n")

python_data.close()
cpp_data.close()
file_test.close()
