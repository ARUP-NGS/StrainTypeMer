import os
import numpy as np

_base = "/Users/331-SimmonkLPTP/Box Sync/ARUP/ACI_fastqs/"
full = os.path.join(_base, "cookson_full.txt")
shared = os.path.join(_base, "cookson_rapid.txt")


def load_matrix(file_path):
    matrix = []
    parse = False
    labels_x = []
    labels_y = []
    prev_line = ""
    for line in open(file_path):
        if '[SIMILARITY TABLE END]' in line:
            parse = False

        if parse is True:
            if line[0] == ",": # x labels
                labels_x = line[1:].strip().split(",")
            else:  # normal line
                values = line.strip().split(",")
                labels_y.append(values[0])
                matrix.append([float(i) for i in values[1:]])

        if '[SIMILARITY TABLE]' in line and "KMER REFERENCE]" not in prev_line:
            parse = True

        prev_line = line

    if labels_y != labels_x:
        raise NameError('something weird with the labels')

    return matrix, labels_x

def find_indices_with_value(matrix, cutoff_value=98):
    indices = set([])
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] >= 98.0:
                indices.add("{0},{1}".format(i, j))
    return indices


def compare_matrices(matrix_1, matrix_2):
    m1 = np.array(matrix_1, )
    m2 = np.array(matrix_2, )

    print(np.average(m1 - m2))

    #find_indinces_over_98%
    s1 = find_indices_with_value(matrix_1)
    s2 = find_indices_with_value(matrix_2)

    indices = s1.intersection(s2)

    total = 0.0
    for index in indices:
        i, j = [int(v) for v in index.split(",")]
        total += abs(m1[i][j] - m2[i][j])

    print(total/len(indices))






np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
full_matrix, full_labels = load_matrix(full)
shared_matrix, sharded_labels = load_matrix(shared)

if full_labels == sharded_labels:
    compare_matrices(full_matrix, shared_matrix)

