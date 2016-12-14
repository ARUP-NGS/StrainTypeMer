from straintypemer.sub_commands.plots import generage_matrix


def plot_output(input_file, output_prefix):
    results = parse_results(input_file)
    for key, result in results.items():
        # if "full" in key:
            # print(result[0])
            # print(result[1])
            # print(len(result[2]), len(result[2][0]))
            #print(result[2][0])
            #generage_matrix(result[0], result[1], result[2], output_prefix + "_" + key, result[3],
                            #identical=98.5, possibly_related=95.0, different=0) # (x_labels, y_labels, matrix, kmer_counts)

            if "inverse" in key:
                pass
                # generage_matrix(result[0], result[1], result[2], output_prefix + "_" + key, result[3],
                #                 identical=98.0, possibly_related=90.0, different=0) # (x_labels, y_labels, matrix, kmer_counts)
            elif "kmer_reference" == key:
                pass
                # generage_matrix(result[0], result[1], result[2], output_prefix + "_" + key, result[3],
                #                 identical=99.9, possibly_related=99.0, different=0)
            else:
                # generage_matrix(result[0], result[1], result[2], output_prefix + "_" + key, result[3],
                #                 identical=99.5, possibly_related=95.0, different=0)
                generage_matrix(result[0], result[1], result[2], "test", result[3],
                                identical=None, possibly_related=None, different=None)



def parse_results(input_file):
    input_text = iter(input_file.readlines())
    kmer_reference_matrix_data = []
    inverse_reference_matrix_data = []
    full_matrix_data = []
    try:
        while input_text:
            line = input_text.__next__()
            if "[KMER REFERENCE]" in line:
                kmer_reference = True
                while kmer_reference:
                    if "[DENOMINATOR TABLE END]" in line:
                        kmer_reference = False
                    else:
                        line = input_text.__next__()
                        kmer_reference_matrix_data.append(line)


            if "[INVERSE KMER REFERENCE]" in line:
                inverse_reference = True
                while inverse_reference:
                    if "[DENOMINATOR TABLE END]" in line:
                        inverse_reference = False
                    else:
                        line = input_text.__next__()
                        inverse_reference_matrix_data.append(line)

            if "[FULL]" in line:
                full_matrix = True
                while full_matrix:
                    if "[DENOMINATOR TABLE END]" in line:
                        full_matrix = False
                    else:
                        line = input_text.__next__()
                        full_matrix_data.append(line)
    except StopIteration:
        pass

    results = {}
    if kmer_reference_matrix_data:
        results.update({"kmer_reference" : parse_tables(kmer_reference_matrix_data)})

    if inverse_reference_matrix_data:
        results.update({"inverse_kmer_reference": parse_tables(inverse_reference_matrix_data)})

    if full_matrix_data:
        results.update({"full": parse_tables(full_matrix_data)})


    return results



def parse_tables(raw_matrix):
    in_identity = True

    x_labels = None
    y_labels = []
    matrix = []
    denominator_matrix = []

    for line in raw_matrix:
        if in_identity:
            if "[SIMILARITY TABLE]" not in line:
                if "[SIMILARITY TABLE END]" in line:
                    in_identity = False
                elif "," == line[0]:  # xlabels
                    x_labels = line.strip().split(",")[1:]

                else:
                    _values = line.strip().split(",")
                    if _values != [""]:
                        y_labels.append(_values[0])
                        matrix.append([float(i) for i in _values[1:]])

        else: # in the denominator table
                if "[" not in line and "]" not in line:
                    if "," == line[0]:  # xlabels
                        pass
                    else:
                        _values = line.strip().split(",")
                        if _values != [""]:
                            denominator_matrix.append(_values[1:])


    #print(len(x_labels), x_labels[-1])
    #print(len(y_labels), y_labels[-1])
    #print(denominator_matrix)
    kmer_counts = {}
    for i in range(len(denominator_matrix)):
        for j in range(len(denominator_matrix)):
            if i == j:
                kmer_counts.update({x_labels[i] : int(denominator_matrix[i][j])})

    #cleanup matrix to remove rescue values
    clean_matrix = matrix.copy()
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i < j:
                clean_matrix[i][j] = matrix[i][j]
            elif i == j:
                clean_matrix[i][j] = matrix[i][j]
            elif i > j:
                clean_matrix[j][i] = matrix[i][j]


    return x_labels, y_labels, clean_matrix, kmer_counts


def main():
    input_file = "/Users/331-SimmonkLPTP/Box Sync/ARUP/info/hauser_straintypemer_remove_poor_kmer_refernce.txt"
    plot_output(open(input_file, "r"), "/Users/331-SimmonkLPTP/Box Sync/ARUP/info/hauser_straintypemer_remove_poor_kmer_reference_gradiant")

if __name__ == "__main__":
    main()