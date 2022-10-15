def get_mapping(map_file):
    # Open the file.
    f = open(map_file, "r")

    # Result is a list of dictionaries.
    mapping_list = []

    # Skip the header on the first line.
    header = f.readline()

    for line in f:
        if line == header:
            continue

        else:
            # strip line of whitespace, split elements into list:
            x = line.strip().split()

            # n = number of dictionaries
            n = len(x) - 1

            # assign first item (Ensembl_ID) as value
            value = x[0]

            # create n dictionaries
            while 0 < n:
                if len(mapping_list) < n:
                    d = {}
                    mapping_list.append(d)
                    n -= 1
                else:
                    break

            # fill n dictionaries with n+1th item as key
            for key in list(range(1, len(x))):

                # select correct (n-th) dictionary in mapping_list
                position = key - 1

                if len(mapping_list) < n:
                    mapping_list.append({})
                elif len(mapping_list) == len(x) - 1:
                    # store key and value in n-th dictionary in mapping_list
                    inner_d = {x[key]: value}
                    mapping_list[position].update(inner_d)

    # Remember to close the file after we're done.
    f.close()

    return mapping_list


mapping_list_hsa = get_mapping("mapping/hsa.map")
mapping_list_rno = get_mapping("mapping/rno.map")

def get_go_terms(mapping_list, go_file):
    # Open the file.
    f = open(go_file, "r")

    # This will be the dictionary that this function returns.
    # Entries will have as a key an Ensembl ID and the value will
    # be a set of GO terms.
    go_dict = dict()

    for line in f:

        if line.startswith("!"):
            continue
            # skip comments

        else:
            # strip line of go-file from whitespace and split into list:
            x = line.strip().split()

            # use ID in go file as search word
            query = x[1]

            #find GO term and save as value
            for item in x:
                if item.startswith("GO:"):
                    value = item
                else:
                    continue

            # search for query in every dictionary in mapping list
            for d in mapping_list:
                if query in d.keys():
                    # if found, save Ensembl_ID as key
                    key = mapping_list[mapping_list.index(d)][query]

                    # update dictionary with key:value
                    if key in go_dict:
                        go_dict[key].add(value)
                    else:
                        go_dict[key] = {value}
                else:
                    continue

    # Remember to close the file after we're done.
    f.close()

    return go_dict


go_dict_hsa = get_go_terms(mapping_list_hsa, "GO/hsa.go")
go_dict_rno = get_go_terms(mapping_list_rno, "GO/rno.go")

def compute_score(alignment_file, go_one_dict, go_two_dict):
    # Open the file.
    f = open(alignment_file, "r")

    # Keep track of the number of proteins we can't map to GO terms
    # and the score.
    unmappable_one = 0
    unmappable_two = 0
    score = 0.0

    for line in f:
        # strip line from whitespace and split into list:
        x = line.strip().split()

        # find go-terms for first protein from go_dictionary
        protein1 = x[0]

        if protein1 in go_one_dict:
            protein1_go = go_one_dict[protein1]
        elif protein1 in go_two_dict:
            protein1_go = go_two_dict[protein1]
        else:
            unmappable_one = unmappable_one + 1

        # find go-terms for second protein from go_dictionary
        protein2 = x[1]

        if protein2 in go_one_dict:
            protein2_go = go_one_dict[protein2]
        elif protein2 in go_two_dict:
            protein2_go = go_two_dict[protein2]
        else:
            unmappable_two = unmappable_two + 1

        # compute jaccard score for 2 proteins
        line_score = len(protein1_go.intersection(protein2_go)) / len(protein1_go.union(protein2_go))

        # add line score to global score for file
        score = score + line_score

    # Remember to close the file after we're done.
    f.close()

    # Return the statistics and the score back so the main code
    # can print it out.
    return unmappable_one, unmappable_two, score

compute_score("alignments/rno-hsa.sif", go_dict_hsa, go_dict_rno)