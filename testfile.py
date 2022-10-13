
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


def get_go_terms(mapping_list, go_file):
    # Open the file.
    f = open(go_file, "r")

    # This will be the dictionary that this function returns.
    # Entries will have as a key an Ensembl ID and the value will
    # be a set of GO terms.
    go_dict = dict()

    for line in f:
        # TODO: PUT YOUR CODE HERE



    # Remember to close the file after we're done.
    f.close()

    return go_dict
