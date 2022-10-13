import sys

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


def compute_score(alignment_file, go_one_dict, go_two_dict):
    # Open the file.
    f = open(alignment_file, "r")

    # Keep track of the number of proteins we can't map to GO terms
    # and the score.
    unmappable_one = 0
    unmappable_two = 0
    score = 0.0

    for line in f:
        # TODO: PUT YOUR CODE HERE


    # Remember to close the file after we're done.
    f.close()

    # Return the statistics and the score back so the main code
    # can print it out.
    return unmappable_one, unmappable_two, score


def main():
    # TODO: PUT YOUR CODE HERE


if __name__ == '__main__':
    main()