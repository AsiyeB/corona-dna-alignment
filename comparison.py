import pickle5 as pickle

global m
m = 4096


def hash(k):
    total_value = 0
    dict_dna = {'A': 0, 'a': 0, 'T': 1, 't': 1, 'C': 2, 'c': 2, 'G': 3, 'g': 3}
    for v in k:
        total_value = 4 * total_value
        if v in dict_dna:
            total_value = total_value + dict_dna[v]

    return total_value


def comparison(cmp_str, ht):
    hash_value = hash(cmp_str)
    index_value = hash_value % m
    llist = None
    try:
        llist = ht[index_value]
    except KeyError:
        return []
    for i in llist:
        if getattr(i, "_Entry__key") == cmp_str and getattr(i, "_Entry__hashcode") == hash_value:
            return getattr(i, "_Entry__indexlist")
    return []


def load_object(filename):
    try:
        with open(filename, "rb") as f:
            return pickle.load(f)
    except Exception as ex:
        print("Error during unpickling object (Possibly unsupported):", ex)


if __name__ == '__main__':
    # f = open("DeltaCoronavirus.fasta", "r")
    f = open("Coronavirus.fasta", "r")
    # [next(f) for i in range(0,19)]
    delta_dna = next(f)
    delta_dna = f.read()
    delta_dna = delta_dna.translate(str.maketrans('', '', ' \n\t\r'))
    obj = load_object("data.pickle")
    ht = getattr(getattr(obj, "_MCLMerFinder__table"), "_HashTable__table")

    print(delta_dna[6751:6770])
    # print()
    # delta_dna2 = delta_dna[22:28]
    # delta_dna3 = delta_dna[28:34]
    # delta_dna4 = delta_dna[34:40]
    # print(delta_dna2)
    #
    #
    print(comparison("GTGTTT", ht))
    # print(delta_dna3)
    # print(comparison(delta_dna3, ht))
    # print(delta_dna4)
    # print(comparison(delta_dna4, ht))
    # print(comparison("TAAACG", ht))
