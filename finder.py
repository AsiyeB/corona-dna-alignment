"""
Author:         Elaine Chan Yun Ru
Assignment:     PA2_HashTables
Date:           2/12/2021
Description:    contains the class MCLMerFinder - to use other classes to search for the most common L-mer string that is in a dna string
"""
# import statements here
from hashtable import HashTable


class MCLMerFinder(object):
    # constructor for MCLMerFinder class
    def __init__(self):
        # initialize attribute __table
        self.__table = HashTable()

    def findMCLMer(self, dna_string, I):
        # split up the dna_string into L-mer strings according to the provided number of character per substring
        dna_len = len(dna_string)
        split_strings = []  # to store each substring of the dna_strand which is split up accordingly
        split_index = []
        index_value_list = []  # to store the index value of each substring (index value of substring is the remainder of the total decimal ASCII value of the substring by the amount of substrings)
        res_index_value_list = []  # to store the index value of L-mer string without duplicates
        value_list = dict()
        m = 4096  # represents the amount of substrings per dna strand
        for index in range(0, dna_len, I):
            split_strings.append(dna_string[index:index + I])
            split_index.append(index)
            # m+=1

        # calculate the index value per substring and store into index_value_list
        size = len(split_strings)
        for i in range(size):
            hash_value = self.__table.hash(split_strings[i])
            index_value = hash_value % m
            index_value_list.append(index_value)
            self.__table.insert(index_value, split_strings[i], hash_value, split_index[i])
            # store the index value of the substrings without duplicates into res_index_value_list
            # [res_index_value_list.append(x) for x in index_value_list if x not in res_index_value_list]
            # f = open("result.txt", "w")
            # for i in range(len(res_index_value_list)):
            #     llist = self.__table.lookup(res_index_value_list[i])
            #     for i in llist:
            #         keyhash = str(getattr(i,"_Entry__key")) + " => " + str(getattr(i,"_Entry__hashcode")) + str(getattr(i,"_Entry__indexlist")) + "\n"
            #         f.write(keyhash)
            # f.close()
