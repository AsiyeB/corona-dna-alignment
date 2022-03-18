"""
Author:         Elaine Chan Yun Ru
Assignment:     PA2_HashTables
Date:           2/12/2021
Description:    contains the class Entry - to set key-value pairs of the occurence of L-mer
"""

class Entry(object):
    #constructor for Entry class
    def __init__(self):
        #initialize __key and __value attributes here
        #__key is the L-mer string
        #__value is the occurence of L-mer string in the DNA strand
        self.__key = ""
        self.__value = 0
        self.__hashcode = 0
        self.__indexlist = [-1]

