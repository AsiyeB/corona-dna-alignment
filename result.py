class Result(object):
    # constructor for Entry class
    def __init__(self, score, align1,symbol, align2,identity):
        # initialize __key and __value attributes here
        # __key is the L-mer string
        # __value is the occurence of L-mer string in the DNA strand
        self.a = 0
        self.score = score
        self.align1 = align1
        self.symbol = symbol
        self.align2 = align2
        self.identity = identity
        self.seqA = ""
        self.seqB = ""
