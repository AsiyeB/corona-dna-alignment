"""
Author:         Elaine Chan Yun Ru
Assignment:     PA2_HashTables
Date:           2/12/2021
Description:    Main function of program is here - to run the program
"""
# import statements here
from finder import MCLMerFinder
from hashtable import HashTable
from comparison import load_object
from comparison import comparison
from nw import needle
import pickle5 as pickle
from result import Result
import asyncio

splitStep = 6
minSimilarityCount = 24
loop = asyncio.get_event_loop()
needleManLength = 102

f = open("Coronavirus.fasta", "r")
next(f)
output = f.read()
coronaVirus = output.translate(str.maketrans('', '', ' \n\t\r'))

# todo first time
obj = load_object("data.pickle")
ht = getattr(getattr(obj, "_MCLMerFinder__table"), "_HashTable__table")


async def main():
    coronaVirusDna = readFile("Coronavirus.fasta")
    hashCoronaVirus(coronaVirusDna)
    delta_dna = readFile("DeltaCoronavirus.fasta")
    coronaVirusDnaLength = (len(coronaVirusDna))
    print("computing....")
    result2 = divideAndConquer(delta_dna, range(0, coronaVirusDnaLength))
    result = await result2
    print("compute done")
    print("saving in file...")
    saveResult("coronaVirusAlignment.txt", result.align1)
    saveResult("deltaCoronaVirusAlignment.txt", result.align2)
    saveResult("symbol.txt", result.symbol)
    print("All done")
    print("Score", result.score)
    identity = float(result.identity) / len(result.align1) * 100
    print(identity - len(result.align1))
    print('Identity =', "%3.3f" % identity, 'percent')


def saveResult(path, str):
    f = open(path, "w")
    step = 70
    for i in range(0, len(str), step):
        f.write(str[i:i + step] + "\n")
    f.close()


async def divideAndConquer(deltaStr, originalCoronaRange):
    originalCoronaLen = originalCoronaRange.stop - originalCoronaRange.start + 1
    if len(deltaStr) <= needleManLength or originalCoronaLen <= needleManLength:
        return needle(deltaStr, coronaVirus[originalCoronaRange.start:originalCoronaRange.stop])

    baseSelectedDeltaIndex = int(len(deltaStr) / 2)
    i = 0
    sign = 1
    while i < len(deltaStr) / 2:
        if sign == -1:
            i += 1
            sign = 1
        else:
            sign = -1
        selectedDeltaIndex = baseSelectedDeltaIndex + sign * i
        if selectedDeltaIndex < 0 or selectedDeltaIndex + splitStep > len(deltaStr):
            continue
        item = deltaStr[selectedDeltaIndex:selectedDeltaIndex + splitStep]
        coronaIndexList = comparison(item, ht)
        coronaNearestIndex = getNearestIndex(originalCoronaRange.start + selectedDeltaIndex, coronaIndexList)
        if coronaNearestIndex == -1:
            continue

        sameRangeDelta = extendSubStr(deltaStr, selectedDeltaIndex, originalCoronaRange, coronaNearestIndex)

        if sameRangeDelta is None:
            continue
        elif sameRangeDelta.stop - sameRangeDelta.start < needleManLength:
            selectBestHundredRange(sameRangeDelta, 0, 0)

        startCorona = coronaNearestIndex - (selectedDeltaIndex - sameRangeDelta.start)
        endCorona = coronaNearestIndex + (sameRangeDelta.stop - selectedDeltaIndex)

        #
        # before = divideAndConquer(deltaStr[0:sameRangeDelta.start], range(originalCoronaRange.start, startCorona))
        # result = needleMan(deltaStr[sameRangeDelta.start:sameRangeDelta.stop + 1], range(startCorona, endCorona))
        # afterr = divideAndConquer(deltaStr[sameRangeDelta.stop + 1:len(deltaStr)],
        #                           range(endCorona, originalCoronaRange.stop))
        asyncOut = await asyncio.gather(
            divideAndConquer(deltaStr[0:sameRangeDelta.start], range(originalCoronaRange.start, startCorona)),
            needleMan(deltaStr[sameRangeDelta.start:sameRangeDelta.stop + 1], range(startCorona, endCorona)),
            divideAndConquer(deltaStr[sameRangeDelta.stop + 1:len(deltaStr)],
                             range(endCorona, originalCoronaRange.stop))
        )
        before = asyncOut[0]
        result = asyncOut[1]
        afterr = asyncOut[2]
        return Result(
            before.score + result.score + afterr.score,
            before.align1 + result.align1 + afterr.align1,
            before.symbol + result.symbol + afterr.symbol,
            before.align2 + result.align2 + afterr.align2,
            before.identity + result.identity + afterr.identity,
        )

    print("not Found!!")
    return Result
    # getattr(getattr(obj, "_MCLMerFinder__table"), "_HashTable__table")


async def needleMan(deltaList, originalCoronaRange):
    return needle(deltaList, coronaVirus[originalCoronaRange.start:originalCoronaRange.stop])


def selectBestHundredRange(oldRange, deltaStr, coronRange):
    # # //TODo Maybe throw exeption
    padding = oldRange.stop - oldRange.start - needleManLength
    # score = 0
    # outputRange = oldRange
    # for i in range(0, padding):
    #     start = oldRange.start - i
    #     end = oldRange.stop + padding - i
    #     r = needleMan(deltaStr[start,end],coronRange)
    #     if r.score > score or score==0:
    #         score = r.score
    #         outputRange = range(start,end)
    #
    # return outputRange
    if oldRange.start - padding / 2 > 0:
        return range(int(oldRange.start - padding / 2),
                     int(oldRange.start - padding / 2 + int(needleManLength)))
    else:
        return range(0, int(needleManLength))


def extendSubStr(deltaList, deltaIndex, originalCoronaRange, coronaIndex):
    upDeltaIndex = deltaIndex
    downDeltaIndex = deltaIndex
    upCoronaIndex = coronaIndex
    downCoronaIndex = coronaIndex
    while True:  # for set Up
        if (upDeltaIndex + splitStep >= len(deltaList)) or (upCoronaIndex + splitStep > originalCoronaRange.stop):
            break

        getUpCoronaList = comparison(deltaList[upDeltaIndex + splitStep: upDeltaIndex + 2 * splitStep], ht)

        if isExistIndex(upCoronaIndex + splitStep, getUpCoronaList):
            upDeltaIndex += splitStep
            upCoronaIndex += splitStep
        else:
            break

    while True:  # for set Down
        if downDeltaIndex - splitStep < 0 or downCoronaIndex - splitStep < originalCoronaRange.start:
            break

        getDownCoronaList = comparison(deltaList[downDeltaIndex - splitStep: downDeltaIndex], ht)
        if isExistIndex(downCoronaIndex - splitStep, getDownCoronaList):
            downDeltaIndex -= splitStep
            downCoronaIndex -= splitStep
        else:
            break
    if upDeltaIndex - downDeltaIndex >= minSimilarityCount:

        return range(downDeltaIndex, upDeltaIndex)
    else:
        return None


def isExistIndex(goalIndex, indexList):
    for i in indexList:
        if i == goalIndex:
            return True
    return False


def getNearestIndex(goalIndex, indexList):
    minDiff = 6000
    selectedIndex = -1
    for i in indexList:
        diff = abs(goalIndex - i)
        if minDiff > diff:
            minDiff = diff
            selectedIndex = i

    return selectedIndex


def readFile(filePath):
    f = open(filePath, "r")
    next(f)
    output = f.read()
    output = output.translate(str.maketrans('', '', ' \n\t\r'))
    # delta_dna_array = spliter(delta_dna, splitStep)
    return output


def splitter(sourceStr, step):
    return [(sourceStr[i:i + step]) for i in range(0, len(sourceStr), step)]


def hashCoronaVirus(sourceDnaStr):
    # hardcode a dna string and the L value
    # create instances of MCLMerFinder and HashTable
    lmer_finder = MCLMerFinder()
    hash_table = HashTable()
    table = dict()

    # make sure each element has an empty list
    setattr(hash_table, '_HashTable__table', table)
    setattr(lmer_finder, '_MCLMerFinder__table', hash_table)

    # find and print most common L-mer in dna
    lmer_finder.findMCLMer(sourceDnaStr, splitStep)

    try:
        with open("data.pickle", "wb") as f:
            pickle.dump(lmer_finder, f, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as ex:
        print("Error during pickling object (Possibly unsupported):", ex)

    # clear all values when done
    getattr(lmer_finder, '_MCLMerFinder__table').clear()


if __name__ == '__main__':
    # s = "TTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCAC "
    # s2 = splitter(s,splitStep)
    #
    # b = comparison(s2[82],ht)
    # print(b)
    # a = getNearestIndex(82,b)
    # print(a)
    # c = extendSubStr(s2,82,range(0,250),a)
    # print(c)
    # print(len(s2))

    loop.run_until_complete(main())
    loop.close()
    # main()
