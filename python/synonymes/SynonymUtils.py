from collections import Counter


def handleCommonExcludeWords(synList, excludeWords, mostCommonCount=66, maxCommonCount=10, addAlphaBeta=False, addHyphenGene=False, removeSyn=None, minSynCount=0):

    synCounter = Counter()
    for synonym in synList:

        for syn in synonym:
            synCounter[syn] += 1

    setCommonWords = set()

    for synWordCount in synCounter.most_common(mostCommonCount):

        if synCounter[synWordCount[0]] <= maxCommonCount:
            continue

        setCommonWords.add(synWordCount[0])
        print(synWordCount[0] + " " + str(synWordCount[1]))

    vPrintSyns = []


    for synonym in synList:

        synonym.removeCommonSynonymes(setCommonWords)
        synonym.removeNumbers()
        synonym.removeSynUpper(excludeWords)


        if addHyphenGene:
            synonym.addHyphenVariants()

        if addAlphaBeta:
            synonym.addAlphaBetaVariants()


        if removeSyn:
            if removeSyn(synonym):
                continue

        if len(synonym) > minSynCount:
            vPrintSyns.append(synonym)

    return vPrintSyns