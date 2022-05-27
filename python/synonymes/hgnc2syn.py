import shlex
from collections import Counter
import os, sys

sys.path.append(os.path.dirname(__file__) + "/../")


from utils.DataFrame import DataFrame

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import printToFile, dataDir, loadExludeWords


if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-t', '--tsv', type=argparse.FileType("r"), required=True, help="input genes file")
    parser.add_argument('-s', '--syn', type=argparse.FileType("w"), required=True, help="output synonym file")
    args = parser.parse_args()


    hgncData = DataFrame.parseFromFile(args.tsv.name)
    hgncIDIdx = hgncData.getColumnIndex("HGNC ID")
    hgncSymIdx= hgncData.getColumnIndex("Approved symbol")
    hgncNameIdx= hgncData.getColumnIndex("Approved name")
    hgncPrevSymIdx= hgncData.getColumnIndex("Previous symbols")
    hgncPrevNameIdx= hgncData.getColumnIndex("Previous name")
    hgncSynsIdx= hgncData.getColumnIndex("Alias symbols")
    hgncNameSynIdx= hgncData.getColumnIndex("Alias names")

    vAllSyns = []
    dAllSyns = {}

    synIDCounter = Counter()

    for result in hgncData:

        hgncID = result[hgncIDIdx]
        mirexploreGeneID = result[hgncSymIdx]

        if mirexploreGeneID == None:
            print("No symbol: " + str(hgncID))
            continue


        if "~WITHDRAWN" in mirexploreGeneID.upper():
            pass

        if hgncID == "HGNC:2434":
            print(hgncID)

        mirexploreGeneID = mirexploreGeneID.upper().replace('~WITHDRAWN', '').replace(':', '_')

        synIDCounter[mirexploreGeneID] += 1


        if len(mirexploreGeneID) == 1:
            mirexploreGeneID = hgncID

        result = [x if not x is None else "" for x in result]
        #syn = Synonyme(hgncID)
        syn = Synonym( mirexploreGeneID )
        if not "WITHDRAWN" in result[hgncIDIdx].upper():
            syn.addTextSyns(result[hgncIDIdx], addBrackets=False)

        if not "WITHDRAWN" in result[hgncSymIdx].upper():
            syn.addTextSyns('"'+result[hgncSymIdx]+'"', addBrackets=False)
        
        if not "WITHDRAWN" in result[hgncNameIdx].upper():
            syn.addTextSyns('"'+result[hgncNameIdx]+'"', addBrackets=False)

        if not "WITHDRAWN" in result[hgncPrevSymIdx].upper():
            syn.addTextSyns(result[hgncPrevSymIdx], addBrackets=False)

        if not "WITHDRAWN" in result[hgncPrevNameIdx].upper():
            syn.addTextSyns(result[hgncPrevNameIdx], captureBrackets=False)
        
        if not "WITHDRAWN" in result[hgncSynsIdx].upper():
            syn.addTextSyns(result[hgncSynsIdx], addBrackets=False)
        
        if not "WITHDRAWN" in result[hgncNameSynIdx].upper():
            syn.addTextSyns(result[hgncNameSynIdx], addBrackets=False)

        vAllSyns.append(syn)

        if not mirexploreGeneID in dAllSyns:
            dAllSyns[mirexploreGeneID] = syn
        else:
            existingSyn = dAllSyns[mirexploreGeneID]

            for x in existingSyn.syns:
                syn.addTextSyns(x, addBrackets=False)

            dAllSyns[mirexploreGeneID] = syn

    vAllSyns = [dAllSyns[x] for x in dAllSyns]

    for x in synIDCounter:
        if synIDCounter[x] > 1:
            print("Attention: double syn id! " + str(x) + " -> " + str(synIDCounter[x]))

    for syn in vAllSyns:

        removeSyns = []
        for synword in syn.syns:

            if len(synword) == 1:
                removeSyns.append(synword)

        syn.removeSyn(removeSyns)


    globalKeywordExcludes = loadExludeWords(cell_co=False, common=False, generic=True, syngrep=False)
    vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=500, maxCommonCount=10, addAlphaBeta=True, addHyphenGene=True, removeSyn=lambda synonym: synonym.id.startswith('MIR') and not synonym.id.endswith('HG'))
    
    printToFile(vPrintSyns, args.syn.name, codec="utf8")

