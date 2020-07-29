import shlex
from collections import Counter

from utils.DataFrame import DataFrame

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import printToFile, dataDir, loadExludeWords

hgncData = DataFrame.parseFromFile(dataDir + "/miRExplore/hgnc.tsv")
hgncIDIdx = hgncData.getColumnIndex("HGNC ID")
hgncSymIdx= hgncData.getColumnIndex("Approved Symbol")
hgncNameIdx= hgncData.getColumnIndex("Approved Name")
hgncPrevSymIdx= hgncData.getColumnIndex("Previous Symbols")
hgncPrevNameIdx= hgncData.getColumnIndex("Previous Name")
hgncSynsIdx= hgncData.getColumnIndex("Synonyms")
hgncNameSynIdx= hgncData.getColumnIndex("Name Synonyms")

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

    #if mirexploreGeneID == 'RN7SL6P':
    #    print(result)


    if len(mirexploreGeneID) == 1:
        mirexploreGeneID = hgncID

    #syn = Synonyme(hgncID)
    syn = Synonym( mirexploreGeneID )
    syn.addTextSyns(result[hgncIDIdx], addBrackets=False)
    syn.addTextSyns('"'+result[hgncSymIdx]+'"', addBrackets=False)
    syn.addTextSyns('"'+result[hgncNameIdx]+'"', addBrackets=False)
    syn.addTextSyns(result[hgncPrevSymIdx], addBrackets=False)
    syn.addTextSyns(result[hgncPrevNameIdx], captureBrackets=False)
    syn.addTextSyns(result[hgncSynsIdx], addBrackets=False)
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
            print(syn.id, syn.syns)

            removeSyns.append(synword)


globalKeywordExcludes = loadExludeWords(cell_co=False, common=False, generic=True, syngrep=False)
vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=500, maxCommonCount=10, addAlphaBeta=True, addHyphenGene=True, removeSyn=lambda synonym: synonym.id.startswith('MIR') and not synonym.id.endswith('HG'))
#printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/hgnc.syn", codec='utf8')
printToFile(vPrintSyns, "/mnt/d/dev/data/pmid_jun2020/synonyms/hgnc.syn", codec="utf8")

