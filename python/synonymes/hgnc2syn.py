import shlex
from collections import Counter

from porestat.utils.DataFrame import DataFrame

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

    mirexploreGeneID = mirexploreGeneID.upper().replace('~WITHDRAWN', '').replace(':', '_')

    synIDCounter[mirexploreGeneID] += 1

    #if mirexploreGeneID == 'RN7SL6P':
    #    print(result)

    #syn = Synonyme(hgncID)
    syn = Synonym( mirexploreGeneID )
    syn.addTextSyns(result[hgncIDIdx])
    syn.addTextSyns('"'+result[hgncSymIdx]+'"')
    syn.addTextSyns('"'+result[hgncNameIdx]+'"')
    syn.addTextSyns(result[hgncPrevSymIdx])
    syn.addTextSyns(result[hgncPrevNameIdx])
    syn.addTextSyns(result[hgncSynsIdx])
    syn.addTextSyns(result[hgncNameSynIdx])

    vAllSyns.append(syn)

    if not mirexploreGeneID in dAllSyns:
        dAllSyns[mirexploreGeneID] = syn
    else:
        existingSyn = dAllSyns[mirexploreGeneID]

        for x in existingSyn.syns:
            syn.addTextSyns(x)

        dAllSyns[mirexploreGeneID] = syn

vAllSyns = [dAllSyns[x] for x in dAllSyns]

for x in synIDCounter:
    if synIDCounter[x] > 1:
        print("Attention: double syn id! " + str(x) + " -> " + str(synIDCounter[x]))

globalKeywordExcludes = loadExludeWords()
vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=0, addAlphaBeta=True, removeSyn=lambda synonym: synonym.id.startswith('MIR') and not synonym.id.endswith('HG'))
printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/hgnc.syn", codec='utf8')
