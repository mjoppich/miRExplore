import shlex
from collections import Counter

from porestat.utils import DataFrame

from synonymes.Synonym import Synonym
from utils.idutils import printToFile, dataDir

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

    #syn = Synonyme(hgncID)
    syn = Synonym( mirexploreGeneID )
    syn.addTextSyns(result[hgncSymIdx])
    syn.addTextSyns(result[hgncSymIdx])
    syn.addTextSyns(result[hgncNameIdx])
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

synCounter = Counter()
for synonym in vAllSyns:

    for syn in synonym:
        synCounter[syn] += 1

setCommonWords = set()

for synWordCount in synCounter.most_common(66):

    setCommonWords.add(synWordCount[0])
    print(synWordCount[0] + " " + str(synWordCount[1]))

vPrintSyns = []

for synonym in vAllSyns:

    synonym.removeCommonSynonymes(setCommonWords)

    if len(synonym) > 0:
        vPrintSyns.append(synonym)

printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonymes/hgnc.syn")
