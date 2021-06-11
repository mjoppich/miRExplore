import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")
sys.path.insert(0, "/mnt/f/dev/git/NERtoolkit/")


from collections import defaultdict
from synonymes.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


infile = sys.argv[1] # dataDir + "miRExplore/doid.obo"
outfile = sys.argv[2] #"/mnt/d/dev/data/pmid_jun2020/synonyms/disease.syn"

celloObo = GeneOntology(infile)

ignoreTerms = set()

ignoreTerms.add("DOID:4")
print("Total terms:", len(celloObo.dTerms), "Ignore terms", len(ignoreTerms))


vAllSyns = []

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name

    if oboID in ignoreTerms:
        continue

    if oboNode.is_obsolete:
        print("skipping", oboName)
        continue

    oboSyns = oboNode.synonym
    oboRels = oboNode.is_a

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)

    if oboSyns != None:
        for x in oboSyns:
            newSyn.addSyn(x.syn)

    #print(str(taxID) + " " + str(newSyn))

    vAllSyns.append(newSyn)

globalKeywordExcludes = loadExludeWords()

vPrintSyns = handleCommonExcludeWords(vAllSyns, None, mostCommonCount=10, maxCommonCount=5)


printToFile(vPrintSyns, outfile, codec='utf8')