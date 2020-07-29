from collections import defaultdict

import sys, os
sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/miRExplore/python/")))


from synonymes.GeneOntology import GeneOntology
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


filepath = sys.argv[1]


fileObo = GeneOntology(filepath)
namespace2syn = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])


allNodes = []
for cellID in fileObo.dTerms:
    oboNode = fileObo.dTerms[cellID]

    allNodes.append(oboNode)

globalKeywordExcludes = loadExludeWords(common=False, cell_co=False, disease=False, generic=False)

for x in globalKeywordExcludes:
    if 'membrane' in globalKeywordExcludes[x]:
        print("Membrane: " + x)

synSet = set()

for node in allNodes:
    newSyn = Synonym(node.id)
    newSyn.addSyn(node.name)

    if node.synonym != None:
        for x in node.synonym:
            if x == None:
                continue
            newSyn.addSyn(x.syn)

    synSet.add(newSyn)


vPrintSyns = handleCommonExcludeWords(synSet, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=5, minSynCount=0)
printToFile(vPrintSyns, sys.argv[2])