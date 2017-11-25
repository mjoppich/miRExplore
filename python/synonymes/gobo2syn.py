from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology
from owlready2 import namespace

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

celloObo = GeneOntology(dataDir + "miRExplore/go/go.obo")
namespace2syn = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    if len(oboNode.namespace)==0:
        print("has no namespace: " + oboNode.id)
        continue

    for ns in oboNode.namespace:
        namespace2syn[ns].add(oboNode)


globalKeywordExcludes = loadExludeWords()

for namespace in namespace2syn:

    allNodes = namespace2syn[namespace]
    synSet = set()

    for node in allNodes:
        newSyn = Synonym(node.id)
        newSyn.addSyn(node.name)

        if node.synonym != None:
            for x in node.synonym:
                newSyn.addSyn(x.syn)

        synSet.add(newSyn)


    vPrintSyns = handleCommonExcludeWords(synSet, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=5)
    printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/go."+namespace.replace(' ', '_') + ".syn")