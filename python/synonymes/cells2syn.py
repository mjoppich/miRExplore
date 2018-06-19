from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

celloObo = GeneOntology(dataDir + "miRExplore/cell_ontology/cl.obo")
vAllSyns = []

allOboNames = defaultdict(set)

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id

    if not oboID.startswith("CL"):
        continue

    oboName = oboNode.name

    allOboNames[oboName].add(oboID)

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id

    if not oboID.startswith("CL"):
        continue

    if oboID == 'CL:1000413':
        print(oboID)
        print(oboNode.name)

    oboName = oboNode.name
    oboSyns = oboNode.synonym
    oboRels = oboNode.is_a

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)

    if oboSyns != None:
        for x in oboSyns:

            if x == None:
                continue

            if x.syn in allOboNames:
                continue

            newSyn.addSyn(x.syn)

    for x in newSyn.syns:

        if 'cell of' in x:
            coidx = x.index('cell of')

            newval = x[:coidx]
            suffix = x[coidx+8:]

            newval = suffix + " " + newval + "cell"

            if newval in allOboNames:
                print(x, newval, oboID, allOboNames[newval])

                if not oboID in allOboNames[newval]:
                    continue


            newSyn.addSyn(newval)


    vAllSyns.append(newSyn)

globalKeywordExcludes = loadExludeWords()

vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=1, maxCommonCount=5)
printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/cells.syn")