from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology
from owlready2 import namespace

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

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id

    if not oboID.startswith("CL"):
        continue

    oboName = oboNode.name

    oboSyns = oboNode.synonym
    oboRels = oboNode.is_a

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)

    if oboSyns != None:
        for x in oboSyns:

            if x == None:
                continue

            newSyn.addSyn(x.syn)

    for x in newSyn.syns:

        if 'cell of' in x:
            coidx = x.index('cell of')

            newval = x[:coidx]
            newval += 'cell'

            newSyn.addSyn(newval)


    vAllSyns.append(newSyn)

globalKeywordExcludes = loadExludeWords()

vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=1, maxCommonCount=5)
printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/cells.syn")