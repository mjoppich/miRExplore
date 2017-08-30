from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

celloObo = GeneOntology(dataDir + "miRExplore/cellosaurus/cellosaurus.obo")
tax2cells = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name

    oboSyns = oboNode.synonym
    oboXRefs = oboNode.xref

    taxID = {'all'}
    if oboXRefs != None:
        for xref in oboXRefs:
            if xref.startswith('NCBI_TaxID'):

                newTaxID = xref.split(' ')[0].split(':')[1]

                if newTaxID in allowedTaxIDs:
                    taxID.add(newTaxID)

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)

    if oboSyns != None:
        for x in oboSyns:
            newSyn.addSyn(x.syn)

    for taxid in taxID:
        tax2cells[taxid].add(newSyn)

    #print(str(taxID) + " " + str(newSyn))

globalKeywordExcludes = loadExludeWords()

for taxid in tax2cells:
    taxSyns = tax2cells[taxid]

    vPrintSyns = handleCommonExcludeWords(taxSyns, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=5)
    printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/celllines."+taxid+".syn")