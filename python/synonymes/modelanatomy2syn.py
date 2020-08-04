import os,sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")



from collections import defaultdict
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID
from synonymes.GeneOntology import GeneOntology

bodypartsObo = GeneOntology(dataDir + "miRExplore/foundational_model_anatomy/fma_obo.obo")
vAllSyns = []

for cellID in bodypartsObo.dTerms:

    oboNode = bodypartsObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name

    oboSyns = oboNode.synonym
    oboRels = oboNode.is_a

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)


    aName = oboName.split(' ')

    if len(aName) > 1 and len(aName) < 5:

        acro = ""
        if aName[-1].upper() == 'CELL':
            acro = "".join([x[0].upper() for x in aName])

        newSyn.addSyn(acro)

    if oboSyns != None:
        for x in oboSyns:
            newSyn.addSyn(x.syn)

    #print(str(taxID) + " " + str(newSyn))

    vAllSyns.append(newSyn)

globalKeywordExcludes = loadExludeWords(cell_co=False)

vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=200, maxCommonCount=5)
printToFile(vPrintSyns, "/mnt/d/dev/data/pmid_jun2020/synonyms/model_anatomy.syn")