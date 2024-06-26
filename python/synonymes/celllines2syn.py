from collections import defaultdict
import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from synonymes.GeneOntology import GeneOntology
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

celloObo = GeneOntology(dataDir + "miRExplore/meta_cells.obo")

ignoreTerms = []
for ignoreTerm in ["META:3", "META:707", "EFO:0000408", "META:403", "EFO:0000408", "BFO:0000017"]:
    ignoreTerms += [x.termid for x in celloObo[ignoreTerm].getAllChildren()]

ignoreTerms = set(ignoreTerms)

ignoreTerms.add("META:340")
ignoreTerms.add("META:10")

for term in celloObo.dTerms:
    oboNode = celloObo.dTerms[term]
    if oboNode.id.startswith(("CHEBI", "PR", "DOID", "UBERON", "IAO")):
        ignoreTerms.add(oboNode.id)

print("Total terms:", len(celloObo.dTerms), "Ignore terms", len(ignoreTerms))

tax2cells = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name

    if oboID.startswith('GO'):
        continue

    if oboID in ignoreTerms:
        continue

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

    if "-" in oboName:
        newSyn.addSyn(oboName.replace("-", ""))

    isHela = False

    if oboSyns != None:
        for x in oboSyns:

            if x == None:
                print(oboNode)

            newSyn.addSyn(x.syn)

            if "-" in x.syn:

                newSyn.addSyn(x.syn.replace("-", ""))


    if oboNode.id.startswith('CL:'):

        toadd = set()
        for syn in newSyn.syns:
            asyn = syn.split(' ')
            if asyn[-1].upper() == 'CELL' and len(asyn) < 4:
                short = "".join([x[0].upper() for x in asyn])
                toadd.add(short)

        if len(toadd) > 0:

            for x in toadd:
                newSyn.addSyn(x)

    """
    if oboNode.name.split(" ")[-1].upper() == 'CELL':
        allsyn = [x for x in newSyn.syns]

        for x in allsyn:
            newSyn.addSyn(x + "s")
    """
    allsyn = [x for x in newSyn.syns]
    for x in allsyn:
        if x[-1] != 's' and not x[-1].isdigit():
            newSyn.addSyn(x + "s")

    if isHela:
        print(newSyn)

    if oboID == "CVCL_B478":
        print(newSyn)

    for taxid in taxID:
        tax2cells[taxid].add(newSyn)

    #print(str(taxID) + " " + str(newSyn))

globalKeywordExcludes = loadExludeWords(common=False,cell_co=False)


for taxid in tax2cells:
    taxSyns = tax2cells[taxid]

    vPrintSyns = handleCommonExcludeWords(taxSyns, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=15)
    printToFile(vPrintSyns, "/mnt/d/dev/data/pmid_jul2020/synonyms/celllines.{}.syn".format(taxid)) #dataDir + "/miRExplore/textmine/synonyms/celllines."+taxid+".syn")

    print("Wrote", dataDir + "/miRExplore/textmine/synonyms/celllines."+taxid+".syn")