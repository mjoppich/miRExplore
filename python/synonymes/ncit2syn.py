from synonymes.GeneOntology import GeneOntology
from textdb.NcitTerm2Symbols import NcitTermSymbolDB
from collections import defaultdict

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

ncitObo = GeneOntology(dataDir + "miRExplore/obodir/ncit.obo")
ncitTerm2Sym = NcitTermSymbolDB.loadFromFolder()

vAllSyns = []

for termID in ncitObo.dTerms:

    oboNode = ncitObo.dTerms[termID]

    oboID = oboNode.id
    oboName = oboNode.name

    oboSyns = oboNode.synonym
    oboRels = oboNode.is_a

    newSyn = Synonym(oboID)
    newSyn.addSyn(oboName)

    if oboSyns != None:
        for x in oboSyns:
            newSyn.addSyn(x.syn)


    allOrgs = [x for x in ncitTerm2Sym.org_term2symbol]

    for org in allOrgs:

        ncitID = oboID[oboID.index(":")+1:]

        if ncitID in ncitTerm2Sym.org_term2symbol[org]:

            orgSyms = ncitTerm2Sym.org_term2symbol[org][ncitID]

            for sym in orgSyms:
                newSyn.addSyn(sym)

    vAllSyns.append(newSyn)

globalKeywordExcludes = loadExludeWords()

vPrintSyns = handleCommonExcludeWords(vAllSyns, None, mostCommonCount=100, maxCommonCount=5) #globalKeywordExcludes
#printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/ncit.syn")
printToFile(vPrintSyns, "/mnt/d/dev/data/pmid_jun2020/synonyms/ncit.syn", codec='utf8')


