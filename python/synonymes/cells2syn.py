import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import defaultdict
from synonymes.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


from collections import defaultdict
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


if __name__ == '__main__':

    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-o', '--obo', type=argparse.FileType("r"), required=True, help="input ontology file")
    parser.add_argument('-s', '--syn', type=argparse.FileType("w"), required=True, help="output synonym file")
    args = parser.parse_args()

    celloObo = GeneOntology(args.obo.name)
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
    printToFile(vPrintSyns, args.syn.name)