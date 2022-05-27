from collections import defaultdict
import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from synonymes.GeneOntology import GeneOntology
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


if __name__ == '__main__':

    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-o', '--obo', type=argparse.FileType("r"), required=True, help="input ontology file")
    parser.add_argument('-s', '--syn', type=str, required=True, help="output synonym FOLDER")
    args = parser.parse_args()

    celloObo = GeneOntology(args.obo.name)
    namespace2syn = defaultdict(set)

    allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

    ignoreTerms = set()
    ignoreTerms.add("GO:0005574") #obsolete DNA/DNA

    print("Total terms:", len(celloObo.dTerms), "Ignore terms", len(ignoreTerms))

    for cellID in celloObo.dTerms:

        oboNode = celloObo.dTerms[cellID]

        if oboNode.id in ignoreTerms:
            continue

        if len(oboNode.namespace)==0:
            print("has no namespace: " + oboNode.id)
            continue

        for ns in oboNode.namespace:
            namespace2syn[ns].add(oboNode)


    globalKeywordExcludes = loadExludeWords(common=False, cell_co=False, disease=False, generic=False)

    for x in globalKeywordExcludes:
        if 'membrane' in globalKeywordExcludes[x]:
            print("Membrane: " + x)

    for namespace in namespace2syn:

        allNodes = namespace2syn[namespace]
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
        #printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/go."+namespace.replace(' ', '_') + ".syn")
        #printToFile(vPrintSyns, dataDir + "/hpyloriDB/tm/go."+namespace.replace(' ', '_') + ".syn")
        printToFile(vPrintSyns, "{}/go.{}.syn".format(args.syn, namespace.replace(' ', '_')))