from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID


def printChildren(node, maxlevel, curlevel, excluded):

    if maxlevel == 0:
        return

    if excluded != None and node.id in excluded:
        return

    children = node.getAllChildren(1)

    print("".join( ["-"]*curlevel + [node.name] ))

    for child in [x.term for x in children]:
        printChildren(child, maxlevel-1, curlevel+1, excluded)



def printTerms( oboLocation, selOboIDs, excludeIDs=None, printDepht = 3):

    ontology = GeneOntology(oboLocation)

    allExclIDs = set()

    if excludeIDs != None:
        for elem in excludeIDs:

            exclElem = ontology.dTerms[elem]
            allChildren = exclElem.getAllChildren()

            for child in allChildren:
                allExclIDs.add(child.term.id)


    for termID in ontology.dTerms:

        oboNode = ontology.dTerms[termID]

        oboID = oboNode.id
        oboName = oboNode.name

        oboSyns = oboNode.synonym
        oboRels = oboNode.is_a

        if oboID in selOboIDs:

            print("Base", oboID, oboName)
            printChildren(oboNode, 5, 1, allExclIDs)


print("Homeostasis")
print()
printTerms(dataDir + "miRExplore/go/go.obo",  ['GO:0042592'], printDepht=5)

print()
print()
print()

print()

# print("Neutrophils")
# print()
# printTerms(dataDir + "miRExplore/foundational_model_anatomy/fma_obo.obo",  ['FMA:62860'])
#
# print()
# print()
# print()
#
# print()
# print("Tissues")
# print()
# printTerms(dataDir + "miRExplore/foundational_model_anatomy/fma_obo.obo",  ['FMA:67498', 'FMA:9637', 'FMA:68646'])
#
#
# print()
# print()
# print()
#
# print()
# print("Acute microbial infections")
# print()
# printTerms(dataDir + "miRExplore/doid.obo",  ['DOID:0050117'], printDepht=5)

print()
print()
print()

print()
print("Chronic sterile inflammatory responses")
print()
printTerms(dataDir + "miRExplore/doid.obo",  ['DOID:4'], excludeIDs=['DOID:0050117', 'DOID:225', 'DOID:150'], printDepht=5)


#printTerms(dataDir + "miRExplore/go/go.obo",  ['GO:0002544', 'GO:0002676'], printDepht=3)
#printTerms(dataDir + "miRExplore/doid.obo",  ['DOID:104'], printDepht=5)

