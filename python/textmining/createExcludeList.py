from utils.idutils import aminoAcids, dataDir, printToFile
from nertoolkit.geneontology.GeneOntology import GeneOntology


addGeneric = True
addNames = True
addGO = True

if addGeneric:
    globalKeywordExcludes = ['MIR'] + ['MIR' + chr(x) for x in range(ord('A'), ord('Z') + 1)] + ['MIR-' + chr(x) for x
                                                                                                 in range(ord('A'),
                                                                                                          ord('Z') + 1)]
    globalKeywordExcludes += aminoAcids + [x.upper() for x in aminoAcids]
    globalKeywordExcludes += ['TYPE II', 'type II', 'type I', 'TYPE I', ] + ['type ' + str(x) for x in range(0, 10)] + [
        'TYPE ' + str(x) for x in range(0, 10)]

    globalKeywordExcludes = sorted(globalKeywordExcludes)
    printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.generic.syn")

if addNames:
    with open(dataDir + "miRExplore/textmine/excludes/names.dmp", 'r') as infile:

        globalKeywordExcludes = set()

        for line in infile:
            aline = [x.strip() for x in line.split('|')]

            taxName = aline[1]

            if taxName[0] == '\"':
                closing = taxName.find("\"", 1)
                taxName = taxName[1:closing]

            if taxName[0] == '\'':
                closing = taxName.find("\'", 1)
                taxName = taxName[1:closing]

            taxName = taxName.strip()

            aTaxName = taxName.split(" ")
            taxName = " ".join(aTaxName[0:min(2, len(aTaxName))])

            try:
                int(taxName)
            except ValueError:
                globalKeywordExcludes.add(taxName)

        globalKeywordExcludes = sorted(globalKeywordExcludes)
        printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.names.syn")


if addGO:

    goObo = GeneOntology(dataDir + "miRExplore/textmine/excludes/go.obo")

    globalKeywordExcludes = set()

    for id in goObo.dTerms:

        child = goObo.dTerms[id]
        if 'cellular_component' in child.namespace:

            globalKeywordExcludes.add(child.name)

            if child.synonym != None:
                for syn in child.synonym:
                    if syn != None:
                        globalKeywordExcludes.add(syn.syn)

    globalKeywordExcludes = sorted(globalKeywordExcludes)
    printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.cell_co.syn")
