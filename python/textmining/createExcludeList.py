import csv
from io import StringIO

from utils.idutils import aminoAcids, dataDir, printToFile
from nertoolkit.geneontology.GeneOntology import GeneOntology


addGeneric = True
addNames = True
addGO = True
addDisease = True

if addDisease:

    csv.register_dialect('phenotypes',  delimiter=',', doublequote=True, escapechar='', lineterminator='\n')
    with open(dataDir + "miRExplore/pharmgkb/phenotypes.tsv", 'r') as infile:

        globalKeywordExcludes = set()
        infile.readline()

        for line in infile:
            aline = [x.strip() for x in line.split('\t')]

            name = aline[1]

            altNames = StringIO()
            altNames.write(aline[2] + "\n")

            names = []

            for line in csv.reader( [aline[2]], dialect='phenotypes'):
                for elem in line:
                    names.append(elem)

            names = names + [name]

            for x in names:

                if x.startswith('[D]') or x.startswith('[X]') or x.startswith('[M]'):
                    x = x[3:]

                globalKeywordExcludes.add(x.strip())

        globalKeywordExcludes = sorted(globalKeywordExcludes)
        printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.disease.syn")
        print("Done: disease")

if addGeneric:
    globalKeywordExcludes = ['MIR'] + ['MIR' + chr(x) for x in range(ord('A'), ord('Z') + 1)] + ['MIR-' + chr(x) for x
                                                                                                 in range(ord('A'),
                                                                                                          ord('Z') + 1)]
    globalKeywordExcludes += aminoAcids + [x.upper() for x in aminoAcids]
    globalKeywordExcludes += ['TYPE I', 'type II', 'type III' ] + ['type ' + str(x) for x in range(0, 10)] + [
        'TYPE ' + str(x) for x in range(0, 10)]

    globalKeywordExcludes += ['ligand', 'murine', 'high grade', 'fetal', 'class I', 'class II', 'class III', 'II', 'I', 'III', 'transcription factor', 'transmembrane', 'Na+', 'K+', 'trans']
    globalKeywordExcludes += ['etc.', 'smooth muscle', 'skeletal muscle', 'generalized', 'G protein-coupled', 'proto-oncogene']
    globalKeywordExcludes += ['membrane bound', 'early onset', 'atopic', 'short form', 'embryonic lethal', 'D-aspartate', 'C-terminal']

    with open(dataDir + "miRExplore/textmine/excludes/common_english_words.syn", 'r') as infile:

        globalKeywordExcludes = set(globalKeywordExcludes)

        for word in infile:
            globalKeywordExcludes.add(word)


    globalKeywordExcludes = sorted(globalKeywordExcludes)
    printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.generic.syn")
    print("Done: generic")

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
        print("Done: names")


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

                        synWord = syn.syn

                        globalKeywordExcludes.add(synWord)

                        if ' ' in synWord:
                            globalKeywordExcludes.add(synWord.split(' ')[0])

    globalKeywordExcludes = sorted(globalKeywordExcludes)
    printToFile(globalKeywordExcludes, dataDir + "miRExplore/textmine/excludes/exclude_words.cell_co.syn")
    print("Done: cell_co")
