import csv
from io import StringIO

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, printToFile, loadExludeWords

globalKeywordExcludes = loadExludeWords()

csv.register_dialect('phenotypes', delimiter=',', doublequote=True, escapechar='', lineterminator='\n')
with open(dataDir + "miRExplore/pharmgkb/phenotypes.tsv", 'r') as infile:
    infile.readline()

    vAllSyns = []
    idCnt = 0

    for line in infile:
        aline = [x.strip() for x in line.split('\t')]

        name = aline[1]

        altNames = StringIO()
        altNames.write(aline[2] + "\n")

        names = []

        for line in csv.reader([aline[2]], dialect='phenotypes'):
            for elem in line:
                names.append(elem)

        names = names + [name]

        newSyn = Synonym( 'DISEASE' + str(len(vAllSyns)+1))
        newSyn.addSyn(name)

        for x in names:

            if x.startswith('[D]') or x.startswith('[X]') or x.startswith('[M]'):
                x = x[3:]

            xsyns = []

            if ' - ' in x:
                xsyns += x.split(' - ')
            else:
                xsyns.append(x)

            for xsyn in xsyns:
                newSyn.addSyn(xsyn.strip())

        vAllSyns.append(newSyn)

    vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=10)
    printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/disease.syn")