import codecs

from utils.idutils import dataDir

resultBase = dataDir + "/miRExplore/textmine/results/"

with codecs.open(resultBase + 'hgnc/medline17n0885.index', 'rb' ) as infile:

    for line in infile:
        print(line.strip())
        print(line.decode('latin1'))