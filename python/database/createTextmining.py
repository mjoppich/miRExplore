from collections import Counter
from porestat.utils.DataFrame import DataFrame

from synonymes.SynfileMap import SynfileMap
from synonymes.SynonymFile import Synfile
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import ltype2label, makeDBGeneID, mirtarbase_exp_type, mirtarbase_function_label, speciesName2TaxID, \
    dataDir
from database.Neo4JInterface import neo4jInterface


resultBase = dataDir + "/miRExplore/textmine/results/"
mirnaSyns = SynfileMap(resultBase + "/mirna/synfile.map")
mirnaSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

hgncSyns = SynfileMap(resultBase + "/hgnc/synfile.map")
hgncSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )


class Cooccurrence:

    def __init__(self):
        self.pubmed = None
        self.idtype = None
        self.gene = None
        self.mirna = None
        self.mirnadesc = None

    def __str__(self):
        return "{pub}\t{type}\t{geneid}\t{mirnaid}\t{mirnaname}".format(pub=self.pubmed, type=self.idtype, geneid=self.gene, mirnaid=self.mirna, mirnaname=self.mirnadesc)

    def __repr__(self):
        return self.__str__()

def findCooccurrences( pubmed, hgncHits, mirnaHits ):

    def checkSynHit(synhit):
        if len(synhit.foundSyn) <= 5:
            return synhit.perfectHit == True

        return True


    setAllGenes = set([x.synonyme for x in hgncHits if checkSynHit(x)])
    setAllMirnas = set([x.synonyme for x in mirnaHits if checkSynHit(x)])

    foundCoocs = []

    for x in setAllMirnas:
        for y in setAllGenes:

            foundCooc = Cooccurrence()
            foundCooc.pubmed = pubmed

            if x.id.startswith("MIPF"):
                foundCooc.idtype="MIRNA_FAMILY"
            elif x.id.startswith('MIMAT'):
                foundCooc.idtype="MIRNA"
            elif x.id.startswith('MI:'):
                foundCooc.idtype='MIRNA_PRE'
            elif x.id.startswith('ORGMIR'):
                foundCooc.idtype='MIRNA_ORG'
            else:
                foundCooc.idtype='UNKNOWN'

            foundCooc.gene = y.id
            foundCooc.mirna = x.id
            foundCooc.mirnadesc=str(x)

            foundCoocs.append(foundCooc)


    return foundCoocs

for splitFileID in range(892, 0, -1):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    hgncFile = resultBase + "/hgnc/medline17n"+fileID+".index"
    mirnaFile = resultBase + "/mirna/medline17n"+fileID+".index"

    hgncHits = SyngrepHitFile(hgncFile, hgncSyns)
    mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)

    if len(mirnaHits) == 0 or len(hgncHits) == 0:
        continue

    for docID in mirnaHits:


        if docID in hgncHits:

            mirnaSynHits = mirnaHits.getHitsForDocument(docID)
            hgncSynHits = hgncHits.getHitsForDocument(docID)

            foundCoocs = findCooccurrences(str(docID), hgncSynHits, mirnaSynHits)

            for x in foundCoocs:
                print(x)



