from collections import Counter, defaultdict
from porestat.utils.DataFrame import DataFrame
import re
from database.ORGMIRs import ORGMIRDB
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

db = neo4jInterface(simulate=False)
db.deleteRelationship('n', ['GENE'], None, 'm', ['PUBMED'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRNA'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRNA_FAMILY'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRNA_ORGMI'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRNA_ORGMIR'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRNA_PRE'], None, ['ST_MENTION'], None, 'r')
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['IS_A_MIRNA'], None, ['ST_MENTION'], None, 'r')
db.deleteNode(["PUBMED"], None)

db.createUniquenessConstraint('PUBMED', 'id')

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

    def getIdTuple(self):
        return (self.gene, self.mirna, self.idtype)

def findCooccurrences( pubmed, hgncHits, mirnaHits ):

    def checkSynHit(synhit):
        if len(synhit.foundSyn) <= 5:
            return synhit.perfectHit == True

        return True


    setAllGenes = set([x.synonym for x in hgncHits if checkSynHit(x)])
    setAllMirnas = set([x.synonym for x in mirnaHits if checkSynHit(x)])

    foundCoocs = []

    for x in setAllMirnas:
        for y in setAllGenes:

            foundCooc = Cooccurrence()
            foundCooc.pubmed = pubmed

            if re.match('MIPF[0-9]+', x.id) != None:
                foundCooc.idtype="MIRNA_FAMILY"
            elif re.match('MIMAT[0-9]+', x.id) != None:
                foundCooc.idtype="MIRNA"
            elif re.match('MI[0-9]+', x.id) != None:
                foundCooc.idtype='MIRNA_PRE'
            elif re.match('ORGMIR[0-9]+', x.id) != None:
                foundCooc.idtype='MIRNA_ORGMIR'
            elif re.match('ORGMI[0-9]+', x.id) != None:
                foundCooc.idtype = 'MIRNA_ORGMIR'
            else:
                foundCooc.idtype='UNKNOWN'

            foundCooc.gene = y.id
            foundCooc.mirna = x.id
            foundCooc.mirnadesc=str(x)

            foundCoocs.append(foundCooc)

    return foundCoocs

coocCounter = Counter()
idTuple2Pubmed = defaultdict(set)
orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

for splitFileID in range(892, 0, -1):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    hgncFile = resultBase + "/hgnc/medline17n"+fileID+".index"
    mirnaFile = resultBase + "/mirna/medline17n"+fileID+".index"

    mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)
    if len(mirnaHits) == 0:
        continue

    hgncHits = SyngrepHitFile(hgncFile, hgncSyns)
    if len(hgncHits) == 0:
        continue

    print("Found something in: " + str(fileID))

    for docID in mirnaHits:

        if docID == '22419229':
            print(docID)

        if docID in hgncHits:

            mirnaSynHits = mirnaHits.getHitsForDocument(docID)
            hgncSynHits = hgncHits.getHitsForDocument(docID)

            #if docID == 'a27229723':
            #    [print(x.synonyme) for x in hgncSynHits]
            #    [print(x.synonyme) for x in mirnaSynHits]

            foundCoocs = findCooccurrences(str(docID), hgncSynHits, mirnaSynHits)

            assocByGene = defaultdict(set)
            for x in foundCoocs:

                geneID = x.gene
                geneLabel = 'GENE'
                mirnaID = x.mirna
                mirnaLabel = x.idtype

                assoc= (geneID, geneLabel, mirnaID, mirnaLabel)
                assocByGene[assoc[0]].add(assoc)

            addDocAsEvidence = False
            assocByTypeForGene = {}
            for gene in assocByGene:
                assocs = assocByGene[gene]

                mimatSet = set()
                miSet = set()
                orgmirSet = set()
                familySet = set()

                for assoc in assocs:
                    if assoc[3] == 'MIRNA':
                        mimatSet.add(assoc)
                    elif assoc[3] == 'MIRNA_PRE':
                        miSet.add(assoc)
                    elif assoc[3] == 'MIRNA_ORGMIR':
                        orgmirSet.add(assoc)
                    elif assoc[3] == 'MIRNA_FAMILY':
                        familySet.add(assoc)
                    else:
                        print("Unknown relation in doc: " + docID)
                        print(assoc)

                if len(mimatSet) > 0 or len(miSet) > 0 or len(orgmirSet) > 0 or len(familySet) > 0:

                    assocByTypeForGene[gene] = (mimatSet, miSet, orgmirSet, familySet)

                #filter assocs here such that if a taxid specific version was found, not the general version is added
            if len(assocByTypeForGene) > 0:

                db.createNodeIfNotExists(['EVIDENCE', 'PUBMED'], {'id': docID})
                print("Adding: " + str(docID) + ": " + str(assocByTypeForGene))
                #
                # TODO first create all unique edges
                # TODO add genes to mirna edges and mirnas to gene edges to keep track from where an edge originates
                # TODO edges should get weights = how many relations have been found
                #

                mirnaEdges = defaultdict(set)
                geneEdges = defaultdict(set)


                for gene in assocByTypeForGene:
                    assocsForGene = assocByTypeForGene[gene] # (mimatSet, miSet, orgmirSet, familySet)

                    for subSet in assocsForGene:
                        for cooc in subSet:

                            geneEdges[ (cooc[0], cooc[1]) ].add( (cooc[2], cooc[3]) )
                            mirnaEdges[(cooc[2], cooc[3])].add(  (cooc[0], cooc[1]) )

                for edge in geneEdges:

                    edgeMirnas = list(geneEdges[edge])
                    db.createRelationship('gene', [edge[1]], {'id': edge[0]}, 'pmid', ['PUBMED'], {'id': docID},
                                          ['ST_MENTION'], {'type': 'GENE_MENTION', 'mirnas': edgeMirnas})

                for edge in mirnaEdges:

                    edgeGenes = list(mirnaEdges[edge])
                    db.createRelationship('pmid', ['PUBMED'], {'id': docID}, 'mi', [edge[1]], {'id': edge[0]},
                                          ['ST_MENTION'], {'type': 'MIRNA_MENTION', 'genes': edgeGenes})



for idTuple in coocCounter:

    cnt = coocCounter[idTuple]
    print(idTuple[0], idTuple[1], str(idTuple[2]), str(cnt), ",".join(idTuple2Pubmed[idTuple]))

    #print(str(idTuple) + " --> " + str(cnt) + " in " + str(idTuple2Pubmed[idTuple]))



