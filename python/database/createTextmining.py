import glob
import sys
import os

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import Counter, defaultdict

import nltk

from porestat.utils.DataFrame import DataFrame
import re


from database.ORGMIRs import ORGMIRDB
from synonymes.SynfileMap import SynfileMap
from synonymes.SynonymFile import Synfile
from synonymes.mirnaID import miRNA, miRNAPART
from textmining.SentenceDB import SentenceDB, RegPos
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import ltype2label, makeDBGeneID, mirtarbase_exp_type, mirtarbase_function_label, speciesName2TaxID, \
    dataDir
from database.Neo4JInterface import neo4jInterface
from utils.parallel import MapReduce
from enum import Enum

resultBase = dataDir + "/miRExplore/textmine/results/"
mirnaSyns = SynfileMap(resultBase + "/mirna/synfile.map")
mirnaSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

hgncSyns = SynfileMap(resultBase + "/hgnc/synfile.map")
hgncSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

db = None

if False:

    db = neo4jInterface(simulate=False)
    db.deleteRelationship('n', ['GENE'], None, 'm', ['PUBMED'], None, ['ST_MENTION'], None, 'r')

    db.deleteRelationship('n', ['PUBMED_AUTHOR'], None, 'm', ['PUBMED'], None, ['IS_AUTHOR'], None, 'r')
    db.deleteRelationship('n', ['DISEASE'], None, 'm', ['PUBMED'], None, ['DISEASE_MENTION'], None, 'r')
    db.deleteRelationship('n', ['PUBMED'], None, 'm', ['MIRTARBASE'], None, ['MIRTARBASE_LITERATURE_SUPPORT'], None, 'r')

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

        self.sameSentence = False
        self.sameParagraph = False
        self.relation = None
        self.mirnaFound = None

    def __str__(self):
        return "{pub}\t{type}\t{geneid}\t{mirnaid}\t{mirnaname}".format(pub=self.pubmed, type=self.idtype, geneid=self.gene, mirnaid=self.mirna, mirnaname=self.mirnadesc)

    def __repr__(self):
        return self.__str__()

    def getIdTuple(self):
        return (self.gene, self.mirna, self.idtype)


neutralAssoc = {'modulate': [RegPos.BEFORE, RegPos.AFTER, RegPos.BETWEEN],
                'involve': [RegPos.BEFORE, RegPos.AFTER],
                'inhibit': [RegPos.BEFORE, RegPos.AFTER]
                }
posRegulation = {'induce': [RegPos.BETWEEN],
                 'promote': [RegPos.BETWEEN],
                 'enhance': [RegPos.BETWEEN],
                 'up-regulate': [RegPos.BETWEEN],
                 'increase': [RegPos.BETWEEN]

                 }

negRegulation = {'reduce': [RegPos.BETWEEN],
                 'inhibit': [RegPos.BETWEEN],
                 'enhance': [RegPos.BETWEEN],
                 'down-regulate': [RegPos.BETWEEN],
                 'decrease': [RegPos.BETWEEN]

                 }

negation = {'don\'t': [RegPos.BEFORE, RegPos.BETWEEN, RegPos.AFTER],
            'doesn\'t': [RegPos.BEFORE, RegPos.BETWEEN, RegPos.AFTER],
            'didn\'t': [RegPos.BEFORE, RegPos.BETWEEN, RegPos.AFTER],
            'couldn\'t': [RegPos.BEFORE, RegPos.BETWEEN, RegPos.AFTER],
            'not': [RegPos.BEFORE, RegPos.BETWEEN, RegPos.AFTER]
            }

def findAssocs( assocs, text, textLoc):

    res = []
    for word in assocs:

        allowedAssocPos = assocs[word]

        if False and not textLoc in allowedAssocPos:
            continue

        if word in text:
            res.append( (word, textLoc) )

    return res



def findRelation(mirnaHit, hgncHit, sentDB):

    sentence = sentDB.get_sentence( mirnaHit.documentID )

    (textBefore, textBetween, textAfter, hitOrder) = sentence.extract_text(mirnaHit.position, hgncHit.position)

    tokens = nltk.RegexpTokenizer(r'[^\s]+').tokenize(sentence.text.lower()) #[mir\-|\w]+
    tokenPos = [x for x in nltk.RegexpTokenizer(r'[^\s]+').span_tokenize(sentence.text.lower())]

    text = nltk.Text(tokens)
    tags = nltk.pos_tag(text)

    test = []
    for idx, word in enumerate(tags):

        wordPos = tokenPos[idx]

        accept = False
        if word[1][0:2] in ['VB', 'NN', 'JJ']:
            accept = True

        """
        if word[1][0:2] in ['JJ']:

            if idx+1 < len(tags):
                if tags[idx+1][1][0:2] in ['VB', 'NN']:
                    accept = True

            beforeCCNN = True
            for i in range(0, 3):
                if idx+i < len(tags):
                    if not tags[idx+i][1][0:2] in ['JJ', 'CC', 'NN']:
                        beforeCCNN = False

            accept = beforeCCNN
        """

        if accept:
            test.append( (word[0], word[1], wordPos[0], wordPos[1]) )

    def findTagWithPos( alltags, pos):

        startIdx = None
        endIdx = None

        for idx,tagpos in enumerate(alltags):

            if tagpos[2] <= pos[0] and pos[0] <= tagpos[3] and startIdx == None:
                startIdx = idx
            if tagpos[2] <= pos[1] and pos[1] <= tagpos[3]:
                endIdx = idx


        if startIdx == None or endIdx == None or not startIdx <= endIdx:
            return None

        return (startIdx, endIdx)


    idxM = findTagWithPos(test, mirnaHit.position)
    idxG = findTagWithPos(test, hgncHit.position)


    # cases to consider

    ## 1: find verb between two elements => take that
    ## 2: no verb between two elements
    ##### take closest verb to mirna

    def findInteraction(alltags, startI, endI):

        foundInteractions = []
        for i in range(min(startI, endI), max(startI, endI)):

            tagpos = alltags[i]

            tagPrev = None if i-1 < 0 else alltags[i-1]
            tagNext = None if i+1 >= len(alltags) else alltags[i+1]

            if tagpos[1][0:2] == 'VB':

                for x in negRegulation:
                    if x in tagpos[0]:
                        foundInteractions.append( (x, 'NEG') )

                for x in posRegulation:
                    if x in tagpos[0]:
                        foundInteractions.append( (x, 'POS') )

                for x in neutralAssoc:
                    if x in tagpos[0]:
                        foundInteractions.append( (x, 'NEU') )


        return foundInteractions

    if idxM == None or idxG == None:
        #print(idxM, idxG)
        return None


    if idxM[0] < idxG[0]:
        # type MG
        betweeenInteractions = findInteraction(test, idxM[1], idxG[0])
    else:
        betweeenInteractions = findInteraction(test, idxM[0], idxG[1])


    if betweeenInteractions != None and len(betweeenInteractions) > 0:
        return betweeenInteractions

    return None


def findCooccurrences( pubmed, hgncHits, mirnaHits, sentDB ):

    def checkSynHit(synhit):
        if len(synhit.foundSyn) <= 5:
            return synhit.perfectHit == True

        return True

    def chekSynHitMirna(synhit):

        if len(synhit.foundSyn) <= 5:
            foundSyn = synhit.foundSyn.lower()
            return foundSyn.startswith('mir') or foundSyn.startswith('micro')

        return True


    setAllGenes = set([x for x in hgncHits if checkSynHit(x)])
    setAllMirnas = set([x for x in mirnaHits if chekSynHitMirna(x)])

    hgncBySent = defaultdict(list)
    mirnaBySent = defaultdict(list)

    hgncToSent = {}
    mirnaToSent = {}

    for hit in hgncHits:
        parSenID = (hit.documentID.parID, hit.documentID.senID)
        hgncBySent[parSenID].append(hit)

        hgncToSent[ hit ] = parSenID

    for hit in mirnaHits:
        parSenID = (hit.documentID.parID, hit.documentID.senID)
        mirnaBySent[parSenID].append(hit)

        mirnaToSent[ hit ] = parSenID


    allCoocs = []

    for x in setAllMirnas:
        for y in setAllGenes:

            foundCooc = Cooccurrence()
            foundCooc.pubmed = pubmed

            if re.match('MIPF[0-9]+', x.synonym.id) != None:
                foundCooc.idtype="MIRNA_FAMILY"
            elif re.match('MIMAT[0-9]+', x.synonym.id) != None:
                foundCooc.idtype="MIRNA"
            elif re.match('MI[0-9]+', x.synonym.id) != None:
                foundCooc.idtype='MIRNA_PRE'
            elif re.match('ORGMIR[0-9]+', x.synonym.id) != None:
                foundCooc.idtype='MIRNA_ORGMIR'
            elif re.match('ORGMI[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = 'MIRNA_ORGMIR'
            else:
                foundCooc.idtype='UNKNOWN'

            foundCooc.gene = y.synonym.id
            foundCooc.mirna = x.synonym.id
            foundCooc.mirnadesc=str(x.synonym)

            foundCooc.mirnaFound = x.hitSyn

            for mirnaSyn in x.synonym.syns:

                if mirnaSyn.startswith("miR-") and not 'mediated' in mirnaSyn:
                    test = miRNA(mirnaSyn)
                    outstr = test.getStringFromParts([miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM])
                    foundCooc.mirnaFound = outstr
                    break

            if True and pubmed == '21752897':
                print(pubmed)


            miRNALoc = mirnaToSent[x]
            hgncLoc = hgncToSent[y]

            if miRNALoc[0] == hgncLoc[0]:
                foundCooc.sameParagraph = True

                if miRNALoc[1] == hgncLoc[1]:
                    foundCooc.sameSentence = True

                    foundCooc.relation = findRelation(x,y, sentDB)

            allCoocs.append(foundCooc)


    coocinfos = defaultdict(list)
    for cooc in allCoocs:
        if cooc.relation == None and not cooc.sameSentence:
            continue

        coocinfos[(cooc.gene, cooc.mirnaFound, cooc.mirna)].append((cooc.pubmed, cooc.relation))

    if len(coocinfos) > 0:
        for gene,mirna, mirnaid in coocinfos:
            print("{gene}\t{mirna}\t{mirnaid}\t{info}".format(gene=gene, mirna=mirna, mirnaid=mirnaid, info=coocinfos[(gene, mirna, mirnaid)]))

    return allCoocs

coocCounter = Counter()
idTuple2Pubmed = defaultdict(set)
orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")


allfiles = glob.glob(resultBase + "/hgnc/pubmed18n*.index")
allfileIDs = [int(os.path.basename(x).replace('pubmed18n', '').replace('.index','')) for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

#allfileIDs = [700]

def analyseFile(splitFileIDs, env):


    for splitFileID in splitFileIDs:

        fileID = "{:>4}".format(splitFileID).replace(" ", "0")

        hgncFile = resultBase + "/hgnc/pubmed18n"+fileID+".index"
        mirnaFile = resultBase + "/mirna/pubmed18n"+fileID+".index"
        sentFile = "/mnt/c/dev/data/pubmed/pubmed18n" + fileID + ".sent"

        mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)
        if len(mirnaHits) == 0:
            return

        hgncHits = SyngrepHitFile(hgncFile, hgncSyns)
        if len(hgncHits) == 0:
            return

        sentDB = SentenceDB(sentFile)

        sys.stderr.write("Found something in: " + str(fileID) + "\n")

        for docID in mirnaHits:

            if docID in hgncHits:

                mirnaSynHits = mirnaHits.getHitsForDocument(docID)
                hgncSynHits = hgncHits.getHitsForDocument(docID)

                #if docID == 'a27229723':
                #    [print(x.synonyme) for x in hgncSynHits]
                #    [print(x.synonyme) for x in mirnaSynHits]

                foundCoocs = findCooccurrences(str(docID), hgncSynHits, mirnaSynHits, sentDB)

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

                    if db:
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

                        edgeMirnas = [x[0] for x in geneEdges[edge]]

                        if db:
                            db.createRelationship('gene', [edge[1]], {'id': edge[0]}, 'pmid', ['PUBMED'], {'id': docID},
                                              ['ST_MENTION'], {'type': 'GENE_MENTION', 'mirnas': edgeMirnas})

                    for edge in mirnaEdges:

                        edgeGenes = [x[0] for x in mirnaEdges[edge]]
                        if db:
                            db.createRelationship('pmid', ['PUBMED'], {'id': docID}, 'mi', [edge[1]], {'id': edge[0]},
                                              ['ST_MENTION'], {'type': 'MIRNA_MENTION', 'genes': edgeGenes})


ll = MapReduce(4)
result = ll.exec( allfileIDs, analyseFile, None, 1, None)

for idTuple in coocCounter:

    cnt = coocCounter[idTuple]
    print(idTuple[0], idTuple[1], str(idTuple[2]), str(cnt), ",".join(idTuple2Pubmed[idTuple]))

    #print(str(idTuple) + " --> " + str(cnt) + " in " + str(idTuple2Pubmed[idTuple]))



