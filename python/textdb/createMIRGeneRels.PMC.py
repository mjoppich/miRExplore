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
from synonymes.SynonymFile import Synfile, AssocSynfile
from synonymes.mirnaID import miRNA, miRNAPART
from textmining.SentenceDB import SentenceDB, RegPos
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import ltype2label, makeDBGeneID, mirtarbase_exp_type, mirtarbase_function_label, speciesName2TaxID, \
    dataDir
from database.Neo4JInterface import neo4jInterface
from utils.parallel import MapReduce
from enum import Enum

resultBase = dataDir + "/miRExplore/textmine/results_pmc/"
mirnaSyns = SynfileMap(resultBase + "/mirna/synfile.map")
mirnaSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

hgncSyns = SynfileMap(resultBase + "/hgnc/synfile.map")
hgncSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

relSyns = SynfileMap(resultBase + "/relations/synfile.map")
relSyns.loadSynFiles(('/home/users/joppich/ownCloud/data/', dataDir) )

relationSyns = AssocSynfile('/mnt/c/ownCloud/data/miRExplore/relations/allrels.csv')

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

def findAssocs( assocs, text, textLoc):

    res = []
    for word in assocs:

        allowedAssocPos = assocs[word]

        if False and not textLoc in allowedAssocPos:
            continue

        if word in text:
            res.append( (word, textLoc) )

    return res



def findRelation(mirnaHit, hgncHit, sentDB, relHits):
    """

    same sentence can be assumed!

    :param mirnaHit:
    :param hgncHit:
    :param sentDB:
    :param relHits:
    :return:
    """

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


    if idxM == None or idxG == None:
        #print(idxM, idxG)
        return None


    def findInteraction(a,b,c):

        return None


    relSyns = relHits.getHitsForDocument(mirnaHit.documentID.docID)
    relSyns = [x for x in relSyns if x.documentID == mirnaHit.documentID]


    if idxM[0] < idxG[0]:
        # type MG
        betweeenInteractions = findInteraction(test, idxM[1], idxG[0])
    else:
        #type GM
        betweeenInteractions = findInteraction(test, idxM[0], idxG[1])


    if betweeenInteractions != None and len(betweeenInteractions) > 0:
        return betweeenInteractions

    return None


def findRelationBySyns(mirnaHit, hgncHit, sentDB, relHits):

    sentence = sentDB.get_sentence( mirnaHit.documentID )
    (textBefore, textBetween, textAfter, hitOrder) = sentence.extract_text(mirnaHit.position, hgncHit.position)


    detRelations = relHits.getHitsForDocument(mirnaHit.documentID.docID)

    if detRelations == None or len(detRelations) == 0:
        return None

    sentHits = []
    for rel in detRelations:
        if rel.documentID == mirnaHit.documentID:
            sentHits.append(rel)


    negatedSentence = any([x in sentence.text for x in ['not', 'n\'t', 'nega']])

    allRelations = []

    for rel in sentHits:
        assocDir = None
        assocType = None
        assocWord = None
        assocSent = None
        assocDirRel = None

        if mirnaHit.position[0] < hgncHit.position[0]:
            assocDir = 'MG'

            if rel.position[0] < mirnaHit.position[0]:
                assocDirRel = 'V'+assocDir
            elif hgncHit.position[0] < rel.position[0]:
                assocDirRel = assocDir + 'V'
            else:
                assocDirRel = 'MVG'

        else:
            assocDir = 'GM'

            if rel.position[0] < hgncHit.position[0]:
                assocDirRel = 'V'+assocDir
            elif mirnaHit.position[0] < rel.position[0]:
                assocDirRel = assocDir + 'V'
            else:
                assocDirRel = 'GVM'

        assocSent = str(mirnaHit.documentID)
        assocWord = rel.synonym.id

        assocType = relationSyns.synid2class[rel.synonym.id]

        allRelations.append(
            (
                assocDir,
                assocDirRel,
                assocType,
                assocWord,
                assocSent,
                negatedSentence,
                mirnaHit.position,
                hgncHit.position,
                rel.position
            )
        )


    return allRelations




def findCooccurrences( pubmed, hgncHits, mirnaHits, sentDB, relHits ):

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


            idx = x.synonym.syns.index(x.hitSyn)
            foundCooc.mirnaFound = None

            if idx >= 0:

                try:
                    test = miRNA(x.synonym.syns[idx])
                    outstr = test.getStringFromParts(
                        [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS,
                         miRNAPART.ARM], normalized=True)
                    foundCooc.mirnaFound = outstr

                except:


                    #sys.stderr.write("cannot parse mirna: " + x.synonym.syns[idx])

                    if __debug__:
                        pass
                        #miRNA(x.synonym.syns[idx])
                        #exit(-1)

                    foundCooc.mirnaFound = None

            if idx < 0 or foundCooc.mirnaFound == None:

                for mirnaSyn in x.synonym.syns:

                    if mirnaSyn.startswith("miR-") and not 'mediated' in mirnaSyn:
                        test = miRNA(mirnaSyn)
                        outstr = test.getStringFromParts([miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM])
                        foundCooc.mirnaFound = outstr
                        break

            miRNALoc = mirnaToSent[x]
            hgncLoc = hgncToSent[y]

            if miRNALoc[0] == hgncLoc[0]:
                foundCooc.sameParagraph = True

                if miRNALoc[1] == hgncLoc[1]:
                    foundCooc.sameSentence = True

                    foundCooc.relation = findRelationBySyns(x,y, sentDB, relHits)

            allCoocs.append(foundCooc)

    return allCoocs

idTuple2Pubmed = defaultdict(set)
orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")


allfiles = glob.glob(resultBase + "/hgnc/*.index")
allfileIDs = [os.path.basename(x).replace(".index", "") for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

#allfileIDs = [894]


def analyseFile(splitFileIDs, env):

    fileCoocs = []


    for splitFileID in splitFileIDs:


        hgncFile = resultBase + "/hgnc/"+splitFileID +".index"
        mirnaFile = resultBase + "/mirna/"+splitFileID +".index"
        relFile = resultBase + "/relations/" + splitFileID + ".index"

        sentFile = "/mnt/c/dev/data/pmc/allsent/"+splitFileID +".sent"

        mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)
        if len(mirnaHits) == 0:
            continue

        hgncHits = SyngrepHitFile(hgncFile, hgncSyns)
        if len(hgncHits) == 0:
            continue

        relHits = SyngrepHitFile(relFile, relSyns)


        sentDB = SentenceDB(sentFile)

        sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

        for docID in mirnaHits:

            if docID in hgncHits:

                mirnaSynHits = mirnaHits.getHitsForDocument(docID)
                hgncSynHits = hgncHits.getHitsForDocument(docID)

                #if docID == 'a27229723':
                #    [print(x.synonyme) for x in hgncSynHits]
                #    [print(x.synonyme) for x in mirnaSynHits]

                foundCoocs = findCooccurrences(str(docID), hgncSynHits, mirnaSynHits, sentDB, relHits)

                fileCoocs += foundCoocs


    sys.stderr.write("Found {cnt} elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs)))

    printed = printStuff(None, fileCoocs, None)

    sys.stderr.write("Found {cnt} (printed: {printed}) elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs), printed=printed))


    return None



threads = 4

print(__debug__, threads)

if __debug__:
    threads = 1
    print("Running on threads:", threads)

def printStuff(old, fileCoocs, env):


    """
    coocinfos = defaultdict(list)
    for cooc in fileCoocs:
        if cooc.relation == None and not cooc.sameSentence:
            continue

        coocinfos[(cooc.gene, cooc.mirnaFound, cooc.mirna)].append((cooc.pubmed, cooc.relation))

    if len(coocinfos) > 0:
        for gene,mirna, mirnaid in coocinfos:
            print("{gene}\t{mirna}\t{mirnaid}\t{info}".format(gene=gene, mirna=mirna, mirnaid=mirnaid, info=coocinfos[(gene, mirna, mirnaid)]))
    """


    setSeenRels = set()

    printed = 0


    for cooc in fileCoocs:


        coocRel = None if cooc.relation == None else tuple(cooc.relation)
        thisCooc = (cooc.gene, cooc.mirnaFound, cooc.mirna, cooc.pubmed, cooc.sameParagraph, cooc.sameSentence, coocRel)

        if thisCooc in setSeenRels:
            continue

        setSeenRels.add(thisCooc)

        print("{gene}\t{mirna}\t{mirnaid}\t{pubmed}\t{sapar}\t{sase}\t{relation}\n".format(
            gene=cooc.gene,
            mirna=cooc.mirnaFound,
            mirnaid=cooc.mirna,
            pubmed=cooc.pubmed,
            sapar=cooc.sameParagraph,
            sase=cooc.sameSentence,
            relation=cooc.relation,
        ), end='', flush=True)

        printed += 1

    return printed




ll = MapReduce(threads)
result = ll.exec( allfileIDs, analyseFile, None, 1, None)

