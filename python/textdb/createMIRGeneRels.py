import argparse

import sys

import glob
import os

import spacy

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


from utils.parallel import MapReduce
from enum import Enum

nlp = spacy.load('en')  # create blank Language class


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
        return "{pub}\t{type}\t{geneid}\t{mirnaid}\t{mirnaname}".format(pub=self.pubmed, type=self.idtype,
                                                                        geneid=self.gene, mirnaid=self.mirna,
                                                                        mirnaname=self.mirnadesc)

    def __repr__(self):
        return self.__str__()

    def getIdTuple(self):
        return (self.gene, self.mirna, self.idtype)


def findAssocs(assocs, text, textLoc):
    res = []
    for word in assocs:

        allowedAssocPos = assocs[word]

        if False and not textLoc in allowedAssocPos:
            continue

        if word in text:
            res.append((word, textLoc))

    return res


def findRelationBySyns(mirnaHit, hgncHit, sentDB, relHits):

    sentHits = relHits[str(mirnaHit.documentID)]

    sentence = sentDB.get_sentence(mirnaHit.documentID)
    #(textBefore, textBetween, textAfter, hitOrder) = sentence.extract_text(mirnaHit.position, hgncHit.position)

    negatedSentence = any([x in sentence.text for x in ['not', 'n\'t', 'nega']])

    allRelations = []

    if len(sentHits) > 0:

        doc = nlp(sentence.text)
        alldeps = [(t.idx, t.text, t.dep_, t.pos_, t.head.text) for t in doc]
        verbDeps = [x for x in alldeps if x[3] == 'VERB']

        assocRelDeps = []

        for rel in sentHits:

            for dep in verbDeps:

                if rel.position[0] <= dep[0] and dep[0]+len(dep[1]) <= rel.position[1]:

                    if dep[3] == 'VERB':
                        assocDep = dep
                        assocRel = rel

                        assocRelDeps.append((assocRel, assocDep))

                        break


        depHits = []

        if len(assocRelDeps) != 0:

            curDist = len(sentence.text) + 1
            curRel = None
            for reldep in assocRelDeps:

                relDist = abs(reldep[0].position[0] - mirnaHit.position[0])
                if relDist < curDist:
                    curDist = relDist
                    curRel = reldep

            if curRel != None:
                depHits = [curRel[0]]
                discoveredBy = "spacy"


        elif len(assocRelDeps) == 0 or len(depHits) == 0:

            curDist = len(sentence.text)+1
            curRel = None
            for rel in sentHits:

                relDist = abs(rel.position[0]-mirnaHit.position[0])
                if relDist < curDist:
                    curDist = relDist
                    curRel = rel
                    depHits = [curRel]

                    discoveredBy = "all_rels"



        for rel in depHits:
            assocDir = None
            assocType = None
            assocWord = None
            assocSent = None
            assocDirRel = None

            if mirnaHit.position[0] < hgncHit.position[0]:
                assocDir = 'MG'

                if rel.position[0] < mirnaHit.position[0]:
                    assocDirRel = 'V' + assocDir
                elif hgncHit.position[0] < rel.position[0]:
                    assocDirRel = assocDir + 'V'
                else:
                    assocDirRel = 'MVG'

            else:
                assocDir = 'GM'

                if rel.position[0] < hgncHit.position[0]:
                    assocDirRel = 'V' + assocDir
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
                    rel.position,
                    discoveredBy
                )
            )

    else:
        assocDir = None

        if mirnaHit.position[0] < hgncHit.position[0]:
            assocDir = 'MG'
        else:
            assocDir = 'GM'

        allRelations.append(
            (
                assocDir,
                None,
                None,
                None,
                None,
                negatedSentence,
                mirnaHit.position,
                hgncHit.position,
                None,
                None
            )
        )

    return allRelations


def findCooccurrences(pubmed, hgncHits, mirnaHits, sentDB, relHits):
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

        hgncToSent[hit] = parSenID

    for hit in mirnaHits:
        parSenID = (hit.documentID.parID, hit.documentID.senID)
        mirnaBySent[parSenID].append(hit)

        mirnaToSent[hit] = parSenID

    allCoocs = []

    pmidRels = relHits.getHitsForDocument(pubmed)

    pmidRelBySent = defaultdict(list)


    if pmidRels != None:
        for rel in pmidRels:
            pmidRelBySent[str(rel.documentID)].append(rel)


    for x in setAllMirnas:
        for y in setAllGenes:

            foundCooc = Cooccurrence()
            foundCooc.pubmed = pubmed

            if re.match('MIPF[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = "MIRNA_FAMILY"
            elif re.match('MIMAT[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = "MIRNA"
            elif re.match('MI[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = 'MIRNA_PRE'
            elif re.match('ORGMIR[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = 'MIRNA_ORGMIR'
            elif re.match('ORGMI[0-9]+', x.synonym.id) != None:
                foundCooc.idtype = 'MIRNA_ORGMIR'
            else:
                foundCooc.idtype = 'UNKNOWN'

            foundCooc.gene = y.synonym.id
            foundCooc.mirna = x.synonym.id
            foundCooc.mirnadesc = str(x.synonym)

            foundCooc.mirnaFound = x.hitSyn

            idx = x.synonym.syns.index(x.hitSyn)
            foundCooc.mirnaFound = None

            if idx >= 0:

                try:
                    test = miRNA(x.synonym.syns[idx])
                    outstr = test.getStringFromParts(
                        [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR,
                         miRNAPART.MATURE_SEQS,
                         miRNAPART.ARM], normalized=True)
                    foundCooc.mirnaFound = outstr

                except:

                    # sys.stderr.write("cannot parse mirna: " + x.synonym.syns[idx])

                    if __debug__:
                        pass
                        # miRNA(x.synonym.syns[idx])
                        # exit(-1)

                    foundCooc.mirnaFound = None

            if idx < 0 or foundCooc.mirnaFound == None:

                for mirnaSyn in x.synonym.syns:

                    if mirnaSyn.startswith("miR-") and not 'mediated' in mirnaSyn:
                        test = miRNA(mirnaSyn)
                        outstr = test.getStringFromParts(
                            [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR,
                             miRNAPART.MATURE_SEQS, miRNAPART.ARM])
                        foundCooc.mirnaFound = outstr
                        break

            miRNALoc = mirnaToSent[x]
            hgncLoc = hgncToSent[y]

            if miRNALoc[0] == hgncLoc[0]:
                foundCooc.sameParagraph = True

                if miRNALoc[1] == hgncLoc[1]:
                    foundCooc.sameSentence = True

                    foundCooc.relation = findRelationBySyns(x, y, sentDB, pmidRelBySent)

            allCoocs.append(foundCooc)

    return allCoocs

def analyseFile(splitFileIDs, env):

    fileCoocs = []

    for splitFileID in splitFileIDs:

        hgncFile = resultBase + "/"+args.folderG+"/" + splitFileID + ".index"
        mirnaFile = resultBase + "/"+args.folderM+"/" + splitFileID + ".index"
        relFile = resultBase + "/relations/" + splitFileID + ".index"

        sentFile = args.sentdir + "/" + splitFileID + ".sent"

        mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)
        if len(mirnaHits) == 0:
            continue

        hgncHits = SyngrepHitFile(hgncFile, hgncSyns)
        if len(hgncHits) == 0:
            continue

        relHits = SyngrepHitFile(relFile, relSyns)

        # only load sentences if there's a hit ...
        sentDB = None

        sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

        for docID in mirnaHits:

            if docID in hgncHits:

                if sentDB == None:
                    sentDB = SentenceDB(sentFile)

                mirnaSynHits = mirnaHits.getHitsForDocument(docID)
                hgncSynHits = hgncHits.getHitsForDocument(docID)

                # if docID == 'a27229723':
                #    [print(x.synonyme) for x in hgncSynHits]
                #    [print(x.synonyme) for x in mirnaSynHits]

                foundCoocs = findCooccurrences(str(docID), hgncSynHits, mirnaSynHits, sentDB, relHits)

                fileCoocs += foundCoocs

    sys.stderr.write("Found {cnt} elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs)))

    printed = printStuff(None, fileCoocs, None)

    thisProcID = str(os.getpid())
    sys.stderr.write("{procID}: Found {cnt} (printed: {printed}) elems in files {ids}\n".format(
        cnt=str(len(fileCoocs)),
        ids=str(splitFileIDs),
        printed=printed,
        procID=thisProcID))

    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='aggregate tm results', add_help=False)
    parser.add_argument('-s', '--sentdir', type=str, help='where are the sentences?', required=True)
    parser.add_argument('-r', '--resultdir', type=str, help='where are all the index-files?', required=True)
    parser.add_argument('-d', '--datadir', type=str, help='where is te miRExplore bsae?', required=True)

    parser.add_argument('-f1', '--folderM', type=str, help='where is te miRExplore bsae?', default="mirna", required=False)
    parser.add_argument('-f2', '--folderG', type=str, help='where is te miRExplore bsae?', default="hgnc", required=False)


    args = parser.parse_args()

    #resultBase = dataDir + "/miRExplore/textmine/results_pmc/"
    resultBase = args.resultdir
    dataDir = args.datadir

    mirnaSyns = SynfileMap(resultBase + "/"+args.folderM+"/synfile.map")
    mirnaSyns.loadSynFiles(('/home/users/joppich/ownCloud/data/', dataDir))

    hgncSyns = SynfileMap(resultBase + "/"+args.folderG+"/synfile.map")
    hgncSyns.loadSynFiles(('/home/users/joppich/ownCloud/data/', dataDir))

    relSyns = SynfileMap(resultBase + "/relations/synfile.map")
    relSyns.loadSynFiles(('/home/users/joppich/ownCloud/data/', dataDir))

    relationSyns = AssocSynfile(args.datadir + '/miRExplore/relations/allrels.csv')


    idTuple2Pubmed = defaultdict(set)
    orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

    allfiles = glob.glob(resultBase + "/"+args.folderG+"/*.index")
    allfileIDs = [os.path.basename(x).replace(".index", "") for x in allfiles]
    allfileIDs = sorted(allfileIDs, reverse=True)


    # allfileIDs = [894]
    threads = 4

    if __debug__:
        threads = 1
        sys.stderr.write("Running on threads:" + str(threads) + "\n")

    sys.stderr.write("Debug Mode? " + str(__debug__) + " and threads " + str(threads) + "\n")


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
            thisCooc = (
                cooc.gene, cooc.mirnaFound, cooc.mirna, cooc.pubmed, cooc.sameParagraph, cooc.sameSentence, coocRel
            )

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
    result = ll.exec(allfileIDs, analyseFile, None, 1, None)