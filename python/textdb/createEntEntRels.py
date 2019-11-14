import argparse

import sys

import glob
import os

import spacy

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import Counter, defaultdict

import nltk

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../../poreSTAT/")


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

nlp = spacy.load('en')  # create blank Language class #en_core_web_lg


class Cooccurrence:

    def __init__(self):
        self.pubmed = None

        self.ent1 = None
        self.ent2 = None
        self.ent1type = None
        self.ent2type = None
        self.ent1found = None
        self.ent2found = None

        self.sameSentence = False
        self.sameParagraph = False
        self.relation = None
        self.mirnaFound = None

    def __str__(self):
        return "{pub}\t{type}\t{ent1}\t{ent2}\t{ent1type}\t{ent2type}".format(pub=self.pubmed, ent1=self.ent1, ent2=self.ent2, ent1type=self.ent1type, ent2type=self.ent2type)

    def __repr__(self):
        return self.__str__()

    def getIdTuple(self):
        return (self.gene, self.mirna)


def findAssocs(assocs, text, textLoc):
    res = []
    for word in assocs:

        allowedAssocPos = assocs[word]

        if False and not textLoc in allowedAssocPos:
            continue

        if word in text:
            res.append((word, textLoc))

    return res


def findRelationBySyns(ent1Hit, ent2Hit, sentDB, relHits):

    sentHits = relHits[str(ent1Hit.documentID)]

    sentence = sentDB.get_sentence(ent1Hit.documentID)
    #(textBefore, textBetween, textAfter, hitOrder) = sentence.extract_text(mirnaHit.position, hgncHit.position)

    negatedSentence = any([x in sentence.text for x in ['not', 'n\'t', 'nega']])

    allRelations = []

    """
    
    TODO we still have to sort out APOE-/- mice out and LDLR-/-
    
    ... what about relation like in neutrophils?
    
    
    """

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

                relDist = abs(reldep[0].position[0] - ent1Hit.position[0])
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

                relDist = abs(rel.position[0]-ent1Hit.position[0])
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

            if ent1Hit.position[0] < ent2Hit.position[0]:
                assocDir = '12'

                if rel.position[0] < ent1Hit.position[0]:
                    assocDirRel = 'V' + assocDir
                elif ent2Hit.position[0] < rel.position[0]:
                    assocDirRel = assocDir + 'V'
                else:
                    assocDirRel = '1V2'

            else:
                assocDir = '21'

                if rel.position[0] < ent2Hit.position[0]:
                    assocDirRel = 'V' + assocDir
                elif ent1Hit.position[0] < rel.position[0]:
                    assocDirRel = assocDir + 'V'
                else:
                    assocDirRel = '2V1'

            assocSent = str(ent1Hit.documentID)
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
                    ent1Hit.position,
                    ent2Hit.position,
                    rel.position,
                    discoveredBy
                )
            )

    else:
        assocDir = None

        if ent1Hit.position[0] < ent2Hit.position[0]:
            assocDir = '12'
        else:
            assocDir = '21'

        allRelations.append(
            (
                assocDir,
                None,
                None,
                None,
                None,
                negatedSentence,
                ent1Hit.position,
                ent2Hit.position,
                None,
                None
            )
        )

    return allRelations

def getMirnaType(mirnaSyn):
    if re.match('MIPF[0-9]+', mirnaSyn.synonym.id) != None:
        return "MIRNA_FAMILY"
    elif re.match('MIMAT[0-9]+', mirnaSyn.synonym.id) != None:
        return "MIRNA"
    elif re.match('MI[0-9]+', mirnaSyn.synonym.id) != None:
        return  'MIRNA_PRE'
    elif re.match('ORGMIR[0-9]+', mirnaSyn.synonym.id) != None:
        return 'MIRNA_ORGMIR'
    elif re.match('ORGMI[0-9]+', mirnaSyn.synonym.id) != None:
        return 'MIRNA_ORGMIR'
    else:
        return 'UNKNOWN'

def handleHarmonizedNameMirna(x):
    idx = x.synonym.syns.index(x.hitSyn)

    if idx >= 0:

        try:
            test = miRNA(x.synonym.syns[idx])
            outstr = test.getStringFromParts(
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR,
                 miRNAPART.MATURE_SEQS,
                 miRNAPART.ARM], normalized=True)

            return outstr

        except:

            # sys.stderr.write("cannot parse mirna: " + x.synonym.syns[idx])

            if __debug__:
                pass
                # miRNA(x.synonym.syns[idx])
                # exit(-1)

    for mirnaSyn in x.synonym.syns:

        if mirnaSyn.startswith("miR-") and not 'mediated' in mirnaSyn:
            test = miRNA(mirnaSyn)
            outstr = test.getStringFromParts(
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR,
                 miRNAPART.MATURE_SEQS, miRNAPART.ARM], normalized=True)

            return outstr

    if __debug__:
        print("Could not match", x.hitSyn)
    return None

def findCooccurrences(pubmed, ent1Hits, ent2Hits, sentDB, relHits):
    def checkSynHit(synhit):
        if len(synhit.foundSyn) <= 5:
            return synhit.perfectHit == True

        return True

    def chekSynHitMirna(synhit):

        if len(synhit.foundSyn) <= 5:
            foundSyn = synhit.foundSyn.lower()
            return foundSyn.startswith('mir') or foundSyn.startswith('micro')

        return True


    setAllEnt1 = set()
    if args.folderType1.upper() == 'MIRNA':
        setAllEnt1 = set([x for x in ent1Hits if chekSynHitMirna(x)])
    else:
        setAllEnt1 = set([x for x in ent1Hits if checkSynHit(x)])

    setAllEnt2 = set()
    if args.folderType2.upper() == 'MIRNA':
        setAllEnt2 = set([x for x in ent2Hits if chekSynHitMirna(x)])
    else:
        setAllEnt2 = set([x for x in ent2Hits if checkSynHit(x)])


    ent1BySent = defaultdict(list)
    ent2BySent = defaultdict(list)

    ent1ToSent = {}
    ent2ToSent = {}

    for hit in ent1Hits:
        parSenID = (hit.documentID.parID, hit.documentID.senID)
        ent1BySent[parSenID].append(hit)
        ent1ToSent[hit] = parSenID

    for hit in ent2Hits:
        parSenID = (hit.documentID.parID, hit.documentID.senID)
        ent2BySent[parSenID].append(hit)
        ent2ToSent[hit] = parSenID

    allCoocs = []

    pmidRels = relHits.getHitsForDocument(pubmed)

    pmidRelBySent = defaultdict(list)


    if pmidRels != None:
        for rel in pmidRels:
            pmidRelBySent[str(rel.documentID)].append(rel)

    ftype1 = args.folderType1.upper()
    ftype2 = args.folderType2.upper()


    for x in setAllEnt1:
        for y in setAllEnt2:

            foundCooc = Cooccurrence()
            foundCooc.pubmed = pubmed

            foundCooc.ent1type = ftype1
            foundCooc.ent2type = ftype2

            foundCooc.ent1 = x.synonym.id
            foundCooc.ent2 = y.synonym.id

            foundCooc.ent1found = None
            foundCooc.ent2found = None

            if foundCooc.ent1type == 'MIRNA':

                foundEnt = handleHarmonizedNameMirna(x)
                foundCooc.ent1found = foundEnt if foundEnt != None else x.hitSyn

            else:
                foundCooc.ent1found = x.hitSyn

            if foundCooc.ent2type == 'MIRNA':

                foundEnt = handleHarmonizedNameMirna(y)
                foundCooc.ent2found = foundEnt if foundEnt != None else y.hitSyn
            else:
                foundCooc.ent2found = y.hitSyn


            if foundCooc.ent1type == 'MIRNA':
                foundCooc.ent1 = foundCooc.ent1found
                foundCooc.ent1found = x.hitSyn

            if foundCooc.ent2type == 'MIRNA':
                foundCooc.ent2 = foundCooc.ent2found
                foundCooc.ent2found = y.hitSyn


            ent1Loc = ent1ToSent[x]
            ent2Loc = ent2ToSent[y]

            if ent1Loc[0] == ent2Loc[0]:
                foundCooc.sameParagraph = True

                if ent1Loc[1] == ent2Loc[1]:
                    foundCooc.sameSentence = True

                    foundCooc.relation = findRelationBySyns(x, y, sentDB, pmidRelBySent)

            allCoocs.append(foundCooc)

    return allCoocs

def analyseFile(splitFileIDs, env):

    fileCoocs = []

    for splitFileID in splitFileIDs:

        ent1File = resultBase + "/"+args.folder1+"/" + splitFileID + ".index"
        ent2File = resultBase + "/"+args.folder2+"/" + splitFileID + ".index"
        relFile = resultBase + "/relations/" + splitFileID + ".index"

        sentFile = args.sentdir + "/" + splitFileID + ".sent"

        ent1Hits = SyngrepHitFile(ent1File, ent1Syns)
        if len(ent1Hits) == 0:
            continue

        ent2Hits = SyngrepHitFile(ent2File, ent2Syns)
        if len(ent2Hits) == 0:
            continue

        relHits = SyngrepHitFile(relFile, relSyns)

        # only load sentences if there's a hit ...
        sentDB = None

        sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

        for docID in ent1Hits:

            if docID in ent2Hits:

                if sentDB == None:
                    sentDB = SentenceDB(sentFile)

                ent1SynHits = ent1Hits.getHitsForDocument(docID)
                ent2SynHits = ent2Hits.getHitsForDocument(docID)

                # if docID == 'a27229723':
                #    [print(x.synonyme) for x in hgncSynHits]
                #    [print(x.synonyme) for x in mirnaSynHits]

                foundCoocs = findCooccurrences(str(docID), ent1SynHits, ent2SynHits, sentDB, relHits)

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

    parser.add_argument('-f1', '--folder1', type=str, help='entity 1: hgnc, mirna', default="hgnc", required=False)
    parser.add_argument('-f2', '--folder2', type=str, help='entity 2: mgi, mirna', default="mirna", required=False)

    parser.add_argument('-ft1', '--folderType1', type=str, help='entity type 1: entity: mirna, gene, lncrna, ...', default="gene", required=False)
    parser.add_argument('-ft2', '--folderType2', type=str, help='entity type 2: entity: mirna', default="mirna", required=False)


    args = parser.parse_args()

    #resultBase = dataDir + "/miRExplore/textmine/results_pmc/"
    resultBase = args.resultdir
    dataDir = args.datadir

    ent1Syns = SynfileMap(resultBase + "/"+args.folder1+"/synfile.map")
    ent1Syns.loadSynFiles(('/mnt/c/ownCloud/data', dataDir))

    ent2Syns = SynfileMap(resultBase + "/"+args.folder2+"/synfile.map")
    ent2Syns.loadSynFiles(('/mnt/c/ownCloud/data', dataDir))

    relSyns = SynfileMap(resultBase + "/relations/synfile.map")
    relSyns.loadSynFiles(('/mnt/c/ownCloud/data', dataDir))

    relationSyns = AssocSynfile(args.datadir + '/miRExplore/relations/allrels.csv')


    idTuple2Pubmed = defaultdict(set)
    orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

    allfiles = glob.glob(resultBase + "/"+args.folder1+"/*.index")
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
                cooc.ent1, cooc.ent1type, cooc.ent1found, cooc.ent2, cooc.ent2type, cooc.ent2found, cooc.pubmed, cooc.sameParagraph, cooc.sameSentence, coocRel
            )

            if thisCooc in setSeenRels:
                continue

            setSeenRels.add(thisCooc)

            print("{ent1}\t{ent1found}\t{ent1type}\t{ent2}\t{ent2found}\t{ent2type}\t{pubmed}\t{sapar}\t{sase}\t{relation}\n".format(
                ent1=cooc.ent1,
                ent2=cooc.ent2,
                ent1found=cooc.ent1found,
                ent2found=cooc.ent2found,
                ent1type=cooc.ent1type,
                ent2type=cooc.ent2type,
                pubmed=cooc.pubmed,
                sapar=cooc.sameParagraph,
                sase=cooc.sameSentence,
                relation=cooc.relation,
            ), end='', flush=True)

            printed += 1

        return printed


    ll = MapReduce(threads)
    result = ll.exec(allfileIDs, analyseFile, None, 1, None)