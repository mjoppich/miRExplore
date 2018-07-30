import glob
import sys
import os

from nertoolkit.geneontology.GeneOntology import GeneOntology, GOTerm

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
diseaseSyns = SynfileMap(resultBase + "/cellline/synfile.map")
diseaseSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )


allfiles = glob.glob(resultBase + "/cellline/*.index")
allfileIDs = [os.path.basename(x).replace(".index", "") for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

#allfileIDs = [894]

celloObo = GeneOntology(dataDir + "miRExplore/meta_cells.obo")


def getTerm(synid, obo):

    if synid in obo.dTerms:
        return obo.getID(synid)

    synid = synid.replace('_', ':', 1)

    return obo.getID(synid)

def analyseFile(splitFileIDs, env):

    fileCoocs = []


    for splitFileID in splitFileIDs:


        diseaseFile = resultBase + "/cellline/"+splitFileID +".index"

        diseaseHits = SyngrepHitFile(diseaseFile, diseaseSyns)
        if len(diseaseHits) == 0:
            continue

        sentFile = "/mnt/c/dev/data/pmc/allsent/"+splitFileID +".sent"
        sentDB = SentenceDB(sentFile)

        sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

        for docID in diseaseHits:

            docHits = diseaseHits.getHitsForDocument(docID)

            allSynIDs = set()
            for hit in docHits:

                if "and " in hit.foundSyn:
                    continue

                allSynIDs.add(hit.synonym.id)


            removeIDs = set()
            for synID in allSynIDs:

                gterm = getTerm(synID, celloObo)
                if gterm == None:
                    sys.stderr.write("Invalid synID: " + synID)
                    removeIDs.add(synID)
                    continue

                allChildren = gterm.getAllChildren(maxLevel=2)

                for x in allChildren:
                    removeIDs.add(x)


            allowedIDs = [x for x in allSynIDs if not x in removeIDs]
            #allowedIDs = allSynIDs.remove(removeIDs)

            allterms = []
            for synID in allowedIDs:

                gterm = getTerm(synID, celloObo)

                allterms.append(gterm)

                fileCoocs.append((docID, gterm.id, gterm.name))



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

    printed = 0

    for cooc in fileCoocs:

        print("{pmid}\t{cl}\t{name}\n".format(
            pmid=cooc[0],
            cl=cooc[1],
            name=cooc[2]
        ), end='', flush=True)

        printed += 1

    return printed




ll = MapReduce(threads)
result = ll.exec( allfileIDs, analyseFile, None, 1, None)

