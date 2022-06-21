from nertoolkit.geneontology.GeneOntology import GeneOntology

from utils.parallel import MapReduce
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


def getSynIDs( oboLocation, selOboIDs):

    ontology = GeneOntology(oboLocation)

    allSynIDs = []


    for termID in ontology.dTerms:

        oboNode = ontology.dTerms[termID]

        oboID = oboNode.id
        oboName = oboNode.name

        oboSyns = oboNode.synonym
        oboRels = oboNode.is_a

        if oboID in selOboIDs:

            allSynIDs.append(oboID)

            allchildren = oboNode.getAllChildren()

            allSynIDs += [x.term.id for x in allchildren]

    allSynIDs = [x.replace(":", '_') for x in allSynIDs]

    return allSynIDs



neutrophilSynIDs = getSynIDs(dataDir + "miRExplore/foundational_model_anatomy/fma_obo.obo",  ['FMA:62860'])
tissueIDs = getSynIDs(dataDir + "miRExplore/foundational_model_anatomy/fma_obo.obo",  ['FMA:67498', 'FMA:9637', 'FMA:68646'])
doidSynIDs = getSynIDs(dataDir + "miRExplore/doid.obo",  ['DOID:104'])
goSynIDs = getSynIDs(dataDir + "miRExplore/textmine/neutrophils.obo",  ['NP:001'])

for x in neutrophilSynIDs:
    if x in tissueIDs:
        tissueIDs.remove(x)

print(neutrophilSynIDs)
print(tissueIDs)
print(doidSynIDs)
print(goSynIDs)


resultBase = dataDir + "/miRExplore/textmine/results/"
fmaSyns = SynfileMap(resultBase + "/model_anatomy/synfile.map")
fmaSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

doidSyns = SynfileMap(resultBase + "/disease/synfile.map")
doidSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

goSyns = SynfileMap(resultBase + "/neutrophils/synfile.map")
goSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )



allfiles = glob.glob(resultBase + "/hgnc/pubmed18n*.index")
allfileIDs = [int(os.path.basename(x).replace('pubmed18n', '').replace('.index','')) for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

#allfileIDs = allfileIDs[0:10]

def testMentioned(allhits, targetSynIDs):

    if allhits == None:
        return False

    for hit in allhits:

        if hit.synonym.id in targetSynIDs:
            return True

    return False

#allfileIDs = [700]

def analyseFile(splitFileIDs, env):


    subject2pmids = defaultdict(set)


    for splitFileID in splitFileIDs:

        fileID = "{:>4}".format(splitFileID).replace(" ", "0")

        fmaFile = resultBase + "/model_anatomy/pubmed18n"+fileID+".index"
        doidFile = resultBase + "/disease/pubmed18n"+fileID+".index"
        goFile = resultBase + "/neutrophils/pubmed18n"+fileID+".index"

        sentFile = "/mnt/c/dev/data/pubmed/pubmed18n" + fileID + ".sent"

        fmaHits = SyngrepHitFile(fmaFile, fmaSyns)
        if len(fmaHits) == 0:
            return

        print("Processing file: ", fileID)

        doidHits = SyngrepHitFile(doidFile, doidSyns)
        goHits = SyngrepHitFile(goFile, goSyns)

        sentDB = SentenceDB(sentFile)

        for docID in fmaHits:

            docFMAHits = fmaHits.getHitsForDocument(docID)
            docDoidHits = doidHits.getHitsForDocument(docID)
            docGOHits = goHits.getHitsForDocument(docID)

            # no tissue hits at all ...
            if len(docFMAHits) == 0:
                continue

            # is neutrophil hit? if no -> continue
            neutrophilMentioned = testMentioned(docFMAHits, neutrophilSynIDs)

            if not neutrophilMentioned:
                continue

            subject2pmids['NEUTROPHIL'].add(int(docID))

            if testMentioned(docFMAHits, tissueIDs):
                subject2pmids['TISSUES'].add(int(docID))

            if testMentioned(docDoidHits, doidSynIDs):
                subject2pmids['DOID'].add(int(docID))

            if testMentioned(docGOHits, goSynIDs):
                subject2pmids['INFLAMM'].add(int(docID))

    return subject2pmids





finalResults = defaultdict(set)

def reduceSets(old, new, env):

    if old == None:
        return new

    for x in new:
        if x in old:
            old[x] = old[x] | new[x]

        else:
            old[x] = new[x]

    return old


ll = MapReduce(4)
result = ll.exec( allfileIDs, analyseFile, None, 1, reduceSets)


with open("/tmp/tm_soehnlein", 'w') as fout:
    for x in result:
        fout.write(str(x) + "\t" + str(result[x]) + "\n")