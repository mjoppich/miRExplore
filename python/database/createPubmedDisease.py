import glob
import os

from mjoppich.geneontology import GeneOntology
from porestat.utils.Parallel import MapReduce

from database.Neo4JInterface import neo4jInterface
from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir, eprint

resultBase = dataDir + "/miRExplore/textmine/results/"
diseaseMap = SynfileMap(resultBase + "/disease/synfile.map")
diseaseMap.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )
diseaseObo = GeneOntology(dataDir + "miRExplore/doid.obo")

db = neo4jInterface(simulate=False)
db.deleteRelationship('n', ['DISEASE'], None, 'm', ['PUBMED'], None, ['DISEASE_MENTION'], None)

allfiles = glob.glob(resultBase + "/hgnc/medline17n*.index")
allfileIDs = [int(os.path.basename(x).replace('medline17n', '').replace('.index','')) for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

addUnknownPubmeds=False

retVal = db.matchNodes(['PUBMED'], None, nodename='n')
relevantPMIDs = set()

for x in retVal:

    nodeData = x['n']
    if 'id' in nodeData.properties:
        pmid = nodeData.properties['id']
        relevantPMIDs.add(pmid)
    else:
        eprint("No data in: ", str(nodeData))

if len(relevantPMIDs) == 0:
    eprint("No RELEVANT PUBMED entries found")

def analyseFile(splitFileID, relPMIDs):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    diseaseHitsFile = resultBase + "/disease/medline17n"+fileID+".index"

    hitsFile = SyngrepHitFile(diseaseHitsFile, diseaseMap)

    if len(hitsFile) == 0:
        return

    print("Document: " + str(fileID))
    print("Start Document: " + str(fileID))

    procDB = neo4jInterface(simulate=False, printQueries=False)

    for docID in hitsFile:

        if not docID in relPMIDs:
            continue

        synHits = hitsFile.getHitsForDocument(docID)

        foundUniqueHits = set()
        for hit in synHits:

            if len(hit.foundSyn) < 5:
                if not hit.perfectHit:
                    continue

            hitSyn = hit.synonym
            foundUniqueHits.add(hitSyn.id.replace('_', ':'))

        for synonymID in foundUniqueHits:

            pubmedExists=False
            if addUnknownPubmeds:
                procDB.createNodeIfNotExists(['EVIDENCE', 'PUBMED'], {'id': docID})
                pubmedExists = True
            else:
                if procDB.nodeExists(['PUBMED'], {'id': docID}):
                    pubmedExists = True


            if pubmedExists:
                res = procDB.createRelationship('disease', ['DISEASE'], {'id': synonymID}, 'pubmed', ['PUBMED'], {'id': docID}, ['DISEASE_MENTION'], None)
                print("Add: ", fileID, docID, synonymID, [x for x in res if res != None])


    print("End Document: " + str(fileID))
    procDB.close()


db.close()

ll = MapReduce(6)
result = ll.exec( allfileIDs, analyseFile, relevantPMIDs, 1, None)

