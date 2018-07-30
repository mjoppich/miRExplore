import glob
import os

from mjoppich.geneontology import GeneOntology
from porestat.utils.Parallel import MapReduce

from database.Neo4JInterface import neo4jInterface
from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir, speciesName2TaxID, eprint

resultBase = dataDir + "/miRExplore/textmine/results/"
celllinesMap = SynfileMap(resultBase + "/cellline/synfile.map")
celllinesMap.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )


knownTaxIDs = set()
knownTaxIDs.add('all')
for org in speciesName2TaxID:
    knownTaxIDs.add(str(speciesName2TaxID[org]))

synfileID2tax = {}
for synfileID in celllinesMap.synfiles:
    synfileName = celllinesMap.synfiles[synfileID]

    hitOrgs = []
    for org in knownTaxIDs:
        if "."+org+"." in synfileName:
            hitOrgs.append(org)

    if len(hitOrgs) != 1:
        print("No or multiple files for org: " + str(synfileName) + " " + str(hitOrgs))
    else:
        org = hitOrgs[0]
        synfileID2tax[synfileID] = org

print(synfileID2tax)

celllinesObo = GeneOntology(dataDir + "miRExplore/cellosaurus/cellosaurus.obo")

db = neo4jInterface(simulate=False, printQueries=False)
db.deleteRelationship('n', ['CELLLINE'], None, 'm', ['PUBMED'], None, ['CELLLINE_MENTION'], None)

addUnknownPubmeds = False
allSet = {'all'}

allfiles = glob.glob(resultBase + "/hgnc/medline17n*.index")
allfileIDs = [int(os.path.basename(x).replace('medline17n', '').replace('.index','')) for x in allfiles]
allfileIDs = sorted(allfileIDs, reverse=True)

retVal = db.matchNodes(['PUBMED'], None, nodename='n')
relevantPMIDs = set()

for x in retVal:

    nodeData = x['n']
    if 'id' in nodeData.properties:
        pmid = nodeData.properties['id']
        relevantPMIDs.add(pmid)
    else:
        eprint("No data in: ", str(nodeData))

db.close()

if len(relevantPMIDs) == 0:
    eprint("No RELEVANT PUBMED entries found")


def analyseFile(splitFileID, relPMIDs):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    diseaseHitsFile = resultBase + "/cellline/medline17n"+fileID+".index"

    hitsFile = SyngrepHitFile(diseaseHitsFile, celllinesMap)

    if len(hitsFile) == 0:
        return

    print("Start Document: " + str(fileID))

    procDB = neo4jInterface(simulate=False, printQueries=False)

    for docID in hitsFile:

        if not docID in relPMIDs:
            continue

        synHits = hitsFile.getHitsForDocument(docID)

        foundUniqueHits = set()
        foundOrgs = set()

        for hit in synHits:

            if len(hit.foundSyn) < 5:
                if not hit.perfectHit:
                    continue
            hitSynFileID = hit.synonymID.synfile
            foundOrgs.add( synfileID2tax[hitSynFileID] )

            hitSyn = hit.synonym
            foundUniqueHits.add(hitSyn.id)

        if len(foundUniqueHits) == 0:
            continue

        for celllineID in foundUniqueHits:

            pubmedExists=False
            if addUnknownPubmeds:
                procDB.createNodeIfNotExists(['EVIDENCE', 'PUBMED'], {'id': docID})
                pubmedExists = True
            else:
                if procDB.nodeExists(['PUBMED'], {'id': docID}):
                    pubmedExists = True


            if pubmedExists:
                res = procDB.createRelationship('cellline', ['CELLLINE'], {'id': celllineID}, 'pubmed', ['PUBMED'], {'id': docID}, ['CELLLINE_MENTION'], None)
                print("Add: ", fileID, docID, celllineID, [x for x in res if res != None])

        foundOrgs = foundOrgs.difference(allSet)

        if len(foundOrgs) == 1:
            pass
            # create relation
            # print('Associate: ' + str(foundOrgs))
        elif len(foundOrgs) == 0:
            pass
        elif len(foundOrgs) > 1:
            # print('Ambiguous pubmed: ' + docID)
            pass

    print("End Document: " + str(fileID))
    procDB.close()



ll = MapReduce(6)
result = ll.exec( allfileIDs, analyseFile, relevantPMIDs, 1, None)