from mjoppich.geneontology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir, speciesName2TaxID

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


allSet = {'all'}
for splitFileID in range(892, 0, -1):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    diseaseHitsFile = resultBase + "/cellline/medline17n"+fileID+".index"

    hitsFile = SyngrepHitFile(diseaseHitsFile, celllinesMap)

    if len(hitsFile) == 0:
        continue

    print("Document: " + str(fileID))

    for docID in hitsFile:
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
            db.createNodeIfNotExists(['EVIDENCE', 'PUBMED'], {'id': docID})
            db.createRelationship('cellline', ['CELLLINE'], {'id': celllineID}, 'pubmed', ['PUBMED'], {'id': docID}, ['CELLLINE_MENTION'], None)

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


db.close()