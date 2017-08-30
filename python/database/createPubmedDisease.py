from mjoppich.geneontology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir

resultBase = dataDir + "/miRExplore/textmine/results/"
diseaseMap = SynfileMap(resultBase + "/disease/synfile.map")
diseaseMap.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )
diseaseObo = GeneOntology(dataDir + "miRExplore/doid.obo")

db = neo4jInterface(simulate=True)

for splitFileID in range(892, 0, -1):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    diseaseHitsFile = resultBase + "/disease/medline17n"+fileID+".index"

    hitsFile = SyngrepHitFile(diseaseHitsFile, diseaseMap)

    if len(hitsFile) == 0:
        continue

    for docID in hitsFile:
        synHits = hitsFile.getHitsForDocument(docID)

        foundUniqueHits = set()
        for hit in synHits:

            if len(hit.foundSyn) < 5:
                if not hit.perfectHit:
                    continue

            hitSyn = hit.synonym
            foundUniqueHits.add(hitSyn.id)

        for synonymID in foundUniqueHits:
            db.createNodeIfNotExists(['EVIDENCE', 'PUBMED'], {'id': docID})
            db.createRelationship('disease', ['DISEASE'], {'id': synonymID}, 'pubmed', ['PUBMED'], {'id': docID}, ['DISEASE_MENTION'], None)


db.close()