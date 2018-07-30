from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

geneOntologyObo = GeneOntology(dataDir + "miRExplore/go/go.obo")
tax2cells = defaultdict(set)

allNameSpaces = set()
for oboNode in geneOntologyObo:

    for ns in oboNode.namespace:
        allNameSpaces.add( ns )

ns2dbns = {}
for ns in allNameSpaces:
    ns2dbns[ns] = ns.upper().replace(' ', '_')


db = neo4jInterface(simulate=False)
db.deleteRelationship('n', ['GO'], None, 'm', ['GO'], None, ['GO_ISA', None])
db.deleteNode(['GO'], None)
db.createUniquenessConstraint('GO', 'id')


createdNodes = set()

for oboNode in geneOntologyObo:

    oboID = oboNode.id
    oboName = oboNode.name
    oboRels = oboNode.is_a

    oboNS = oboNode.namespace

    dbNS = [ns2dbns[x] for x in oboNS]

    db.createNodeIfNotExists(['GO'] + dbNS, {'id': oboID, 'name': oboName})
    createdNodes.add( oboID )

for oboNode in geneOntologyObo:

    oboID = oboNode.id
    oboRels = oboNode.is_a

    for child in oboRels:
        db.createRelationship('gosrc', ['GO'], {'id': oboID}, 'gotgt', ['GO'], {'id': child.id},
                          ['GO_ISA'], None)

db.close()