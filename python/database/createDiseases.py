from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

diseaseObo = GeneOntology(dataDir + "miRExplore/doid.obo")
tax2cells = defaultdict(set)

id2node = {}
id2derived_from = defaultdict(set)


for cellID in diseaseObo.dTerms:

    oboNode = diseaseObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name
    oboRels = oboNode.is_a

    id2node[oboID] = {'id': oboID, 'name': oboName}

    if oboRels != None:
        for rel in oboRels:
            term = rel.term
            id2derived_from[oboID].add(term.id)

db = neo4jInterface(simulate=False)
db.deleteRelationship('n', ['DISEASE'], None, 'm', ['DISEASE'], None, ['DISEASE_DERIVED_FROM', None])
db.deleteNode(['DISEASE'], None)
db.createUniquenessConstraint('DISEASE', 'id')


for id in id2node:
    node = id2node[id]
    db.createNodeIfNotExists(['DISEASE'], node)

for id in id2derived_from:

    allDerivatives = id2derived_from[id]

    for deriv in allDerivatives:

        if not deriv in id2node:
            continue

        db.createRelationship('disease', ['DISEASE'], {'id': id}, 'other', ['DISEASE'], {'id': deriv}, ['DISEASE_DERIVED_FROM'], None)

db.close()