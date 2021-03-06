from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID, eprint

celloObo = GeneOntology(dataDir + "miRExplore/cellosaurus/cellosaurus.obo")
tax2cells = defaultdict(set)

id2node = {}
id2species = defaultdict(set)
id2derived_from = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]

    oboID = oboNode.id
    oboName = oboNode.name

    oboSyns = oboNode.synonym
    oboXRefs = oboNode.xref
    oboRels = oboNode.is_a

    taxID = set()
    if oboXRefs != None:
        for xref in oboXRefs:
            if xref.startswith('NCBI_TaxID'):

                newTaxID = xref.split(' ')[0].split(':')[1]

                if newTaxID in allowedTaxIDs:
                    taxID.add(newTaxID)

    id2node[oboID] = {'id': oboID, 'name': oboName}

    for tax in taxID:
        id2species[oboID].add(tax)

    if oboRels != None:
        for rel in oboRels:
            term = rel.term

            id2derived_from[oboID].add(term.id)

db = neo4jInterface(simulate=False, printQueries=False)
db.deleteRelationship('n', ['CELLLINE'], None, 'm', ['CELLLINE'], None, ['CELLINE_DERIVED_FROM'], None)
db.deleteRelationship('n', ['TAX'], None, 'm', ['CELLLINE'], None, ['HAS_CELLLINE'], None)
db.deleteNode(['CELLLINE'], None)
db.createUniquenessConstraint('CELLLINE', 'id')


for id in id2node:
    node = id2node[id]
    db.createNodeIfNotExists(['CELLLINE'], node)

    allSpecies = id2species.get(id, set())
    cellLineUnique = len(allSpecies) == 1

    for species in allSpecies:
        try:
            taxID = int(species)
            db.createRelationship('tax', ['TAX'], {'id': taxID}, 'cell', ['CELLLINE'], node, ['HAS_CELLLINE'],
                                  {'unique': cellLineUnique})
        except:
            eprint(str(species) + "is not a valid tax id in database")
            continue


for id in id2derived_from:

    allDerivatives = id2derived_from[id]

    for deriv in allDerivatives:

        if not deriv in id2node:
            eprint("Not in id2node: " + str(deriv))
            continue

        db.createRelationship('id', ['CELLLINE'], {'id': id}, 'other', ['CELLLINE'], {'id': deriv}, ['CELLINE_DERIVED_FROM'], None)

db.close()