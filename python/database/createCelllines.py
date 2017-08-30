from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

from database.Neo4JInterface import neo4jInterface
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID

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
    id2species[oboID].union(taxID)

    if oboRels != None:
        for rel in oboRels:
            term = rel.term

            id2derived_from[oboID].add(term.id)

db = neo4jInterface(simulate=True)

for id in id2node:
    node = id2node[id]
    db.createNodeIfNotExists(['CELLLINE'], node)

    allSpecies = id2species.get(id, set())
    cellLineUnique = len(allSpecies) == 1

    for species in allSpecies:
        db.createRelationship('tax', ['TAXID'], {'id': species}, 'cell', ['CELLLINE'], node, ['IS_ORGANISM'], {'unique': cellLineUnique})

for id in id2derived_from:

    allDerivatives = id2derived_from[id]

    for deriv in allDerivatives:

        if not deriv in id2node:
            continue

        db.createRelationship('id', ['CELLLINE'], {'id': id}, 'other', ['CELLLINE'], {'id': deriv}, ['CELLINE_DERIVED_FROM'], None)

db.close()