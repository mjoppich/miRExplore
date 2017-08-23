from collections import Counter
from porestat.utils.DataFrame import DataFrame
from utils.idutils import ltype2label, makeDBGeneID, mirtarbase_exp_type, mirtarbase_function_label, speciesName2TaxID, \
    dataDir
from database.Neo4JInterface import neo4jInterface

mirtarbaseEvidences = DataFrame.parseFromFile(dataDir + "/miRExplore/miRTarBase.csv", bConvertTextToNumber=False)

print(mirtarbaseEvidences.getHeader())

experimentTypes = Counter()
supportTypes = Counter()
referencesWithComma = Counter()

db = neo4jInterface(simulate=True, printQueries=False)

for mirnaEvidence in mirtarbaseEvidences:

    # 		Species (miRNA)	Target Gene	Target Gene (Entrez Gene ID)	Species (Target Gene)	Experiments	Support Type	References (PMID)


    mirtarID = mirnaEvidence['miRTarBase ID']
    mirtarMIRNA = mirnaEvidence['miRNA']
    mirtarMIRNASpecies = mirnaEvidence['Species (miRNA)']
    mirtarGENE = mirnaEvidence['Target Gene'].upper()
    mirtarGENESpecies = mirnaEvidence['Species (Target Gene)']
    mirtarExperiment = mirnaEvidence['Experiments']

    mirtarExperiment = mirtarExperiment.split("/") if mirtarExperiment != None else []
    mirtarExperimentNew = []

    for x in [y.split(";") for y in mirtarExperiment if len(y) > 0]:
        for elem in x:
            mirtarExperimentNew.append( mirtarbase_exp_type(elem) )

    mirtarExperiment = mirtarExperimentNew

    mirtarSupport = mirtarbase_function_label(mirnaEvidence['Support Type'])
    mirtarRefs = int(mirnaEvidence['References (PMID)'])

    supportTypes[mirtarSupport] += 1

    for exp in mirtarExperiment:
        experimentTypes[exp] += 1

    mirnaSpeciesID = speciesName2TaxID(mirtarMIRNASpecies)
    geneSpeciesID = speciesName2TaxID(mirtarGENESpecies)

    if mirnaSpeciesID == geneSpeciesID:
        commonTaxID = mirnaSpeciesID
    else:
        commonTaxID = None

    db.createNodeIfNotExists(['PUBMED', 'EVIDENCE'], {'id': mirtarRefs})
    db.createNodeIfNotExists(['MIRTARBASE', 'EVIDENCE'], {'id': mirtarID, 'tax_gene': geneSpeciesID, 'tax_mirna': mirnaSpeciesID})

    # TODO add relation props?
    db.createRelationship('mtb', ['MIRTARBASE'], 'pb', {'id': mirtarID}, ['PUBMED'], {'id': mirtarRefs}, ['EVIDENCE_SUPPORT', 'LITERATURE_SUPPORT', mirtarSupport] + mirtarExperiment, {})

    if commonTaxID != None:
        db.createRelationship('mtb', ['MIRTARBASE'], 'taxid', ['TAX'], {'id': mirtarID}, {'id': mirnaSpeciesID},
                              ['ORGANISM_SUPPORT'], {})

    db.createRelationshipIfNotExists('gene', ['GENE'], {'id': mirtarGENE}, 'mtb', ['MIRTARBASE'], {'id': mirtarID}, ['GENE_MENTION'], {'tax': geneSpeciesID})
    db.createRelationshipIfNotExists('mtb', ['MIRTARBASE'], {'id': mirtarID}, 'mirna', ['MIRNA'], {'name': mirtarMIRNA}, ['MIRNA_MENTION'], {'tax': mirnaSpeciesID})

print(referencesWithComma)
print(supportTypes)
print(experimentTypes)

db.close()