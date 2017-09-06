from collections import Counter, defaultdict
from porestat.utils.DataFrame import DataFrame
from utils.idutils import ltype2label, makeDBGeneID, mirtarbase_exp_type, mirtarbase_function_label, speciesName2TaxID, \
    dataDir
from database.Neo4JInterface import neo4jInterface

mirtarbaseEvidences = DataFrame.parseFromFile(dataDir + "/miRExplore/miRTarBase.csv", bConvertTextToNumber=False)

print(mirtarbaseEvidences.getHeader())

experimentTypes = Counter()
supportTypes = Counter()
referencesWithComma = Counter()

db = neo4jInterface(simulate=False, printQueries=False)
db.deleteRelationship('n', ['GENE'], None, 'm', ['MIRTARBASE'], None, ['GENE_MENTION'], None, 'r')
db.deleteRelationship('n', ['MIRTARBASE'], None, 'm', ['MIRNA'], None, ['MIRNA_MENTION'], None, 'r')
db.deleteRelationship('n', ['MIRTARBASE'], None, 'm', ['PUBMED'], None, ['MIRTARBASE_LITERATURE_SUPPORT'], None, 'r')
db.deleteRelationship('n', ['MIRTARBASE_SUPPORT'], None, 'm', ['MIRTARBASE'], None, ['MIRTARBASE_FUNCTIONAL_SUPPORT'], None, 'r')
db.deleteRelationship('n', ['MIRTARBASE_EXPERIMENT'], None, 'm', ['MIRTARBASE'], None, ['MIRTARBASE_EXPERIMENT_SUPPORT'], None, 'r')
db.deleteRelationship('n', ['MIRTARBASE'], None, 'm', ['TAX'], None, ['ORGANISM_SUPPORT'], None, 'r')

db.deleteNode(["MIRTARBASE"], None)
db.deleteNode(["MIRTARBASE_SUPPORT"], None)
db.deleteNode(["MIRTARBASE_EXPERIMENT"], None)
db.createUniquenessConstraint('MIRTARBASE', 'id')


dbcreatedExpTypes = set()
dbcreatedSupportTypes = set()
dbcreatedMIRTIDs = set()

dbcreatedMIRT2TAX = defaultdict(set)
dbcreatedMIRT2PUBMED = defaultdict(set)
dbcreatedMIRT2ExpTypes = defaultdict(set)
dbcreatedMIRT2SupportTypes = defaultdict(set)

dbcreatedMIRT2MIRNA = defaultdict(set)

print("Starting Evidences")
lastCheckedMIRT = 0

for mirnaEvidence in mirtarbaseEvidences:

    # 		Species (miRNA)	Target Gene	Target Gene (Entrez Gene ID)	Species (Target Gene)	Experiments	Support Type	References (PMID)


    mirtarID = mirnaEvidence['miRTarBase ID']

    try:
        if int(mirtarID.replace('MIRT', '')) > lastCheckedMIRT:
            lastCheckedMIRT += 1000
            print(mirtarID)
    except:
        pass

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
    mirtarRefs = mirnaEvidence['References (PMID)']

    supportTypes[mirtarSupport] += 1

    for exp in mirtarExperiment:
        experimentTypes[exp] += 1

    mirnaSpeciesID = speciesName2TaxID.get(mirtarMIRNASpecies, None)
    geneSpeciesID = speciesName2TaxID.get(mirtarGENESpecies, None)

    if mirnaSpeciesID == None and geneSpeciesID == None:
        continue

    if mirnaSpeciesID == geneSpeciesID:
        commonTaxID = mirnaSpeciesID
    else:
        commonTaxID = None

    db.createNodeIfNotExists(['PUBMED', 'EVIDENCE'], {'id': mirtarRefs})

    if not mirtarID in dbcreatedMIRTIDs:
        dbcreatedMIRTIDs.add(mirtarID)
        db.createNodeIfNotExists(['MIRTARBASE', 'EVIDENCE'], {'id': mirtarID, 'tax_gene': geneSpeciesID, 'tax_mirna': mirnaSpeciesID})

    if not mirtarSupport in dbcreatedSupportTypes:
        db.createNodeIfNotExists(['MIRTARBASE_SUPPORT'], {'id': mirtarSupport})
        dbcreatedSupportTypes.add(mirtarSupport)

    if not mirtarSupport in dbcreatedMIRT2SupportTypes[mirtarID]:
        dbcreatedMIRT2SupportTypes[mirtarID].add(mirtarSupport)
        db.createRelationship('ms', ['MIRTARBASE_SUPPORT'], {'id': mirtarSupport}, 'mtb', ['MIRTARBASE'], {'id': mirtarID}, ['MIRTARBASE_FUNCTIONAL_SUPPORT'], None)

    for expType in mirtarExperiment:

        if not expType in dbcreatedExpTypes:
            db.createNodeIfNotExists(['MIRTARBASE_EXPERIMENT'], {'id': expType})
            dbcreatedExpTypes.add(expType)

        if not expType in dbcreatedMIRT2ExpTypes[mirtarID]:
            dbcreatedMIRT2ExpTypes[mirtarID].add(expType)
            db.createRelationship('me', ['MIRTARBASE_EXPERIMENT'], {'id': expType}, 'mtb', ['MIRTARBASE'], {'id': mirtarID}, ['MIRTARBASE_EXPERIMENT_SUPPORT'], None)

    # TODO add relation props?
    if not mirtarRefs in dbcreatedMIRT2PUBMED[mirtarID]:
        dbcreatedMIRT2PUBMED[mirtarID].add(mirtarRefs)
        db.createRelationship('pb', ['PUBMED'], {'id': mirtarRefs}, 'mtb', ['MIRTARBASE'], {'id': mirtarID}, ['MIRTARBASE_LITERATURE_SUPPORT'], {})

    if commonTaxID != None:

        if not mirnaSpeciesID in dbcreatedMIRT2TAX[mirtarID]:
            dbcreatedMIRT2TAX[mirtarID].add(mirnaSpeciesID)
            db.createRelationship('mtb', ['MIRTARBASE'], {'id': mirtarID}, 'taxid', ['TAX'], {'id': mirnaSpeciesID}, ['ORGANISM_SUPPORT'], {})


    db.createRelationshipIfNotExists('gene', ['GENE'], {'id': mirtarGENE}, 'mtb', ['MIRTARBASE'], {'id': mirtarID}, ['GENE_MENTION'], {'tax': geneSpeciesID})

    if not mirtarMIRNA in dbcreatedMIRT2MIRNA[mirtarID]:
        dbcreatedMIRT2MIRNA[mirtarID].add(mirtarMIRNA)
        db.createRelationshipIfNotExists('mtb', ['MIRTARBASE'], {'id': mirtarID}, 'mirna', ['MIRNA'],
                                         {'name': mirtarMIRNA}, ['MIRNA_MENTION'], {'tax': mirnaSpeciesID})

print(referencesWithComma)
print(supportTypes)
print(experimentTypes)

db.close()