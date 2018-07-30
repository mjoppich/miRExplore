from collections import Counter
from neo4j.v1 import GraphDatabase, basic_auth
from porestat.utils.DataFrame import DataFrame
from utils.idutils import ltype2label, makeDBGeneID, dataDir
from database.Neo4JInterface import neo4jInterface

hgncGenes = DataFrame.parseFromFile(dataDir + "/miRExplore/hgnc_ensembl_entrez.tsv", bConvertTextToNumber=False)

allStatus = Counter()

db = neo4jInterface(simulate=False)
db.createUniquenessConstraint('GENE', 'id')

db.deleteRelationship('n', None, None, 'm', None, None, ['HAS_GENE'], None, 'r')
db.deleteNode(["GENE"], None)


for gene in hgncGenes:

    hgncID = gene['HGNC ID']
    hgncSym = gene['Approved Symbol']

    hgncName = gene['Approved Name']
    hgncEnsembl = gene['Ensembl ID(supplied by Ensembl)']
    hgncEntrez = gene['Entrez Gene ID(supplied by NCBI)']

    hgncStatus = gene['Status']
    hgncLocusType = gene['Locus Type']

    hgncDBLType = [ltype2label[ hgncLocusType ]] if hgncLocusType in ltype2label else []

    for x in hgncDBLType:
        allStatus[x] += 1

    if hgncStatus != 'Approved':
        continue

    dbID = makeDBGeneID(hgncSym)
    db.createNodeIfNotExists(['GENE']+hgncDBLType, {'id': dbID, 'name':hgncName})

    orgProps = {}
    if hgncEnsembl != None and len(hgncEnsembl) > 0:
        orgProps['ensembl'] = hgncEnsembl
    if hgncEntrez != None and len(hgncEntrez) > 0:
        orgProps['entrez'] = int(hgncEntrez)

    #print(hgncID + " " + hgncSym + " " + hgncEnsembl)
    if len(orgProps) > 0:
        db.createRelationshipIfNotExists("tax",['TAX'], {'id': 9606}, "gene", ['GENE'], {'id': dbID}, ['HAS_GENE'], orgProps)

db.close()