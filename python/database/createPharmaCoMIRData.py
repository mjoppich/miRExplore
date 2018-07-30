import csv

from collections import defaultdict

from database.Neo4JInterface import neo4jInterface
from utils.idutils import dataDir

inputFile = dataDir + "/miRExplore/pharmacomir/pharmacomir_VERSE_DB.csv"

db = neo4jInterface(simulate=False, printQueries=True)

db.deleteRelationship('n', ['PUBMED'], None, 'm', ['PHARMACOMIR'], None, ['PHARMACOMIR_PUBMED_SUPPORT'], None)
db.deleteRelationship('n', ['GENE'], None, 'm', ['PHARMACOMIR'], None, ['PHARMACOMIR_GENE_MENTION'], None)
db.deleteRelationship('n', ['PHARMACOMIR'], None, 'm', ['MIRNAS'], None, ['PHARMACOMIR_MIRNA_MENTION'], None)

db.deleteNode(['PHARMACOMIR'], None)


hgncAlt2ID = {}

with open(dataDir+"/miRExplore/hgnc.tsv", 'r') as infile:

    hgncData = csv.DictReader(infile, delimiter='\t', quotechar='\"')
    for line in hgncData:

        sym = line['Approved Symbol']

        if line['Previous Symbols'] == None:
            continue

        altsym = line['Previous Symbols'].split(', ')
        altsyn = line['Synonyms'].split(', ') if 'Synonyms' in line and line['Synonyms'] != None else []

        altsym += altsyn

        altsym = [x.replace('-', '').upper() for x in altsym if len(x) > 0]

        for x in altsym:
            hgncAlt2ID[x] = sym


if 'RAS' in hgncAlt2ID:
    del hgncAlt2ID['RAS']
    #PUBLICATION WITHDRAWN !!


hgncAlt2ID['SERT'] = 'SLC6A4'
hgncAlt2ID['D3R'] = 'DRD3'
hgncAlt2ID['M(1)'] = 'CHRM1'
hgncAlt2ID['Spry2'] = 'SPRY2'
hgncAlt2ID['CDKN4'] = 'CDKN1B'


lines2add = defaultdict(set)

with open(inputFile, 'r') as infile:
    pharmaCoMIR = csv.DictReader(infile, delimiter=',', quotechar='\"')


    for line in pharmaCoMIR:

        mirnaID = line['miRNA']
        geneID = line['Gene']
        drugName = line['Drug']
        pmid = line['PubMed Id']


        geneExists = db.nodeExists(['GENE'], {'id': geneID})

        if geneID in hgncAlt2ID:
            geneExists = True
            geneID = hgncAlt2ID[geneID]

        if not geneExists:
            continue

        lines2add[ (mirnaID, geneID) ].add( (drugName, pmid) )


addedNode = 1
for relation in lines2add:

    mirnaID = relation[0]
    geneID = relation[1]

    assocs = lines2add[relation]
    drugList = list(set([x[0] for x in assocs]))

    nodeID = 'PCM' + str(addedNode)
    addedNode+= 1
    db.createNodeIfNotExists(['PHARMACOMIR', 'EVIDENCE'], {'id': nodeID, 'drugs': drugList})

    db.createRelationshipIfNotExists('n', ['GENE'], {'id': geneID}, 'm', ['PHARMACOMIR'], {'id': nodeID},
                                     ['PHARMACOMIR_GENE_MENTION'], None)

    db.runInDatabase(
        "MATCH (n:PHARMACOMIR), (m:MIRNAS) WHERE n.id='{pharmaid}' and m.name STARTS WITH 'hsa-' and m.name CONTAINS '{mirna}' CREATE (n)-[r:PHARMACOMIR_MIRNA_MENTION]->(m) RETURN r".format(
            pharmaid=nodeID, mirna=mirnaID))

    for assoc in assocs:

        drugname = assoc[0]
        pmid = assoc[1]

        db.createNodeIfNotExists(['PUBMED'], {'id': pmid})
        db.createRelationshipIfNotExists('n', ['PUBMED'], {'id': pmid}, 'm', ['PHARMACOMIR'], {'id': nodeID}, ['PHARMACOMIR_PUBMED_SUPPORT'], None)


db.close()