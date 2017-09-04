from database.Neo4JInterface import neo4jInterface
from pubmed.CoCitationStore import CoCitationStore
from utils.idutils import eprint, dataDir

db = neo4jInterface(simulate=False, printQueries=True)
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['PUBMED'], None, ['PUBMED_CITED_BY'], None)

createCitationLists = False

citationFile = dataDir + "/miRExplore/pubmed_citations.tsv"
citedByFile = dataDir + "/miRExplore/pubmed_citedby.tsv"

if createCitationLists:
    retVal = db.matchNodes(['PUBMED'], None, nodename='n')

    pmids = set()

    for x in retVal:
        nodeData = x['n']
        if 'id' in nodeData.properties:
            pmid = nodeData.properties['id']
            pmids.add(pmid)

        else:
            eprint("No data in: ", str(nodeData))

    print(len(pmids))

    store = CoCitationStore()
    foundCitations = store.getCites(pmids)
    foundCitedBy = store.getCitedBy(pmids)


    with open(citationFile, 'w') as outfile:
        for pmid in foundCitations:
            outfile.write(str(pmid) + "\t" + str(",".join([str(x) for x in foundCitations[pmid]])) + "\n")

    with open(citedByFile, 'w') as outfile:
        for pmid in foundCitedBy:
            outfile.write(str(pmid) + "\t" + str(",".join([x for x in foundCitedBy[pmid]])) + "\n")

with open(citationFile, 'r') as infile:

    for line in infile:

        if len(line.strip()) == 0:
            continue

        aline = [x.strip() for x in line.split('\t')]

        if len(aline) < 2 or aline[1]== '':
            continue

        targetPMID = aline[0]
        sourcePMIDs = set([x.strip() for x in aline[1].split(",")])

        for pmid in sourcePMIDs:
            db.createRelationship('source', ['PUBMED'], {'id': pmid}, 'target', ['PUBMED'], {'id': targetPMID}, ['PUBMED_CITED_BY'], None)

with open(citedByFile, 'r') as infile:

    for line in infile:

        if len(line.strip()) == 0:
            continue

        aline = [x.strip() for x in line.split('\t')]

        if len(aline) < 2 or aline[1]== '':
            continue

        sourcePMID = aline[0]
        targetPMIDs = set([x.strip() for x in aline[1].split(",")])

        for pmid in targetPMIDs:
            db.createRelationship('source', ['PUBMED'], {'id': targetPMID}, 'target', ['PUBMED'], {'id': pmid}, ['PUBMED_CITED_BY'], None)