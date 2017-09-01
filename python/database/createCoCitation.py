from database.Neo4JInterface import neo4jInterface
from pubmed.CoCitationStore import CoCitationStore
from utils.idutils import eprint, dataDir

db = neo4jInterface(simulate=False, printQueries=True)

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


with open(dataDir + "/miRExplore/pubmed_citations.tsv", 'w') as outfile:
    for pmid in foundCitations:
        outfile.write(str(pmid) + "\t" + str(",".join([str(x) for x in foundCitations[pmid]])) + "\n")

with open(dataDir + "/miRExplore/pubmed_citedby.tsv", 'w') as outfile:
    for pmid in foundCitedBy:
        outfile.write(str(pmid) + "\t" + str(",".join([x for x in foundCitedBy[pmid]])) + "\n")