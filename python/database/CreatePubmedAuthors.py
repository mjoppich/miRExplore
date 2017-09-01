from Bio import Entrez

from database.Neo4JInterface import neo4jInterface
from pubmed.CoCitationStore import CoCitationStore
from utils.idutils import eprint, dataDir

db = neo4jInterface(simulate=False, printQueries=True)

retVal = db.matchNodes(['PUBMED'], None, nodename='n')

pmids = set()

for x in retVal:
    if len(pmids) == 10:
        break

    nodeData = x['n']
    if 'id' in nodeData.properties:
        pmid = nodeData.properties['id']
        pmids.add(pmid)

    else:
        eprint("No data in: ", str(nodeData))

handle = Entrez.efetch(db="pubmed", id=pmids, retmode='XML')
record = Entrez.read(handle)

for elem in record:
    print(elem)