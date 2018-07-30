import glob
from collections import defaultdict

from porestat.utils.Parallel import MapReduce

from database.Neo4JInterface import neo4jInterface
from pubmed.CoCitationStore import CoCitationStore
from utils.idutils import eprint, dataDir

db = neo4jInterface(simulate=False, printQueries=True)
db.deleteRelationship('n', ['PUBMED'], None, 'm', ['PUBMED'], None, ['PUBMED_CITED_BY'], None, 'r')
retVal = db.matchNodes(['PUBMED'], None, nodename='n')

pmids = set()

for x in retVal:
    if len(pmids) == -1:
        break

    nodeData = x['n']
    if 'id' in nodeData.properties:
        pmid = nodeData.properties['id']
        pmids.add(pmid)

    else:
        eprint("No data in: ", str(nodeData))

if len(pmids) == 0:
    eprint("No RELEVANT PUBMED entries found")
else:
    print("Relevant PMIDs found: ", len(pmids))

allfiles = glob.glob('/local/storage/pubmed/*.citation')

db.close()

def analyseFile(file, envPMIDs):

    allpmids = defaultdict(set)

    print("Starting file: ", file)

    procDB = neo4jInterface(simulate=False, printQueries=True)


    with open(file, 'r') as infile:

        for line in infile:

            aline = line.strip().split('\t')

            pmid_cites = aline[0]
            pmid_cited_by = aline[1]

            if pmid_cites in envPMIDs and pmid_cited_by in envPMIDs:
                allpmids[pmid_cites].add(pmid_cited_by)

    if len(allpmids) > 0:

        for pmid in allpmids:

            for opmid in allpmids[pmid]:

                if not opmid in envPMIDs:
                    continue

                procDB.createRelationship('cpmid', ['PUBMED'], {'id': pmid}, 'opmid', ['PUBMED'], {'id': opmid}, ['PUBMED_CITED_BY'], None)

    procDB.close()

ll = MapReduce(6)
result = ll.exec( allfiles, analyseFile, pmids, 1, None)