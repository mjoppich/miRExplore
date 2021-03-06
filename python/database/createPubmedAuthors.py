import glob
from collections import defaultdict

from database.Neo4JInterface import neo4jInterface
from pubmed.CoCitationStore import CoCitationStore
from utils.idutils import eprint, dataDir

db = neo4jInterface(simulate=False, printQueries=False)
db.deleteRelationship('n', ['PUBMED_AUTHOR'], None, 'm', ['PUBMED'], None, ['IS_AUTHOR'], None, 'r')
db.deleteNode(["PUBMED_AUTHOR"], None)
db.createUniquenessConstraint('PUBMED_AUTHOR', 'id')

retVal = db.matchNodes(['PUBMED'], None, nodename='n')
pmids = set()
for x in retVal:

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

class PubmedAuthor:

    def __init__(self, firstname, lastname, initials):
        self.firstname = firstname
        self.lastname = lastname
        self.initials = initials

    def makeNodeID(self):
        base = self.firstname.upper() + "_" + self.lastname.upper()
        base = base.replace(" ", "_")
        return base

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "<author {fn} {ln} {init}/>".format(fn=self.firstname, ln=self.lastname, init=self.initials)

    def to_node_props(self):
        return {'firstname': self.firstname, 'lastname': self.lastname, 'id': self.makeNodeID()}

    def __eq__(self, other):
        return self.makeNodeID() == other.makeNodeID()

    def __hash__(self):
        return hash(self.makeNodeID())


allfiles = glob.glob('/local/storage/pubmed/*.author')
allpmids = defaultdict(set)

for file in allfiles:

    with open(file, 'r') as infile:

        for line in infile:

            aline = line.split('\t')

            pmid = aline[0]
            fname = aline[1]
            iname = aline[2]
            lname = aline[3]

            if pmid in pmids:
                allpmids[pmid].add( PubmedAuthor(fname, lname, iname) )

createdAuthors = set()
for pmid in allpmids:

    for author in allpmids[pmid]:

        if not author in createdAuthors:
            authProps = author.to_node_props()
            db.createNode(['PUBMED_AUTHOR'], authProps)

        authID = author.makeNodeID()
        db.createRelationship('author', ['PUBMED_AUTHOR'], {'id': authID}, 'pmid', ['PUBMED'], {'id': pmid}, ['IS_AUTHOR'], None)


"""
handle = Entrez.efetch(db="pubmed", id=",".join(pmids), retmode='XML')
record = Entrez.read(handle)

author2pmid = defaultdict(set)
authorid2author = {}

for article in record['PubmedArticle']:

    pmid = article['MedlineCitation']['PMID'] if 'PMID' in article['MedlineCitation'] else None

    if pmid == None:
        continue

    artInfo = article['MedlineCitation']['Article']
    articleTitle = artInfo['ArticleTitle']

    db.mergeNode(['PUBMED'], {'id': pmid}, {'title': articleTitle})

    if 'AuthorList' in artInfo and artInfo['AuthorList'] != None:
        for author in artInfo['AuthorList']:

            if not 'LastName' in author:
                continue

            lastname = author['LastName']
            firstname = author['ForeName'] if 'ForeName' in author else ''
            articleAuthorInitials = author['Initials'] if 'Initials' in author else ''

            auth = PubmedAuthor(firstname, lastname, articleAuthorInitials)

            author2pmid[auth].add(pmid)

            if auth.makeNodeID() in authorid2author:
                print(auth)
                print(authorid2author[auth.makeNodeID()])
            else:
                authorid2author[auth.makeNodeID()] = auth


seenAuthorIDs = set()
for x in author2pmid:
    if not x.makeNodeID() in seenAuthorIDs:
        seenAuthorIDs.add(x.makeNodeID())
    else:
        eprint("Author ID duplicate: " + x.makeNodeID())


id = 0
for author in author2pmid:
    authID = author.makeNodeID()

    authProps = author.to_node_props()
    db.createNode(['PUBMED_AUTHOR'], authProps)

    pmids = author2pmid[author]

    for pmid in pmids:
        db.createRelationship('author', ['PUBMED_AUTHOR'], {'id': authID}, 'pmid', ['PUBMED'], {'id': pmid}, ['IS_AUTHOR'], None)

db.close()

"""