from collections import defaultdict


class DocOrganismDB:

    def __init__(self):

        self.docid2orgs = defaultdict(set)


    def getDocOrgs(self, docid):

        return self.docid2orgs.get(docid, None)

    @classmethod
    def loadFromFile(cls, filepath):

        retDB = DocOrganismDB()

        with open(filepath, 'r') as fin:

            for line in fin:

                aline = line.strip().split()

                docid = aline[0]
                orgs = set(aline[1].split(","))

                for x in orgs:
                    retDB.docid2orgs[docid].add(x)


        return retDB

