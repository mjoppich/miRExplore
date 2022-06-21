import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import defaultdict
from pymongo import MongoClient
from pprint import pprint

class PMID2PMCDBMongo:


    def __init__(self, databaseName, mainDatabase="miRExplore"):

        self.databaseName = databaseName
        self.mainDatabase = mainDatabase
        self.client=MongoClient()
        self.db =self.client[mainDatabase]
        self.tables = []

        table = None
        
        if databaseName in self.db.list_collection_names():
            table = self.db[databaseName]

        if not table is None:
            self.tables.append(table)
        
        self.pmid2pmc = defaultdict(set)
        self.pmc2pmid = defaultdict(set)

    def has_database(self):
        return self.databaseName in self.db.list_collection_names()

    def create_database(self):
        self.tables.append( self.db[self.databaseName] )

    def add_database(self, obase):

        assert(isinstance(obase, PMID2XDBMongo))

        for x in obase.tables:
            self.tables.append(x)

        for x in obase.all_term_names:
            if not x in self.all_term_names:
                self.all_term_names.append(x)
        #self.all_term_names += obase.all_term_names

    def insert_into_database(self, data):
        self.tables[0].insert_one(data)


    def hasID(self, docid):

        # TODO: search string and number here!
        for table in self.tables:
            for key in ["pmid", "pmc"]:
                res=table.find_one({key: docid})
                if len(res) > 0:
                    return True
        
        return False


    def getDOC(self, docid, default=None):

        res = []

        for table in self.tables:
            res+=table.find({"docid": docid})

        if len(res) == 0:
            return default
        return res


    def hasID(self, sid):

        if sid in self.pmid2pmc:
            return True

        if str(sid) in self.pmid2pmc:
            return True

        return False

    def getID(self, sid, default=None):

        if sid in self.pmid2pmc:
            return self.pmid2pmc[sid]

        if str(sid) in self.pmid2pmc:
            return self.pmid2pmc[str(sid)]

        return default


    def getAllPMIDs(self):
        return [x for x in self.pmid2pmc]

    def getAllPMCs(self):
        return [x for x in self.pmc2pmid]

    @classmethod
    def loadFromFile(cls, filepath, PMC2PMID=True):

        if databaseName is None:
            databaseName = os.path.basename(filepath).split(".")[0]
            print("Assigned databaseName", databaseName)

        ret = PMID2PMCDBMongo(databaseName)

        if not ret.has_database():           
            print("Creating new database", databaseName)
            ret.create_database()


        with open(filepath, 'r') as fin:

            for line in fin:

                aline = line.strip().split('\t')

                assert(len(aline) == 2)

                if PMC2PMID:
                    pmid = aline[1]
                    pmc = aline[0]
                else:
                    pmid = aline[0]
                    pmc = aline[1]

                assert(pmc.startswith("PMC"))

                info = {"PMID": str(pmid), "PMC": str(pmc)}

                ret.insert_into_database(info)


        return ret