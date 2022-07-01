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

    def has_database(self):
        return self.databaseName in self.db.list_collection_names()

    def create_database(self):
        self.tables.append( self.db[self.databaseName] )

    def add_database(self, obase):

        assert(isinstance(obase, PMID2PMCDBMongo))

        for x in obase.tables:
            self.tables.append(x)

        for x in obase.all_term_names:
            if not x in self.all_term_names:
                self.all_term_names.append(x)

    def insert_into_database(self, data):
        self.tables[0].insert_one(data)


    def hasID(self, docid):

        # TODO: search string and number here!
        for table in self.tables:
            for key in ["pmid", "pmc"]:
                res=table.find_one({key: str(docid)})
                if len(res) > 0:
                    return True
        
        return False


    def getDOC(self, docid, default=None):

        res = []

        for table in self.tables:
            res+=table.find({"docid": str(docid)})

        if len(res) == 0:
            return default

        return res

    def getID(self, docid, idtype="pmid", default=None):

        returnIDs = set()
        # TODO: search string and number here!
        for table in self.tables:
            for key in ["pmid", "pmc"]:
                res=table.find({key: str(docid)})
                
                for x in res:
                    returnIDs.add(x[idtype])
                                
        if len(returnIDs) > 0:
            return returnIDs

        return default


    def getAllPMIDs(self):
        
        allDocIDs = set()

        for table in self.tables:
            print("find in table")
            res = table.find({},{ "_id": 0})
            print("adding to set")
            allDocIDs = allDocIDs.union([x["pmid"] for x in res])

        return allDocIDs

    def getAllPMCs(self):
        allDocIDs = set()

        for table in self.tables:
            res = table.find({},{ "_id": 0})
            allDocIDs = allDocIDs.union([x["pmid"] for x in res])
            
        return allDocIDs


    @classmethod
    def loadFromFile(cls, filepath, PMC2PMID=True, databaseName=None):

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

                    info = {"pmid": str(pmid), "pmc": str(pmc)}

                    ret.insert_into_database(info)

            for table in ret.tables:
                usefulIndices = ["pmid", "pmc"]
                print("Creating indices")
                for idx in usefulIndices:
                    print("Creating index", idx)
                    table.create_index(idx)

        return ret


if __name__ == "__main__":

    mDB = PMID2PMCDBMongo.loadFromFile('/mnt/w/miRExplore_pmid_pmc/aggregated_pmc/pmc2pmid', PMC2PMID=True)
    print(len(mDB.getAllPMIDs()))