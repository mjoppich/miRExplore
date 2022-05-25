from collections import defaultdict
from pymongo import MongoClient
from pprint import pprint
import os

class PMID2XDBMongo:


    def __init__(self, assocObo, databaseName, mainDatabase="miRExplore"):

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
        

        self.all_term_names = None

        if assocObo != None:
            self.all_term_names = []
            for termid in assocObo.dTerms:
                oterm = assocObo.dTerms[termid]
                self.all_term_names.append((oterm.name, oterm.id))


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


    def hasDOC(self, docid):

        # TODO: search string and number here!
        res=self.table.find_one({"docid": docid})
        if len(res) > 0:
            return True


    def getDOC(self, docid, default=None):
        res=self.table.find_all({"docid": docid})

        if len(res) == 0:
            return default
        return res

    def getTerms(self):
        return self.all_term_names

    @classmethod
    def loadFromFile(cls, filepath, assocObo, databaseName=None, reqDocIDs=None):

        if databaseName is None:
            databaseName = os.path.basename(filepath).split(".")[0]
            print("Assigned databaseName", databaseName)

        ret = PMID2XDBMongo(assocObo, databaseName)

        if not ret.has_database():
            
            print("Creating new database", databaseName)
            ret.create_database()

            with open(filepath, 'r') as fin:

                # 29113155        DOID:1389       polyneuropathy  [('29113155.2.1', 0, 14), ('29113155.2.2', 139, 153)]
                for line in fin:

                    line = line.strip()
                    aline = line.split('\t')
                    pmid = aline[0]

                    if reqDocIDs != None and not pmid in reqDocIDs:
                        continue

                    termid = aline[1]
                    termname = aline[2]

                    loc = eval(aline[3])

                    info = {
                        'docid': pmid,
                        'termid': termid,
                        'termname': termname,
                        'evidences': loc
                    }

                    ret.insert_into_database(info)


        return ret

