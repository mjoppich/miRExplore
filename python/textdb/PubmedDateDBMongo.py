import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

import glob
from textmining.SentenceID import SentenceID
from collections import defaultdict
import datetime

from pymongo import MongoClient
from pprint import pprint

class PubmedDateDBMongo:


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

        assert(isinstance(obase, PubmedDateDBMongo))

        for x in obase.tables:
            self.tables.append(x)

    def insert_into_database(self, data):
        self.tables[0].insert_one(data)


    def get_document_date(self, docid):

        docDate = self.get_document(docid)

        if docDate == None:
            return None

        ddate = [docDate["year"], docDate["month"], docDate["day"]]
    
        if ddate[0] == 0:
            return None
            
        if ddate[1] == 0:
            ddate[1] = 1
        if ddate[2] == 0:
            ddate[2] = 1
    
        dtDate = datetime.datetime.strptime("{}-{}-{}".format(*ddate), '%Y-%m-%d')

        return dtDate

    def get_document_timestamp(self, docid):

        dt = self.get_document_date(docid)

        if dt == None:
            return None

        return datetime.datetime.timestamp(dt)

    def get_document(self, docid):

        docid = str(docid)

        for table in self.tables:
            fResult = table.find_one({"docid": docid})

            if not fResult is None:
                return fResult

        return None



    @classmethod
    def loadFromFile(cls, basepath, databaseName, requiredDocuments=None, dbPrefix=None):

        if dbPrefix != None:
            databaseName = "{}_{}".format(dbPrefix, databaseName)

        print("Assigned databaseName", databaseName)

        ret = PubmedDateDBMongo(databaseName)

        if not ret.has_database():
            
            print("Creating new database", databaseName)
            ret.create_database()

            for filename in glob.glob(basepath + "*.date"):

                with open(filename, 'r') as fin:

                    print(filename)

                    for line in fin:

                        #14172228        1964    6       0
                        line = line.strip().split('\t')

                        docID = line[0]

                        if not docID in requiredDocuments:
                            continue

                        year = int(line[1])
                        month = int(line[2])
                        day = int(line[3])

                        ret.insert_into_database({"docid":docID, "year":year, "month":month, "day":day})


        print("Loading Dates Finished", file=sys.stderr)

        for table in ret.tables:
            usefulIndices = ["docid", "year"]
            print("Creating indices")
            for idx in usefulIndices:
                print("Creating index", idx)
                table.create_index(idx)

        return ret


if __name__ == '__main__':

    sentenceLocation = "/mnt/w/PubMed/"

    relevantDocsPath = "/mnt/w/miRExplore_pmid_pmc/aggregated_pmid/relevant_pmids.list"
    relevantDocs = set()
    with open(relevantDocsPath) as fin:
        for line in fin:
            line = line.strip()
            relevantDocs.add(line)

    dateDB = PubmedDateDBMongo.loadFromFile(sentenceLocation, databaseName="dates", dbPrefix="pmid", requiredDocuments=relevantDocs)
    date = dateDB.get_document("16141076") # 1979 7 0
    print(date)

    print(dateDB.get_document("test"))


