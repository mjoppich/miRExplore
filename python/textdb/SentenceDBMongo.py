import codecs
import glob

import os
from pymongo import MongoClient
from pprint import pprint

from collections import defaultdict, deque
import progressbar
import sys


class SentenceDBMongo:

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

        assert(isinstance(obase, SentenceDBMongo))

        for x in obase.tables:
            self.tables.append(x)


    def get_sentence(self, docid):

        res = []
        for table in self.tables:
            res+=[table.find_one({"sentid": docid})]

        if len(res) == 0:
            return None
        
        res = res[0]
        
        return (res["sentid"], res["sentence"])

    def insert_into_database(self, data):
        self.tables[0].insert_one(data)

    @classmethod
    def loadFromFile(cls, basepath, databaseName, requiredDocuments=None, dbPrefix=None):

        
        if dbPrefix != None:
            databaseName = "{}_{}".format(dbPrefix, databaseName)

        print("Assigned databaseName", databaseName)

        ret = SentenceDBMongo(databaseName)

        if not ret.has_database():
            
            print("Creating new database", databaseName)
            ret.create_database()

            for filename in glob.glob(basepath + "*.sent"):

                print(filename)

                with codecs.open(filename, 'r') as infile:

                    for line in infile:
                        #line = line.decode('latin1')
                        line = line.strip()
                        if len(line) == 0:
                            continue

                        aline = line.split("\t")

                        if len(aline) != 2:
                            #print("Invalid input line for sent:", line, file=sys.stderr)
                            continue

                        sentID = aline[0]
                        sentText = aline[1]
                        docID = sentID.split(".")[0]

                        if not docID in requiredDocuments:
                            continue

                        info = {"docid": docID, "sentid": sentID, "sentence": sentText}
                        ret.insert_into_database(info)
        
        for table in ret.tables:
            usefulIndices = ["docid", "sentid"]
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

    sentDB = SentenceDBMongo.loadFromFile(sentenceLocation, dbPrefix="pmid", databaseName="sentences", requiredDocuments=relevantDocs)

    senttxt = sentDB.get_sentence("34879371.1.1")
    print(senttxt)

    senttxt = sentDB.get_sentence("34879371.2.9")
    print(senttxt)


