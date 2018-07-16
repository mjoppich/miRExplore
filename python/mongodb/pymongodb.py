import json
import re
import sys
from bson import json_util
import pymongo
from pymongo import MongoClient
from pymongo import cursor


class MongoDB:
    def __init__(self, dbname):
        '''
        Constructs a new MongoClient and establishes the DB connection
        :param dbname: string
        '''
        self._conn = pymongo.MongoClient('127.0.0.1', 27017)
        self._db = self._conn[dbname]

    def createCollection(self, collection_name):
        '''
        Creates a new collection <collection_name> in the database
        :param collection_name: string
        '''
        self._db = self._db.create_collection(collection_name)

    def dropCollection(self, collection_name):
        '''
        Deletes an existing collection <collection_name> from the database
        :param collection_name: string
        '''
        self._db = self._db.drop_collection(collection_name)

    def getCollection(self, collection_name):
        '''
        Returns collection object for given <collection_name>
        :param collection_name: string
        :return: collection object
        '''
        return self._db.get_collection(collection_name)

    def countDictionariesInCollectionForValue(self, collection, key, value):
        '''
        Returns count of dictionaries within a collection, based on given search criteria
        :param collection: object
        :param key: string
        :param value: regex
        :return: int
        '''
        return collection.count({key: {'$regex': value}})

    def countDictionariesInCollectionAll(self, collection):
        '''
        Returns count of all dictionaries within a collection
        :param collection: object
        :return: int
        '''
        return collection.count()

    def findDictionariesInCollectionForValue(self, collection, key, value):
        '''
        Returns all dictionaries within a collection
        :param collection: object
        :param key: string
        :param value: regex
        :return: cursor object
        '''
        return collection.find({key: {'$regex': value}})

    def findDictionariesInCollectionAll(self, collection):
        '''
        Returns all dictionaries within a collection
        :param collection: object
        :return: cursor object
        '''
        return collection.find({})

    def insertJSONIntoCollection(self, collection, json_file):
        '''
        Appends all dictionaries (entries) from a json file to an existing collection.
        The json file should be an array of objects:
        [
            { name: "ex1", msg: "example 1" },
            { name: "ex2", msg: "example 2" },
            { name: "ex3", msg: "example 3" }
        ]
        :param collection: object
        :param json_file: bson file (compatible with pymongo)
        :return: insert object
        '''
        return collection.insert(json_file)

    def insertDictionaryIntoCollection(self, collection, dic):
        '''
        Appends a single dictionary (entry) from a json file to an existing collection
        :param collection: object
        :param dic: dictionary (single entry from json file)
        :return: insert_one object
        '''
        return collection.insert_one(dic)

    def getJsonObjectForGeneID(self, collection, geneID):
        '''
        Given a collection name and a gene_id, return the Json object containing interactions corresponding to the gene_id
        :param collection: object
        :param geneID: string (example: LNC_GE_mm10_00500789)
        '''
        docs_list  = list(collection.find({'gene_id':geneID}))
        return json.dumps(docs_list, default=json_util.default)

