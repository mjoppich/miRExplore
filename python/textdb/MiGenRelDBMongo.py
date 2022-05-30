import copy
import os, sys

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from synonymes.SynonymFile import Synfile
from synonymes.mirnaID import miRNA
from textdb.DocOrganismDB import DocOrganismDB

from pymongo import MongoClient
from pprint import pprint

import random
from collections import defaultdict

import os

from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor


class MiRGeneRel(DataBaseRel):

    def __init__(self, intuple, docid, sameParagraph, sameSentence, lent, rent, datasource, dataID):

        self.assocCat = None

        self.assocFound = None
        self.assocSent = None
        self.assocNeg = None
        self.assocPass = None
        self.lPOS = None
        self.rPOS = None
        self.assocPos = None

        self.lent=lent
        self.rent=rent

        self.lontent = lent
        self.rontent = rent

        self.data_source=datasource
        self.data_id = dataID
        self.trustv = None

        self.orig_names = None


        if intuple != None:
            #('12', '1V2', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (6, 16))
            #self.assocDir = None #intuple[0]
            #self.assocDirV = None #intuple[1]
            self.assocCat = None #intuple[2]
            self.assocInt = ""


            self.assocFound = intuple[3]
            self.assocSent = intuple[4]
            self.assocNeg = intuple[5]
            self.rPOS = intuple[6]
            self.lPOS = intuple[7]
            self.assocPos = intuple[8]

        self.pubID = docid
        self.same_sentence = sameSentence
        self.same_paragraph = sameParagraph

        self.orgs = None

    @property
    def lid(self):
        return self.lent[0]

    @property
    def rid(self):
        return self.rent[0]

    @property
    def ltype(self):
        return self.lent[1]

    @property
    def rtype(self):
        return self.rent[1]

    @property
    def docid(self):
        return self.pubID

    def __str__(self):
        return str(self.toJSON())

    def toJSON(self):

        trust = {}

        if self.trustv != None and len(self.trustv) == 4:
            trust['stack'] = self.trustv[-4]
            trust['verb'] = self.trustv[-3]
            trust['conj'] = self.trustv[-2]
            trust['relex'] = self.trustv[-1]

        retJSON =  {
            #'rel_direction': self.assocDir,
            #'rel_direction_verb': self.assocDirV,
            'rel_category': self.assocCat,
            'rel_interaction': self.assocInt,
            'rel_verb': self.assocFound,
            'rel_sentence': self.assocSent,
            'rel_negated': self.assocNeg,
            'rel_passive': self.assocPass,
            'lpos': self.lPOS,
            'rpos': self.rPOS,
            'rel_pos': self.assocPos,
            'same_paragraph': self.same_paragraph,
            'same_sentence': self.same_sentence,
            'docid': self.docid,
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,
            'lontid': self.lontent[0],
            'rontid': self.rontent[0],
            'data_source': self.data_source,
            'data_id': self.data_id,
            'trust': trust,
            'orig_ids': self.orig_names
        }

        if self.orgs != None:
            retJSON['orgs'] = tuple(self.orgs)

        return retJSON

from pymongo import MongoClient


class MiGenRelDBMongo(DataBaseDescriptor):


    def __init__(self, ltype, rtype, databaseName, mainDatabase="miRExplore"):

        super().__init__()

        self.ltyped = ltype
        self.rtyped = rtype

        self.databaseName = databaseName
        self.mainDatabase = mainDatabase
        self.client = MongoClient()
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

        assert(isinstance(obase, MiGenRelDBMongo))

        for x in obase.tables:
            self.tables.append(x)

        for x in obase.all_term_names:
            if not x in self.all_term_names:
                self.all_term_names.append(x)
        #self.all_term_names += obase.all_term_names

    def insert_into_database(self, data):
        self.tables[0].insert_one(data)

    def get_evidence_docids(self):

        docIDs = []
        for table in self.tables:
            allDocIDs = table.distinct("docid")
            docIDs += allDocIDs

        docIDs = set(docIDs)

        return docIDs

    @property
    def ltype(self):
        return self.ltyped

    @property
    def rtype(self):
        return self.rtyped

    def types(self):
        return set([self.ltype, self.rtype])


    def find_relations(self, ltypes, rtypes):

        ltypeGiven = ltypes != None and len(ltypes) > 0
        rtypeGiven = rtypes != None and len(rtypes) > 0

        filterQuery = {}
        if ltypeGiven and rtypeGiven:
            filterQuery = {"ltype": self.ltype, "lid": {"$in": ltypes}, "rtype": self.rtype, "rid": {"$in": rtypes}}
        elif ltypeGiven:
            filterQuery = {"ltype": self.ltype, "lid": {"$in": ltypes}}
        elif rtypeGiven:
            filterQuery = {"rtype": self.rtype, "rid": {"$in": rtypes}}
            

        allRes = []
        for table in self.tables:
            tableRes = table.find(filterQuery, {'_id': False}, sort=[("docid", -1)])

            for x in tableRes:
                x["lpos"] = tuple(x["lpos"])
                x["rpos"] = tuple(x["rpos"])
                x["rel_pos"] = tuple(x["rel_pos"])

                if x not in allRes:
                    allRes.append(x)

        return allRes


    @classmethod
    def loadFromFile(cls, filepath, ltype, rtype, dbtype='pmid', databaseName=None, dbPrefix=None, normGeneSymbols=None, lontology=None, rontology=None, switchLR=False, lReplaceSc=True, rReplaceSc=True, ignoreDocIDs=None, stopAfter=-1, getDocs=False, excludeIDs=None):

        if databaseName is None:
            databaseName = "_".join(os.path.basename(filepath).split(".")[0:-1])

        if dbPrefix != None:
            databaseName = "{}_{}".format(dbPrefix, databaseName)

        print("Assigned databaseName", databaseName)


        if switchLR:
            ontTmp = lontology
            typeTmp = ltype

            lontology = rontology
            ltype = rtype

            rontology = ontTmp
            rtype = typeTmp


        ret = MiGenRelDBMongo(ltype, rtype, databaseName)
        ret.lontology = lontology
        ret.rontology = rontology


        if not ret.has_database():
            
            print("Creating new database", databaseName)
            ret.create_database()


            file_base = os.path.basename(filepath)
            fileDir = os.path.dirname(filepath)


            docOrgDB = None
            if os.path.exists(fileDir + "/organism."+dbtype):
                docOrgDB = DocOrganismDB.loadFromFile(fileDir + "/organism."+dbtype)

            seenRels = set()

            geneSymbolsNormalized = 0

            seenMirnas = {}
            seenHarmMirnas = set()
            seenGenes = set()

            ignoredDocs = set()
            takenDocs = set()

            addedCount = 0

            with open(filepath, 'r') as fin:


                #miR-182 mir-182 MIRNA   CBFA3   RUNX3   GENE    29054094        True    True    [('21', 'V21', 'POS', 'express', '29054094.2.3', False, (92, 99), (82, 87), (62, 72), 'all_rels')]
                #MIR_497	miR-497	MIRNA	INSR	insulin receptor	GENE	26300412	True	True	[('21', '2V1', 'DOWN', '', '26300412.2.8', False, (39, 46), (166, 182), (0, 0), 'spacy', 1, 1, 0, 0, True, False, False, 'MIR_GENE', 'DOWN')]

                for lineIdx, line in enumerate(fin):

                    if stopAfter != -1:
                        if addedCount >= stopAfter:
                            break

                    line = line.strip()
                    aline = line.split('\t')

                    turnEvs= False
                    docid = aline[6]

                    if excludeIDs != None:
                        if docid in excludeIDs:
                            continue


                    if switchLR:
                        turnEvs = True

                        tmp = aline[0:3]

                        if len(aline) < 4:
                            print(aline)

                        aline[0] = aline[3]
                        aline[1] = aline[4]
                        aline[2] = aline[5]

                        aline[3] = tmp[0]
                        aline[4] = tmp[1]
                        aline[5] = tmp[2]

                    lIDIdx = 0
                    rIDIdx = 3

                    if ltype == "mirna":
                        lIDIdx += 1

                    if rtype == "mirna":
                        rIDIdx += 1

                    lid = aline[lIDIdx]
                    rid = aline[rIDIdx]

                    olid = lid
                    orid = rid

                    if ret.l_ont_based and lReplaceSc:
                        lid = lid.replace('_', ':', 1)

                    if ret.r_ont_based and rReplaceSc:
                        rid = rid.replace('_', ':', 1)

                    if ltype == rtype and lid > rid:
                        tmp = lid
                        lid = rid
                        rid = tmp

                        tmp = olid
                        olid = orid
                        orid = tmp

                        turnEvs = not turnEvs

                    sameParagraph = aline[7] == 'True'
                    sameSentence = aline[8] == 'True'

                    if not sameSentence:
                        continue


                    if normGeneSymbols!= None:

                        if ltype == 'gene':
                            lid = lid.upper()
                            if lid in normGeneSymbols:
                                lid = normGeneSymbols[lid]
                                geneSymbolsNormalized += 1

                        elif rtype == 'gene':
                            rid = rid.upper()
                            if rid in normGeneSymbols:
                                rid = normGeneSymbols[rid]
                                geneSymbolsNormalized += 1

                    org = None
                    if ltype == 'mirna':

                        if lid in seenMirnas:
                            org, lid = seenMirnas[lid]
                        else:
                            origRid = lid
                            org, lid = cls.harmonizeMIRNA(lid)

                            
                            if lid == "None" and aline[lIDIdx-1] == "LET_7":
                                lid="let-7"
                                origLid = lid
                                print("NONE LET7 CASE", rid)

                            if lid == None or lid == "None":
                                continue

                            seenMirnas[origLid] = (org, lid)
                            seenHarmMirnas.add(lid)
                        
                        #lid = olid #WHY???

                    elif rtype == 'mirna':

                        if rid in seenMirnas:
                            org, rid = seenMirnas[rid]
                        else:
                            origRid = rid
                            org, rid = cls.harmonizeMIRNA(rid)

                            if rid == "None" and aline[rIDIdx-1] == "LET_7":
                                rid="let-7"
                                origRid = rid
                                print("NONE LET7 CASE", rid)

                            if rid == None or rid == "None":
                                continue

                            seenMirnas[origRid] = (org, rid)
                            seenHarmMirnas.add(rid)

                        #rid = orid # WHYYY?


                    if ltype == "gene":
                        seenGenes.add(lid)
                    elif rtype == "gene":
                        seenGenes.add(rid)

                    docOrgs = set()
                    if org != None:
                        docOrgs.add(org)


                    if ignoreDocIDs != None and docid in ignoreDocIDs:
                        ignoredDocs.add(docid)
                        continue

                    takenDocs.add(docid)

                    docOrgRes = None
                    
                    if docOrgDB != None:
                        docOrgRes = docOrgDB.getDocOrgs(docid)

                    if docOrgRes != None:
                        docOrgs = docOrgs.union(docOrgRes)


                    tmRelations = None if aline[9] == 'None' else eval(aline[9])


                    if ltype == rtype and lid > rid:
                        tmp = lid
                        lid = rid
                        rid = tmp

                    if tmRelations != None:
                        allrels = set()

                        for relIdx, rel in enumerate(tmRelations):

                            newrel = MiRGeneRel(rel, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), dbtype, None)
                            newrel.orig_names = (olid, orid)

                            #('21', '2V1', 'DOWN', '', '26300412.2.8', False, (39, 46), (166, 182), (0, 0), 'spacy', 1, 1, 0, 0, True, False, False, 'MIR_GENE', 'DOWN')

                            if len(rel) > 10:
                                stackEv = rel[10]
                                verbEv = rel[11]

                                if not stackEv and not verbEv:
                                    continue
                            
                                """
                                if ltype == "mirna":
                                    if rel[17] == "MIR_GENE":
                                        newrel.assocDir = "12"
                                        newrel.assocDirV = "1V2"
                                    else:
                                        newrel.assocDir = "21"
                                        newrel.assocDirV = "2V1"

                                else:
                                
                                    if rel[17] == "GENE_MIR":
                                        newrel.assocDir = "12"
                                        newrel.assocDirV = "1V2"
                                    else:
                                        newrel.assocDir = "21"
                                        newrel.assocDirV = "2V1"
                                """

                                newrel.assocCat = rel[18]
                                newrel.assocInt = rel[17]

                                newrel.assocPass = rel[15]
                                newrel.assocNeg = rel[16]


                            """
                            if turnEvs:

                                tmp = newrel.lPOS
                                newrel.lPOS = newrel.rPOS
                                newrel.rPOS = tmp

                                if newrel.assocDirV != None:
                                    nAssocDirV =""
                                    subs = {'1': '2', '2': '1', 'V': 'V'}
                                    for x in newrel.assocDirV:
                                        nAssocDirV += subs.get(x, x)

                                    newrel.assocDirV = nAssocDirV
                                    newrel.assocDir = nAssocDirV.replace('V', '')
                            """

                            if isinstance(stackEv, int):
                                newrel.trustv = (rel[10], rel[11], rel[12], rel[13])

                            relCopy = copy.deepcopy(newrel)

                            dataID = file_base + "_" + str(lineIdx) + "_" + str(relIdx)

                            if relCopy in seenRels:
                                continue

                            seenRels.add(relCopy)

                            newrel.data_id = dataID
                            allrels.add(newrel)

                        relations = allrels
                    else:
                        relations = set([MiRGeneRel(None, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), dbtype, None)])


                    if len(relations) == 0:
                        continue

                    for relIdx, rel in enumerate(relations):
                        if len(docOrgs) > 0:
                            rel.orgs = tuple(docOrgs)

                        ret.insert_into_database(rel.toJSON())


                    addedCount += 1

            for x in ret.ltype2rel:
                for elem in ret.ltype2rel[x]:
                    if elem.data_id == None:
                        print("Elem with no data id!")

            print("Gene Symbols Normalized", geneSymbolsNormalized)
            print("Loaded file", filepath)
            print("Accepted Doc IDs", len(takenDocs))
            print("Rejected Doc IDs", len(ignoredDocs))
            print("Seen genes", len(seenGenes))
            print("Seen miRNAs", len(seenMirnas))
            print("Seen Harm. miRNAs", len(seenHarmMirnas))

        if getDocs:
            return ret, ret.all_documnets()

        """
          {
            _id: ObjectId("6290ce44580ba7435e5856a7"),
            rel_category: 'NEU',
            rel_interaction: 'MIR_GENE',
            rel_verb: '',
            rel_sentence: '34858199.2.11',
            rel_negated: false,
            rel_passive: false,
            lpos: [ 30, 35 ],
            rpos: [ 84, 91 ],
            rel_pos: [ 0, 0 ],
            same_paragraph: true,
            same_sentence: true,
            docid: '34858199',
            ltype: 'gene',
            rtype: 'mirna',
            rid: 'miR-204',
            lid: 'NLRP3',
            lontid: 'NLRP3',
            rontid: 'miR-204',
            data_source: 'pmid',
            data_id: 'mirna_gene.hsa.pmid_4397_0',
            trust: { stack: 1, verb: 1, conj: 0, relex: 0 },
            orig_ids: [ 'NLRP3', 'miR-204' ],
            orgs: [ 'mmu', 'hsa' ]
        }

        """

        for table in ret.tables:
            usefulIndices = ["ltype", "rtype", "rid", "lid", "docid"]
            print("Creating indices")
            for idx in usefulIndices:
                print("Creating index", idx)
                table.create_index(idx)

        return ret


if __name__ == '__main__':

    pmidBase = "/mnt/w/miRExplore_pmid_pmc/aggregated_pmid/"
    MiGenRelDBMongo.loadFromFile(pmidBase + "/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene")