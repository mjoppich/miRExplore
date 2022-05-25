import copy
import os, sys

from synonymes.SynonymFile import Synfile
from synonymes.mirnaID import miRNA
from textdb.DocOrganismDB import DocOrganismDB

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


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



    @property
    def ltype(self):
        return self.ltyped

    @property
    def rtype(self):
        return self.rtyped

    def types(self):
        return set([self.ltype, self.rtype])

    def getOntologyForType(self, ontType):

        if ontType == self.ltype:
            return self.lontology

        if ontType == self.rtype:
            return self.rontology

        return None

    @classmethod
    def choices(cls, elems, k=1, excluded = []):

        resvec = []

        while len(resvec) < k:

            relem = random.choice(elems)

            if not relem in resvec and not relem in excluded:
                resvec.append(relem)

        return resvec


    def get_rels(self, etype, eid):

        if etype == self.ltype:
            return self.get_lid_rels(eid)
        elif etype == self.rtype:
            return self.get_rid_rels(eid)

        return None



    def get_lid_rels(self, entID):

        if not self.l_ont_based:

            if self.ltype.upper() == "MIRNA":
                miFoundRels = set()
                for strMirna in self.ltype2rel:
                    oMirna = miRNA(strMirna)
                    if oMirna.accept(entID):
                        miFoundRels = miFoundRels.union(self.ltype2rel[strMirna])
                return miFoundRels

            else:
                return self.ltype2rel.get(entID, None)
        else:
            allRels = self.ltype2rel.get(entID, None)

            if allRels == None:
                return None

            rootObo = self.lontology.dTerms.get(entID, None)

            if rootObo != None:

                allChildren = rootObo.getAllChildren()

                lenBefore = len(allRels)
                for child in allChildren:
                    allRels = allRels.union(self.ltype2rel.get(child.term.id, []))

                print("Added", len(allRels)-lenBefore, "for ontology")


            return self.undoOntID(allRels)

    def get_rid_rels(self, entID):

        if not self.r_ont_based:

            if self.rtype.upper() == "MIRNA":
                miFoundRels = set()
                for strMirna in self.rtype2rel:
                    oMirna = miRNA(strMirna)
                    if oMirna.accept(entID):
                        miFoundRels = miFoundRels.union(self.rtype2rel[strMirna])
                return miFoundRels
            else:
                return self.ltype2rel.get(entID, None)

        else:
            allRels = self.rtype2rel.get(entID, None)

            if allRels == None:
                return None

            rootObo = self.rontology.dTerms.get(entID, None)

            if rootObo != None:

                allChildren = rootObo.getAllChildren()

                lenBefore = len(allRels)
                for child in allChildren:
                    allRels = allRels.union(self.rtype2rel.get(child.term.id, []))

                print("Added", len(allRels) - lenBefore, "for ontology")

            return self.undoOntID(allRels)


    def undoOntID(self, rels):

        for rel in rels:

            lid = rel.lid
            rid = rel.rid

            if self.l_ont_based:
                if lid in self.lontology.dTerms:
                    rel.lent = (self.lontology.dTerms[lid].name, rel.lent[1])

            if self.r_ont_based:
                if rid in self.rontology.dTerms:
                    rel.rent = (self.rontology.dTerms[rid].name, rel.rent[1])

        return rels



    @classmethod
    def loadFromFile(cls, filepath, ltype, rtype, dbtype='pmid', databaseName=None, normGeneSymbols=None, lontology=None, rontology=None, switchLR=False, lReplaceSc=True, rReplaceSc=True, ignoreDocIDs=None, stopAfter=-1, getDocs=False, excludeIDs=None):

        if databaseName is None:
            databaseName = os.path.basename(filepath).split(".")[0]
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
                                origlid = lid
                                #print(lid)

                            if lid == None or lid == "None":
                                continue

                            seenMirnas[origRid] = (org, lid)
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
                                #print(rid)

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

                        #rel.data_id = dataID
                        #ltype2rel/rtype2rel 
                        ret.ltype2rel[lid].add(rel)
                        ret.rtype2rel[rid].add(rel)


                    #all_ltypes und all_rtypes kann weg??
                    ret.all_ltypes.add(lid)
                    ret.all_rtypes.add(rid)

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

        return ret


if __name__ == '__main__':

    pmidBase = "/home/mjoppich/ownCloud/data/miRExplore/textmine/aggregated_pmid"
    MiGenRelDBMongo.loadFromFile(pmidBase + "/mirna_gene.spacy.pmid", ltype="gene", rtype="mirna")
