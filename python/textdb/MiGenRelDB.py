import copy
import os, sys

from synonymes.SynonymFile import Synfile
from textdb.DocOrganismDB import DocOrganismDB

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import random
from collections import defaultdict

import os

from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor


class MiRGeneRel(DataBaseRel):

    def __init__(self, intuple, docid, sameParagraph, sameSentence, lent, rent, datasource, dataID):

        self.assocDir = None
        self.assocDirV = None
        self.assocCat = None
        self.assocFound = None
        self.assocSent = None
        self.assocNeg = None
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


        if intuple != None:
            #('12', '1V2', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (6, 16))
            self.assocDir = intuple[0]
            self.assocDirV = intuple[1]
            self.assocCat = intuple[2]
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


    def toJSON(self):

        trust = {}

        if self.trustv != None and len(self.trustv) == 4:
            trust['stack'] = self.trustv[-4]
            trust['verb'] = self.trustv[-3]
            trust['conj'] = self.trustv[-2]
            trust['relex'] = self.trustv[-1]

        retJSON =  {
            'rel_direction': self.assocDir,
            'rel_direction_verb': self.assocDirV,
            'rel_category': self.assocCat,
            'rel_verb': self.assocFound,
            'rel_sentence': self.assocSent,
            'rel_negated': self.assocNeg,
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
            'trust': trust
        }

        if self.orgs != None:
            retJSON['orgs'] = tuple(self.orgs)

        return retJSON

class MiGenRelDB(DataBaseDescriptor):


    def __init__(self, ltype, rtype, dbfilepath):

        super().__init__()

        self.ltyped = ltype
        self.rtyped = rtype

        self.dbfilepath = dbfilepath

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


    def get_lid_rels(self, geneID):

        if not self.l_ont_based:
            return self.ltype2rel.get(geneID, None)
        else:
            allRels = self.ltype2rel.get(geneID, None)

            if allRels == None:
                return None

            rootObo = self.lontology.dTerms.get(geneID, None)

            if rootObo != None:

                allChildren = rootObo.getAllChildren()

                lenBefore = len(allRels)
                for child in allChildren:
                    allRels = allRels.union(self.ltype2rel.get(child.term.id, []))

                print("Added", len(allRels)-lenBefore, "for ontology")


            return self.undoOntID(allRels)

    def get_rid_rels(self, geneID):

        if not self.r_ont_based:
            return self.rtype2rel.get(geneID, None)
        else:
            allRels = self.rtype2rel.get(geneID, None)

            if allRels == None:
                return None

            rootObo = self.rontology.dTerms.get(geneID, None)

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


    def toSQLite(self):

        if os.path.isfile(self.dbfilepath):
            os.remove(self.dbfilepath)

        conn, c = self.getConnCursor()

        c.execute("create table migenrel (lid text, ltype text, rid text, rtype text, docid text, evs text, sameparagraph bool, samesentence bool)")
        c.execute("create index migenrel_lid on migenrel (lid)")
        c.execute("create index migenrel_rid on migenrel (rid)")
        c.execute("create index migenrel_docid on migenrel (docid)")

        insertSQL = '''INSERT INTO migenrel(lid,ltype,rid, rtype, docid, evs, sameparagraph, samesentence) VALUES (?,?,?,?,?,?,?,?)'''

        for x in self.ltype2rel:
            rels = self.ltype2rel[x]

            for rel in rels:
                c.execute(insertSQL, (rel.lid, rel.ltype, rel.rid, rel.rtype, rel.pubID, rel.same_paragraph, rel.same_sentence, str(rel.toJSON()), ))

        conn.commit()
        c.close()


    @classmethod
    def loadFromFile(cls, filepath, ltype, rtype, dbtype='pmid', normGeneSymbols=None, lontology=None, rontology=None, switchLR=False, lReplaceSc=True, rReplaceSc=True, ignoreDocIDs=None):


        if switchLR:
            ontTmp = lontology
            typeTmp = ltype

            lontology = rontology
            ltype = rtype

            rontology = ontTmp
            rtype = typeTmp


        ret = MiGenRelDB(ltype, rtype, filepath + ".db")
        ret.lontology = lontology
        ret.rontology = rontology


        file_base = os.path.basename(filepath)
        fileDir = os.path.dirname(filepath)

        docOrgDB = DocOrganismDB.loadFromFile(fileDir + "/organism."+dbtype)

        seenRels = set()

        geneSymbolsNormalized = 0

        seenMirnas = {}

        ignoredDocs = set()
        takenDocs = set()

        with open(filepath, 'r') as fin:


            #miR-182 mir-182 MIRNA   CBFA3   RUNX3   GENE    29054094        True    True    [('21', 'V21', 'POS', 'express', '29054094.2.3', False, (92, 99), (82, 87), (62, 72), 'all_rels')]

            for lineIdx, line in enumerate(fin):

                line = line.strip()

                aline = line.split('\t')

                turnEvs= False


                if switchLR:
                    turnEvs = True

                    tmp = aline[0:3]

                    aline[0] = aline[3]
                    aline[1] = aline[4]
                    aline[2] = aline[5]

                    aline[3] = tmp[0]
                    aline[4] = tmp[1]
                    aline[5] = tmp[2]


                lid = aline[0]
                rid = aline[3]

                if ret.l_ont_based and lReplaceSc:
                    lid = lid.replace('_', ':', 1)

                if ret.r_ont_based and rReplaceSc:
                    rid = rid.replace('_', ':', 1)

                if ltype == rtype and lid > rid:
                    tmp = lid
                    lid = rid
                    rid = tmp

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
                        seenMirnas[origRid] = (org, lid)

                elif rtype == 'mirna':

                    if rid in seenMirnas:
                        org, rid = seenMirnas[rid]
                    else:
                        origRid = rid
                        org, rid = cls.harmonizeMIRNA(rid)
                        seenMirnas[origRid] = (org, rid)


                docOrgs = set()
                if org != None:
                    docOrgs.add(org)

                docid = aline[6]

                if ignoreDocIDs != None and docid in ignoreDocIDs:
                    ignoredDocs.add(docid)
                    continue

                takenDocs.add(docid)

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

                        if isinstance(rel[-4], int):
                            newrel.trustv = rel[-4:]

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

                for relIdx, rel in enumerate(relations):

                    if len(docOrgs) > 0:
                        rel.orgs = tuple(docOrgs)

                    #rel.data_id = dataID

                    ret.ltype2rel[lid].add(rel)
                    ret.rtype2rel[rid].add(rel)

                ret.all_ltypes.add(lid)
                ret.all_rtypes.add(rid)

        for x in ret.ltype2rel:
            for elem in ret.ltype2rel[x]:
                if elem.data_id == None:
                    print("Elem with no data id!")

        print("Gene Symbols Normalized", geneSymbolsNormalized)

        addTesters = False
        if addTesters:


            prepareTMAnalysis = False
            prepareVerbAnalysis = True

            if prepareTMAnalysis:
                allRelations = list()

                for lid in ret.ltype2rel:
                    allLidRels = ret.ltype2rel[lid]

                    allLidRels = [x for x in allLidRels if x.assocSent != None]

                    allRelations += allLidRels


                allSeen = []

                for tester in ['Tester1', 'Tester2', 'Tester3', 'all']:

                    relems = cls.choices(allRelations, 200, allSeen)

                    allSeen += relems

                    for relem in relems:
                        print(tester + "\t"+ relem.lid + "\t" + relem.rid + "\t" + relem.assocSent)




        print("Loaded file", filepath)
        print("Accepted Doc IDs", len(takenDocs))
        print("Rejected Doc IDs", len(ignoredDocs))

        return ret


if __name__ == '__main__':

    pmidBase = "/home/mjoppich/ownCloud/data/miRExplore/textmine/aggregated_pmid"
    MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.spacy.pmid", ltype="gene", rtype="mirna")
