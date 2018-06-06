import os, sys
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

        self.data_source=datasource
        self.data_id = dataID


        if intuple != None:
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

        return {
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

            'data_source': self.data_source,
            'data_id': self.data_id
        }

class MiGenRelDB(DataBaseDescriptor):


    def __init__(self, ltype, rtype):

        super().__init__()

        self.ltyped = ltype
        self.rtyped = rtype

    @property
    def ltype(self):
        return self.ltyped

    @property
    def rtype(self):
        return self.rtyped

    @classmethod
    def choices(cls, elems, k=1, excluded = []):

        resvec = []

        while len(resvec) < k:

            relem = random.choice(elems)

            if not relem in resvec and not relem in excluded:
                resvec.append(relem)

        return resvec


    @classmethod
    def loadFromFile(cls, filepath, ltype='gene', rtype='mirna'):


        ret = MiGenRelDB(ltype, rtype)

        file_base = os.path.basename(filepath)

        rel2did = {}

        with open(filepath, 'r') as fin:

            #PCBP1	miR-3978	ORGMIR4299	29138420	True	True	[('GM', 'GVM', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (6, 16)), ('GM', 'GMV', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (35, 45))]

            for lineIdx, line in enumerate(fin):

                line = line.strip()

                aline = line.split('\t')

                lid = aline[0]
                rid = aline[1]
                (org, rid) = cls.harmonizeMIRNA(rid)

                docid = aline[3]

                sameParagraph = aline[4] == 'True'
                sameSentence = aline[5] == 'True'

                if not sameSentence:
                    continue

                tmRelations = None if aline[6] == 'None' else eval(aline[6])

                if tmRelations != None:
                    allrels = set()

                    for relIdx, rel in enumerate(tmRelations):
                        newrel = MiRGeneRel(rel, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', None)

                        dataID = file_base + "_" + str(lineIdx) + "_" + str(relIdx)
                        rel2did[newrel] = dataID

                        allrels.add(newrel)

                    relations = allrels
                else:
                    relations = set([MiRGeneRel(None, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', None)])

                for relIdx, rel in enumerate(relations):


                    #rel.data_id = dataID

                    ret.ltype2rel[lid].add(rel)
                    ret.rtype2rel[rid].add(rel)

                ret.all_ltypes.add(lid)
                ret.all_rtypes.add(rid)

        for x in ret.ltype2rel:
            for elem in ret.ltype2rel[x]:
                if elem in rel2did:
                    elem.data_id = rel2did[elem]

        for x in ret.rtype2rel:
            for elem in ret.rtype2rel[x]:
                if elem in rel2did:
                    elem.data_id = rel2did[elem]

        if __name__ == '__main__':
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


        return ret


if __name__ == '__main__':

    pmidBase = "/home/mjoppich/ownCloud/data/miRExplore/textmine/aggregated_pmid"
    MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.spacy.pmid", ltype="gene", rtype="mirna")
