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
    def loadFromFile(cls, filepath, ltype='gene', rtype='mirna'):


        ret = MiGenRelDB(ltype, rtype)

        file_base = os.path.basename(filepath)

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
                        newrel = MiRGeneRel(rel, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', file_base+"_"+str(lineIdx)+"_"+str(relIdx))

                        allrels.add(newrel)

                    relations = allrels
                else:
                    relations = set([MiRGeneRel(None, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', file_base+"_"+str(lineIdx)+"_0")])

                for rel in relations:

                    ret.ltype2rel[lid].add(rel)
                    ret.rtype2rel[rid].add(rel)

                ret.all_ltypes.add(lid)
                ret.all_rtypes.add(rid)

        return ret

