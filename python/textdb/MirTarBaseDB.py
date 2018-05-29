from collections import defaultdict

import os
from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor
from utils.DataFrame import DataFrame


class MirTarBaseRel(DataBaseRel):

    def __init__(self, lent, rent, datasource, dataID, expSupport, funcType, pubmedRef, organism):

        self.lent = lent
        self.rent = rent

        self.data_source=datasource
        self.data_id = dataID

        self.expSupport = tuple(expSupport)
        self.funcType = funcType
        self.pubmedRef = pubmedRef

        self.organism = organism

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
        return self.pubmedRef

    def toJSON(self):

        retJSON= {
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,

            'exp_support': self.expSupport,
            'functional_type': self.funcType,

            'data_source': self.data_source,
            'data_id': self.data_id
        }

        if self.pubmedRef != None:
            retJSON['docid']= self.pubmedRef

        if self.organism != None:
            retJSON['organism']= self.organism

        return retJSON

class MirTarBaseDB(DataBaseDescriptor):


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


        ret = MirTarBaseDB(ltype, rtype)
        file_base = os.path.basename(filepath)

        mirtarbaseEvidences = DataFrame.parseFromFile(filepath,
                                                      bConvertTextToNumber=False)

        for mirtEntry in mirtarbaseEvidences:

            lid = mirtEntry['Target Gene']
            rid = mirtEntry['miRNA']
            org, rid = cls.harmonizeMIRNA(rid)

            if not org in ['hsa', 'mmu']:
                continue

            organism = mirtEntry['Species (miRNA)']
            dataID = mirtEntry['miRTarBase ID']
            dataSource = 'miRTarBase'

            #Experiments     Support Type    References (PMID)

            docID = mirtEntry['References (PMID)']
            expSupport = mirtEntry['Experiments'].split("//")
            supType = mirtEntry['Support Type']

            relations = set([
                MirTarBaseRel((lid, ltype), (rid, rtype), dataSource, dataID, expSupport, supType, docID, organism)
                ])

            for rel in relations:

                ret.ltype2rel[lid].add(rel)
                ret.rtype2rel[rid].add(rel)

            ret.all_ltypes.add(lid)
            ret.all_rtypes.add(rid)

        return ret









