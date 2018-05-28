import re
from collections import defaultdict

from openpyxl import load_workbook

from textdb.AbstractDBClasses import DataBaseDescriptor, DataBaseRel


class MiRecordRel(DataBaseRel):

    def __init__(self, lent, rent, docid, relID):

        self.lent = lent
        self.rent = rent

        self.relID = relID
        self.pmid = docid

        self.data_source='mirecords'

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
        return self.pmid

    def toJSON(self):

        return {
            'docid': self.docid,
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,
            'data_source': self.data_source,
            'data_id': self.relID
        }

class miRecordDB(DataBaseDescriptor):


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
    def loadFromFile(cls, filelocation=None, ltype='gene', rtype='mirna'):

        ret = miRecordDB(ltype, rtype)

        wb = load_workbook(filelocation)
        ws = wb.active

        for row in ws:
            if row[0].row < 2:
                continue

            pubmed = str(row[0].value)
            gene = row[2].value

            if gene == None or not row[2].data_type=='s':
                continue

            gene = gene.upper()

            (org, mirna) = cls.harmonizeMIRNA(row[6].value)

            targetSitePosition = row[16].value

            if mirna == None:
                continue

            props = {
                'TARGET_SITE_POSITION': targetSitePosition
            }

            relations = set([MiRecordRel((gene, ltype), (mirna, rtype), pubmed, 'MIRECORDS_'+str(row[0].row))])


            for rel in relations:
                ret.ltype2rel[gene].add(rel)
                ret.rtype2rel[mirna].add(rel)

            ret.all_ltypes.add(gene)
            ret.all_rtypes.add(mirna)


        return ret


if __name__ == '__main__':

    mirecords = miRecordDB.loadFromFile()

    for x in mirecords.elems:
        print(x)

