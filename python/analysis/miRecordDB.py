import re
from collections import defaultdict

from openpyxl import load_workbook

from synonymes.SynonymFile import Synfile
from textdb.AbstractDBClasses import DataBaseDescriptor, DataBaseRel


class MiRecordRel(DataBaseRel):

    def __init__(self, lent, rent, docid, relID):

        self.lent = lent
        self.rent = rent

        self.relID = relID
        self.pmid = docid

        self.data_source='mirecords'

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
        return self.pmid

    def toJSON(self):

        retJSON = {
            'docid': self.docid,
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,
            'data_source': self.data_source,
            'data_id': self.relID
        }

        if self.orgs != None:
            retJSON['orgs'] = tuple(self.orgs)

        return retJSON

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
    def loadFromFile(cls, filelocation=None, ltype='gene', rtype='mirna', normGeneSymbols=None, getDocs=False):

        ret = miRecordDB(ltype, rtype)
        print(filelocation)

        wb = load_workbook(filelocation)
        ws = wb.active

        geneSymbolsNormalized = 0

        docs=set()
        for row in ws:
            if row[0].row < 2:
                continue

            pubmed = str(row[0].value)
            gene = row[2].value

            if gene == None or not row[2].data_type=='s':
                continue

            gene = gene.strip().upper()

            if gene == 'UNKNOWN':
                continue

            (org, mirna) = cls.harmonizeMIRNA(row[6].value)

            if not org in ['hsa', 'mmu']:
                continue

            if gene in normGeneSymbols:
                gene = normGeneSymbols[gene]
                geneSymbolsNormalized += 1
            elif gene[-1].isdigit():
                agenename = gene[:-2] + " " + gene[-1]
                if agenename in normGeneSymbols:
                    gene = normGeneSymbols[agenename]
                    geneSymbolsNormalized += 1


            targetSitePosition = row[16].value

            if mirna == None:
                continue

            props = {
                'TARGET_SITE_POSITION': targetSitePosition
            }

            if pubmed != None and len(pubmed) > 0:
                docs.add(pubmed)

            relation = MiRecordRel((gene, ltype), (mirna, rtype), pubmed, 'MIRECORDS_'+str(row[0].row))
            if org != None:
                relation.orgs = tuple(set([org]))

            ret.ltype2rel[gene].add(relation)
            ret.rtype2rel[mirna].add(relation)

            ret.all_ltypes.add(gene)
            ret.all_rtypes.add(mirna)

        print("Gene Symbols Normalized", geneSymbolsNormalized)

        if getDocs:
            return ret, docs

        return ret


if __name__ == '__main__':

    mirecords = miRecordDB.loadFromFile()

    for x in mirecords.elems:
        print(x)

