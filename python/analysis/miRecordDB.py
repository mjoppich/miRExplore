import re
from collections import defaultdict

from openpyxl import load_workbook


class MiRecordRel:

    def __init__(self, lent, rent, docid, relID):

        self.ltype = lent[1]
        self.rtype = rent[1]

        self.lid = lent[0]
        self.rid = rent[0]

        self.relID = relID
        self.docid = docid

        self.data_source='mirecords'


    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):

        return hash(tuple([(x, self.__dict__[x]) for x in sorted(self.__dict__)]))


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

class miRecordDB :


    def __init__(self, ltype, rtype):

        self.ltype = ltype
        self.rtype = rtype

        self.ltype2rel = defaultdict(set)
        self.rtype2rel = defaultdict(set)

        self.all_ltypes = set()
        self.all_rtypes = set()

    def get_evidence_docids(self):

        docIDs = set()

        for lid in self.ltype2rel:
            for ev in self.ltype2rel[lid]:
                docIDs.add(ev.docid)

        return docIDs

    def get_rels(self, etype, eid):

        if etype == self.ltype:
            return self.get_lid_rels(eid)
        elif etype == self.rtype:
            return self.get_rid_rels(eid)

        return None

    def get_lid_rels(self, geneID):
        return self.ltype2rel.get(geneID, None)

    def get_rid_rels(self, mirnaID):
        return self.rtype2rel.get(mirnaID, None)



    @classmethod
    def from_xslx(cls, filelocation="/mnt/c/ownCloud/data/miRExplore/miRecords/mirecords_v4.xlsx", ltype='gene', rtype='mirna'):

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
            mirna = row[6].value

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

    mirecords = miRecordDB.from_xslx()

    for x in mirecords.elems:
        print(x)

