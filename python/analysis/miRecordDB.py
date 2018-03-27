import re
from collections import defaultdict

from openpyxl import load_workbook


class miRecordDB :

    def __init__(self):

        self.elems = []



    @classmethod
    def from_xslx(cls, filelocation="/mnt/c/ownCloud/data/miRExplore/miRecords/mirecords_v4.xlsx"):

        mirecords = miRecordDB()

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

            mirecords.elems.append((gene, mirna, pubmed, props))


        return mirecords


if __name__ == '__main__':

    mirecords = miRecordDB.from_xslx()

    for x in mirecords.elems:
        print(x)

