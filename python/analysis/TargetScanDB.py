import re
from collections import defaultdict

from openpyxl import load_workbook


class TargetScanDB :

    def __init__(self):

        self.elems = []
        self.gene2mirnas = defaultdict(list)

    def make_dictionary(self):
        for elem in self.elems:
            self.gene2mirnas[elem[0]].append(elem)


    @classmethod
    def from_tsv(cls, filelocation="/mnt/c/ownCloud/data/miRExplore/targetscan/targetscan_ws_85.tsv"):

        tsdb = TargetScanDB()

        with open(filelocation, 'r') as fin:



            for idx, row in enumerate(fin):

                if idx == 0:
                    continue

                arow = row.strip().split('\t')

                gene = arow[0].upper()
                mirna = arow[1]
                score = float(arow[2])
                percentile = int(arow[3])

                mirna = mirna.replace('mmu-', '').replace('hsa-', '')


                tsdb.elems.append((gene, mirna, score, percentile))


        return tsdb


if __name__ == '__main__':

    tsdb = TargetScanDB.from_tsv()

    for x in tsdb.elems:
        print(x)

