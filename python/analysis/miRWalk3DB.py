import re
from collections import defaultdict

from openpyxl import load_workbook

import pickle


class miRWalk3DB :

    def __init__(self, pickleFile='/tmp/mirwalk.db'):

        self.elems = []

        if pickleFile != None:

            try:
                with open(pickleFile, 'rb') as fin:
                    self.elems = pickle.load(fin)

            except:
                pass




    #hsa_miRWalk_3UTR.short.tsv  hsa_miRWalk_5UTR.short.tsv  hsa_miRWalk_CDS.short.tsv

    @classmethod
    def from_xslx(cls, folderPath="/mnt/c/ownCloud/data/miRExplore/mirwalk/", files={'hUTR3': "hsa_miRWalk_3UTR.short.tsv",
                                                                                     'hUTR5': 'hsa_miRWalk_5UTR.short.tsv',
                                                                                     'hCDS':'hsa_miRWalk_CDS.short.tsv',
                                                                                     'mUTR3': "mmu_miRWalk_3UTR.short.tsv",
                                                                                     'mUTR5': 'mmu_miRWalk_5UTR.short.tsv',
                                                                                     'mCDS': 'mmu_miRWalk_CDS.short.tsv',
                                                                                     }, targetGenes=None, minScore=1.0):


        targetGenes = [x.upper() for x in targetGenes] if targetGenes != None else None
        mirwalk = miRWalk3DB(pickleFile=None)

        for fileType in files:

            filelocation = folderPath + files[fileType]

            with open(filelocation, 'r') as fin:

                for idx, line in enumerate(fin):

                    if idx == 0:
                        continue

                    aline = line.split('\t')

                    mirna = aline[0]
                    gene = aline[1]
                    score = float(aline[2])

                    if not score >= minScore:
                        continue

                    if targetGenes != None and not gene.upper() in targetGenes:
                        continue

                    mirwalk.elems.append((gene, mirna, fileType))

        with open('/tmp/mirwalk.db', 'wb') as fout:

            pickle.dump(mirwalk.elems, fout)


        return mirwalk




if __name__ == '__main__':

    mirecords = miRWalk3DB.from_xslx()


