from collections import defaultdict

from utils.DataFrame import DataFrame


class SymbolEnsemblDB:

    def __init__(self):

        self.org2convert = defaultdict(lambda: defaultdict(set))


    def get_all_genes(self, geneid):

        if geneid == None:
            return None


        geneid = geneid.upper()

        ret = {"query": geneid}
        for org in self.org2convert:

            if geneid in self.org2convert[org]:

                ret[org] = self.org2convert[org][geneid]

        return ret



    @classmethod
    def loadFromFile(cls, basePath):



        ret = SymbolEnsemblDB()


        hsa=True
        mmu=True

        if mmu:
            MGIdata = DataFrame.parseFromFile(basePath + "/MRK_Sequence.rpt")

            transript2gene = {}

            with open(basePath + "/mmu_gene_transcript.txt") as fin:

                for line in fin:
                    aline = line.strip().split("\t")

                    egene = aline[0]
                    etran = aline[1]

                    transript2gene[etran] = egene

            for row in MGIdata:

                symbol = row['Marker Symbol'].upper()
                etranscripts = row['Ensembl transcript IDs'].split("|")

                if len(etranscripts) == 0:
                    continue

                egenes = set()

                for etran in etranscripts:
                    egeneid = transript2gene.get(etran, None)

                    if egeneid != None:
                        egenes.add(egeneid)

                for geneid in egenes:
                    ret.org2convert['mmu'][symbol].add(geneid)

        if hsa:

            hgncData = DataFrame.parseFromFile(basePath + "/hgnc_ext.tsv")

            for row in hgncData:

                symbol = row['Approved Symbol']

                if symbol == None:
                    continue

                symbol = symbol.upper()

                ensemblGeneID = row['Ensembl ID(supplied by Ensembl)']

                if ensemblGeneID == None:
                    continue

                ret.org2convert['hsa'][symbol].add(ensemblGeneID)


        return ret


if __name__ == '__main__':

    db = SymbolEnsemblDB.loadFromFile("/home/mjoppich/ownCloud/data/miRExplore/obodir/sym2ens")

    ret = db.get_all_genes('CXCR1')

    print(ret)