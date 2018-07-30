from collections import defaultdict

from utils.DataFrame import DataFrame


class SymbolEnsemblDB:

    def __init__(self):

        self.org2convert = defaultdict(lambda: defaultdict(set))
        self.org2ens2symbol = defaultdict(lambda: dict())


    def get_all_genes(self, geneid):

        if geneid == None:
            return None


        geneid = geneid.upper()

        ret={}

        for org in self.org2convert:

            if geneid in self.org2convert[org]:

                ret[org] = self.org2convert[org][geneid]

        return ret

    def get_symbol_from_ens(self, org, ens):

        if not org in self.org2convert:
            return None


        if ens in self.org2ens2symbol[org]:
            return self.org2ens2symbol[org][ens]

        return None

    def get_symbol_for_ensembl(self, ensid):

        for org in self.org2ens2symbol:

            if ensid in self.org2ens2symbol[org]:
                return self.org2ens2symbol[org][ensid]

        return None



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
                    ret.org2ens2symbol['mmu'][geneid] = symbol

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
                ret.org2ens2symbol['hsa'][ensemblGeneID] = symbol


        return ret


if __name__ == '__main__':

    db = SymbolEnsemblDB.loadFromFile("/home/mjoppich/ownCloud/data/miRExplore/obodir/sym2ens")

    ret = db.get_symbol_for_ensembl('ENSMUSG00000024171')
    ret = db.get_all_genes("CXCR4")

    print(ret)