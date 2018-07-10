from collections import defaultdict

from utils.DataFrame import DataFrame


class NcitTermSymbolDB:


    def __init__(self):

        self.org_term2symbol = defaultdict(lambda: defaultdict(set))


    @classmethod
    def loadFromFolder(cls):

        baseFolder = "/mnt/c/ownCloud/data/miRExplore/obodir/map_ncit_syms/"


        ncit2swissprotDF = DataFrame.parseFromFile(baseFolder + "NCIt-SwissProt_Mapping.txt")
        hsa_mmu_orthologuesDF = DataFrame.parseFromFile(baseFolder + "human_mouse_orthologues_ensembl.tsv")
        ensembl2hgncDF = DataFrame.parseFromFile(baseFolder + "ensembl_hgnc_uniprot.txt")
        ensembl2mgiDF = DataFrame.parseFromFile(baseFolder + "ensembl_mgi_hgnc.tsv")


        print("ncit2swissprot")
        print(ncit2swissprotDF.getHeader())

        print("hsa_mmu_orthologues")
        print(hsa_mmu_orthologuesDF.getHeader())

        print("ensembl2hgnc")
        print(ensembl2hgncDF.getHeader())

        print("ensembl2mgi")
        print(ensembl2mgiDF.getHeader())


        ncit2swissprot = {}

        for row in ncit2swissprotDF:
            ncit2swissprot[row['NCIt Code']] = row['SwissProt ID']


        swissprot2ensembl = defaultdict(set)
        for row in ensembl2hgncDF:

            ensemblGeneID = row['Gene stable ID']
            swissprotID = row['UniProtKB/Swiss-Prot ID']

            swissprot2ensembl[swissprotID].add(ensemblGeneID)


        hsaEnsembl2Sym = {}
        for row in ensembl2hgncDF:
            ensemblGeneID = row['Gene stable ID']
            hgncSymbol = row['HGNC symbol']

            if hgncSymbol == None or len(hgncSymbol) == 0:
                continue

            hsaEnsembl2Sym[ensemblGeneID] = hgncSymbol

        mmuEnsembl2Sym = {}
        for row in ensembl2mgiDF:
            ensemblGeneID = row['Gene stable ID']
            mgiSymbol = row['MGI symbol']

            if mgiSymbol == None or len(mgiSymbol) == 0:
                continue

            mmuEnsembl2Sym[ensemblGeneID] = mgiSymbol

        hsaEnsembl2mmuEnsembl = defaultdict(set)
        for row in hsa_mmu_orthologuesDF:

            hsaID = row['Gene stable ID']
            mmuID = row['Mouse gene stable ID']

            hsaEnsembl2mmuEnsembl[hsaID].add(mmuID)

        retDB = NcitTermSymbolDB()

        unknownSwissprots = set()

        for ncitID in ncit2swissprot:

            swissprotID = ncit2swissprot[ncitID]

            hsaEnsembls = swissprot2ensembl.get(swissprotID, None)

            if hsaEnsembls == None:
                unknownSwissprots.add(swissprotID)
                print("no ensembl for swissprot", swissprotID)
                continue

            for hsaEnsembl in hsaEnsembls:
                hsaSym = hsaEnsembl2Sym.get(hsaEnsembl, None)

                if hsaSym == None:
                    print("no hsa symols for", hsaEnsembl)
                    continue

                retDB.org_term2symbol['hsa'][ncitID].add(hsaSym)

                mmuIDs = hsaEnsembl2mmuEnsembl.get(hsaEnsembl, None)

                if mmuIDs == None:
                    print("no mmu ortholog for", hsaEnsembl)
                    continue

                for mmuID in mmuIDs:
                    mmuSym = mmuEnsembl2Sym.get(mmuID)

                    if mmuSym == None:
                        print("no mmu symbol for", mmuID)

                    retDB.org_term2symbol['mmu'][ncitID].add(mmuSym)

        print(unknownSwissprots)

        return

        for elem in swissprot2ensembl:
            if len(swissprot2ensembl[elem]) > 1:
                print(elem, swissprot2ensembl[elem])




if __name__ == '__main__':

    NcitTermSymbolDB.loadFromFolder()