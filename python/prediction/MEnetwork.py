from collections import defaultdict

from synonymes.mirnaID import miRNA, miRNAPART
from textdb.MiGenRelDB import MiGenRelDB
from utils.tmutils import normalize_gene_names

if __name__ == '__main__':

    pmidBase = "/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmid/"

    normGeneSymbols = normalize_gene_names(path="/mnt/c/ownCloud/data/miRExplore/obodir/hgnc_no_withdrawn.syn")


    mirelPMIDhsa = MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene",
                                           normGeneSymbols=normGeneSymbols, switchLR=True)
    mirelPMIDmmu = MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene",
                                           normGeneSymbols=normGeneSymbols, switchLR=True)


    # ret.ltype2rel[lid].add(rel)

    gene2mirna = defaultdict(set)

    for mirelPMID in [mirelPMIDhsa, mirelPMIDmmu]:

        for gene in mirelPMID.ltype2rel:

            allMirRels = mirelPMID.ltype2rel[gene]

            allNormedMirs = set()

            for rel in allMirRels:

                relMir = rel.rid

                try:
                    testMirna = miRNA(relMir)
                    nMirna = testMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])
                    allNormedMirs.add(nMirna)
                except:
                    print("error loading miRNA", relMir)

            for x in allNormedMirs:
                gene2mirna[gene].add(x)



    with open("/mnt/c/ownCloud/data/mirpredict/mirexplore_rels.tsv", 'w') as outfile:
        for gene in gene2mirna:

            allmirs = gene2mirna[gene]
            for mirna in allmirs:

                print(gene, mirna, "repress", "mirexplore", sep="\t", file=outfile)