from collections import defaultdict

from synonymes.mirnaID import miRNA, miRNAPART
from utils.tmutils import normalize_gene_names


def harmonizeMIRNA( mirna):
    """

    :param mirna:
    :return: tries to return a normalized name ...

    9761 microRNA
   7958 MicroRNA
   2311 MiRNA
   2191 miRNA
   1844 hsa
   1440 let
    578 miRNAS
    437 MICRORNA
    299 microRNAS
    256 MIRNA
    155 mmu
    125 Micro
    116 micro
    """

    possibleOrgStarts = ['mmu', 'hsa']

    recOrg = None

    for x in possibleOrgStarts:

        if mirna.startswith(x + "-"):
            recOrg = x
            mirna = mirna.replace(x + "-", "", 1)
            break

    possibleStarts = ['microRNA', 'MicroRNA', 'MiRNA', 'miRNA', 'MICRORNA', 'microRNA', 'MIRNA']

    for x in possibleStarts:

        if mirna.startswith(x + "-"):
            mirna = mirna.replace(x, "miR", 1)
            break

    try:
        oMirna = miRNA(mirna)

        mirna = oMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

    except:
        pass

    return (recOrg, mirna)


if __name__ == '__main__':

    base = "/mnt/c/ownCloud/data/mirpredict/"

    normGeneSymbols = normalize_gene_names(path="/mnt/c/ownCloud/data/miRExplore/obodir/hgnc_no_withdrawn.syn")

    mirna2gene = defaultdict(set)
    with open(base + "transmir.hsa.tsv", 'r') as fin:

        for line in fin:

            line = line.strip().split("\t")

            gene = line[0]
            mirna = line[1]

            gene = gene.upper()
            if gene in normGeneSymbols:
                gene = normGeneSymbols[gene]

            org, mirna = harmonizeMIRNA(mirna)

            rtype = line[4]

            mirna2gene[mirna].add( (gene, rtype) )


    with open("/mnt/c/ownCloud/data/mirpredict/transmir_rels.tsv", 'w') as outfile:
        for mirna in mirna2gene:

            allrels = mirna2gene[mirna]
            for gene, rtype in allrels:

                print(mirna, gene, rtype, "transmir.hsa", sep="\t", file=outfile)