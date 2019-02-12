import argparse

from synonymes.GeneOntology import GeneOntology
from synonymes.mirnaID import miRNA, miRNAPART
from textdb.MiGenRelDB import MiGenRelDB
from textdb.PMID2XDB import PMID2XDB
from textdb.AbstractDBClasses import DataBaseDescriptor
from utils.tmutils import normalize_gene_names

normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-p', '--pmidBase', type=str, help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)

    args = parser.parse_args()

    mirelPMIDhsa = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene", normGeneSymbols=normGeneSymbols, switchLR=True)
    mirelPMIDmmu = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene", normGeneSymbols=normGeneSymbols, switchLR=True)

    relDBs = [mirelPMIDhsa, mirelPMIDmmu]

    requiredPMIDs = set()
    for rdb in relDBs:

        assert (isinstance(rdb, DataBaseDescriptor))

        for rpmid in rdb.get_evidence_docids():
            requiredPMIDs.add(rpmid)

    diseaseObo = GeneOntology(args.obodir + "/doid.obo")
    pmid2disease = PMID2XDB.loadFromFile(args.pmidBase + "/disease.pmid", diseaseObo, requiredPMIDs)


    # number of genes with interaction
    allGenes = set()

    for rdb in relDBs:
        allGenes = allGenes.union(set(rdb.all_ltypes))

    print("Number of genes with interaction", len(allGenes))

    # number of miRNAs with interaction
    ## restrict to miR-x
    allMirnas = set()
    for rdb in relDBs:
        for mirna in rdb.all_rtypes:

            mirObj = miRNA(mirna)
            allMirnas.add(mirObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID]))

    print("Number of mirnas with interaction", len(allMirnas))


    # number of gene-mirna interactions with disease association
    totalInteractions = 0
    totalInteractionsWithDisease = 0
    distinctInteractions = set()
    distinctInteractionsWithDisease = set()

    for rdb in relDBs:

        for gene in rdb.ltype2rel:

            for rel in rdb.ltype2rel[gene]:

                totalInteractions += 1

                baseMirna = "miR-" + rel.rid.split("-")[1]
                intTuple = (rel.lid, baseMirna)

                distinctInteractions.add( intTuple )

                docID = rel.docid
                retVal = pmid2disease.getDOC(docID)

                if retVal != None:

                    totalInteractionsWithDisease += 1
                    distinctInteractionsWithDisease.add( intTuple )

    print("Total miRNA-Gene interactions", totalInteractions)
    print("Total miRNA-Gene interactions with disease", totalInteractionsWithDisease)

    print("Distinct miRNA-Gene interactions", len(distinctInteractions))
    print("Distinct miRNA-Gene interactions", len(distinctInteractionsWithDisease))


    # number of gene-mirna interactions with cv association

