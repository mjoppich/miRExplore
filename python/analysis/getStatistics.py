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

    print("Loading hsa")
    mirelPMIDhsa = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene", normGeneSymbols=normGeneSymbols, switchLR=True, stopAfter=-1)
    print("Loading mmu")
    mirelPMIDmmu = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene", normGeneSymbols=normGeneSymbols, switchLR=True, stopAfter=-1)

    relDBs = [mirelPMIDhsa, mirelPMIDmmu]

    requiredPMIDs = set()
    for rdb in relDBs:

        assert (isinstance(rdb, DataBaseDescriptor))

        for rpmid in rdb.get_evidence_docids():
            requiredPMIDs.add(rpmid)

    diseaseObo = GeneOntology(args.obodir + "/doid.obo")


    #{'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
    #{'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
    elemTerm = diseaseObo['DOID:2349']
    elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]

    cvTerm = diseaseObo['DOID:1287']
    cvTerms = [x.term.id for x in cvTerm.getAllChildren()] + [cvTerm.id] + elemTerms

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


    distinctInteractionsWithDisease = set()
    interactionsWithDisease = set()

    interactionsWithAthero = set()
    distinctInteractionsWithAthero = set()

    interactionsWithCV = set()
    distinctInteractionsWithCV = set()

    totalDistinctMirnas = set()
    totalMirnas = set()

    totalInteractions = set()
    totalDistinctInteractions = set()

    for rdb in relDBs:

        for gene in rdb.ltype2rel:

            for rel in rdb.ltype2rel[gene]:

                relOrgs = rel.orgs
                if relOrgs == None:
                    relOrgs = set()
                if not ('mmu' in relOrgs or 'hsa' in relOrgs):
                    continue

                baseMirna = "miR-" + rel.rid.split("-")[1]

                origTuple = (rel.lid, rel.rid)
                intTuple = (rel.lid, baseMirna)

                totalInteractions.add(intTuple)
                totalDistinctInteractions.add( origTuple )

                totalMirnas.add(intTuple[1])
                totalDistinctMirnas.add(origTuple[1])

                docID = rel.docid
                retVal = pmid2disease.getDOC(docID)

                if retVal != None:


                    distinctInteractionsWithDisease.add( origTuple )
                    interactionsWithDisease.add(intTuple)

                    for docDisease in retVal:

                        if docDisease['termid'] in elemTerms:
                            interactionsWithAthero.add(intTuple)
                            distinctInteractionsWithAthero.add(origTuple)

                        if docDisease['termid'] in cvTerms:
                            interactionsWithCV.add(intTuple)
                            distinctInteractionsWithCV.add(origTuple)


    print("Different miRNAs", len(totalMirnas))
    print("Distinct miRNAs", len(totalDistinctMirnas))

    print("total Interactions", len(totalInteractions))
    print("total distinct Interaction", len(totalDistinctInteractions))

    print("total interactions with disease", len(interactionsWithDisease))
    print("total distinct interactions with disease", len(distinctInteractionsWithDisease))

    print("total interactions with athero", len(interactionsWithAthero))
    print("total distinct interactions with athero", len(distinctInteractionsWithAthero))

    print("total interactions with cv", len(interactionsWithCV))
    print("total distinct interactions with cv", len(distinctInteractionsWithCV))

    # number of gene-mirna interactions with cv association

    print(allMirnas.difference(totalMirnas))
