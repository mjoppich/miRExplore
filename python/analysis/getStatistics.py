import argparse, os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


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

    #{'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
    #{'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
    #{'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
    elemTerm = diseaseObo['DOID:1936']
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

            try:
                mirObj = miRNA(mirna)
                allMirnas.add(mirObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]))

            except:

                print(mirna)
                exit(-1)

    print("Number of mirnas with interaction", len(allMirnas))


    # number of gene-mirna interactions with disease association


    distinctInteractionsWithDisease = set()
    interactionsWithDisease = set()

    interactionsWithAthero = set()
    distinctInteractionsWithAthero = set()
    atheroPubmeds = set()
    cvPubmeds = set()
    allPubmeds = set()

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

                mirObj = miRNA(rel.rid)
                baseMirna = mirObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

                origTuple = (rel.lid, rel.rid)
                intTuple = (rel.lid, baseMirna)

                totalInteractions.add(intTuple)
                totalDistinctInteractions.add( origTuple )

                totalMirnas.add(intTuple[1])
                totalDistinctMirnas.add(origTuple[1])

                docID = rel.docid
                retVal = pmid2disease.getDOC(docID)

                allPubmeds.add(docID)

                if retVal != None:


                    distinctInteractionsWithDisease.add( origTuple )
                    interactionsWithDisease.add(intTuple)

                    for docDisease in retVal:

                        if docDisease['termid'] in elemTerms:
                            interactionsWithAthero.add(intTuple)
                            distinctInteractionsWithAthero.add(origTuple)
                            atheroPubmeds.add(docID)

                        if docDisease['termid'] in cvTerms:
                            interactionsWithCV.add(intTuple)
                            distinctInteractionsWithCV.add(origTuple)
                            cvPubmeds.add(docID)


    print("Different miRNAs", len(totalMirnas))
    print("Distinct miRNAs", len(totalDistinctMirnas))

    print("total Interactions", len(totalInteractions))
    print("total distinct Interaction", len(totalDistinctInteractions))
    print("total mirnas in interactions", len(set([x[1] for x in totalInteractions])))


    print("total interactions with disease", len(interactionsWithDisease))
    print("total distinct interactions with disease", len(distinctInteractionsWithDisease))
    print("total mirnas in interactions with disease", len(set([x[1] for x in interactionsWithDisease])))

    atheroMirnas = set([x[1] for x in interactionsWithAthero])
    atheroGenes = set([x[0] for x in interactionsWithAthero])

    print("total interactions with athero", len(interactionsWithAthero))
    print("total distinct interactions with athero", len(distinctInteractionsWithAthero))
    print("total mirnas in interactions with athero", len(atheroMirnas))
    print("total genes in interactions with athero", len(atheroGenes))
    print("total pubmeds for interactions with athero", len(atheroPubmeds))
    print("total pubmeds for interactions with cv", len(cvPubmeds))
    print("total pubmeds for interactions", len(allPubmeds))

    print("total interactions with cv", len(interactionsWithCV))
    print("total distinct interactions with cv", len(distinctInteractionsWithCV))
    print("total mirnas in interactions with cv", len(set([x[1] for x in interactionsWithCV])))

    # number of gene-mirna interactions with cv association

    print(allMirnas.difference(totalMirnas))

    print()
    print()
    print()

    print("Athero miRNAs")
    for x in atheroMirnas:
        print(x)

    print()
    print()
    print()

    print("Athero genes")
    for x in atheroGenes:
        print(x)