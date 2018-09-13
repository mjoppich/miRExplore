from synonymes.SynonymFile import Synfile


def normalize_gene_names(path="/home/mjoppich/ownCloud/data/miRExplore/textmine/synonyms/hgnc_no_withdrawn.syn"):
    geneNameSynFile = Synfile(path)
    normalizeGeneNames = {}

    for sid in geneNameSynFile.mSyns:

        synonym = geneNameSynFile.mSyns[sid]

        for syn in synonym.syns:

            psyn = syn.upper()

            if psyn == sid:
                continue

            normalizeGeneNames[psyn] = sid

            if psyn == 'LIN28':
                print(psyn, sid)

    return normalizeGeneNames
