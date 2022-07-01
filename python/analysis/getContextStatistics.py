import argparse
import datetime
import os, sys
import pickle
from collections import defaultdict
from itertools import chain, combinations

from synonymes.GeneOntology import GeneOntology
from textdb.DIANATarbaseDB import DIANATarbaseDB
from textdb.MiGenRelDB import MiGenRelDB
from textdb.MirTarBaseDB import MirTarBaseDB
from textdb.PMID2XDB import PMID2XDB
from utils.tmutils import normalize_gene_names
from analysis.miRecordDB import miRecordDB

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


def easyPMIDFinder(filepath):

    foundPMIDs = set()

    with open(filepath, 'r') as fin:

        for line in fin:
            line = line.strip()

            aline = line.split('\t')

            pmid = aline[0]

            if pmid != None and len(pmid) > 0:
                foundPMIDs.add(pmid)

    return foundPMIDs


def easyPMIDFinderCells(filepath):
    foundPMIDs = set()

    with open(filepath, 'r') as fin:

        for line in fin:
            line = line.strip()

            aline = line.split('\t')

            pmid = aline[0]
            cellID = aline[1]

            if pmid != None and len(pmid) > 0:
                if cellID != None and len(cellID) > 3 and cellID.startswith("CL:"):
                    foundPMIDs.add(pmid)

    return foundPMIDs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-p', '--pmidBase', type=str,
                        help='Base for Textmining. Includes aggregated_ and results folder',
                        required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)

    args = parser.parse_args()


    if False:
        dbs2pmids = defaultdict(set)

        testStop = -1


        print("Loading hsa")
        mirelPMIDhsa, takenDocsHSA = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene",
                                               normGeneSymbols=normGeneSymbols, switchLR=True, stopAfter=testStop, getDocs=True)

        mirelPMIDhsa=None
        print("Loading mmu")
        mirelPMIDmmu, takenDocsMMU = MiGenRelDB.loadFromFile(args.pmidBase + "/mirna_gene.mmu.pmid", ltype="mirna", rtype="gene",
                                               normGeneSymbols=normGeneSymbols, switchLR=True, stopAfter=testStop, getDocs=True)

        mirelPMIDmmu = None

        print(datetime.datetime.now(), "Loading miRecords")
        recordsDB, mirecordsDocs = miRecordDB.loadFromFile(filelocation=args.obodir + "/mirecords_v4.xlsx", normGeneSymbols=normGeneSymbols, getDocs=True)
        recordsDB = None
        print(datetime.datetime.now(), "Loading miRTarBase")
        mirtarbaseDB, mirtarbaseDocs = MirTarBaseDB.loadFromFile(filepath=args.obodir + "/miRTarBase.csv", normGeneSymbols=normGeneSymbols, getDocs=True)
        mirtarbaseDB = None
        print(datetime.datetime.now(), "Loading hsa_mmu.diana")
        #dianaDB, celllInfos = DIANATarbaseDB.loadFromFile(args.obodir + "/hsa_mmu.diana.csv", normGeneSymbols=normGeneSymbols)

        dbs2pmids["TM"] = dbs2pmids["TM"].union(takenDocsHSA)
        dbs2pmids["TM"] = dbs2pmids["TM"].union(takenDocsMMU)
        dbs2pmids["TM"] = dbs2pmids["TM"].union(mirecordsDocs)
        dbs2pmids["TM"] = dbs2pmids["TM"].union(mirtarbaseDocs)

        print(datetime.datetime.now(), "Loading ontologies")

        print(datetime.datetime.now(), "Loading ontologies finished")
        print(datetime.datetime.now(), "Loading GO")
        goPMIDs = easyPMIDFinder(args.pmidBase + "/go.pmid")
        dbs2pmids["GO"] = goPMIDs

        print(datetime.datetime.now(), "Loading Disease")
        diseasePMIDs = easyPMIDFinder(args.pmidBase + "/disease.pmid")
        dbs2pmids["DOID"] = diseasePMIDs

        pmid2disease = None

        pmid2fma=None
        #print(datetime.datetime.now(), "Loading FMA")
        #pmid2fma = PMID2XDB.loadFromFile(args.pmidBase + "/model_anatomy.pmid", fmaObo)
        print(datetime.datetime.now(), "Loading cellline")
        cellPMIDs = easyPMIDFinderCells(args.pmidBase + "/celllines.pmid")
        dbs2pmids["CELLS"] = cellPMIDs


        print(datetime.datetime.now(), "Loading ncit")
        ncitPMIDs = easyPMIDFinder(args.pmidBase + "/ncit.pmid")
        dbs2pmids["NCIT"] = ncitPMIDs

        with open("/mnt/d/pmidsindims.pickle", 'wb') as fout:
            pickle.dump(dbs2pmids, fout)
    else:

        with open("/mnt/d/pmidsindims.pickle", 'rb') as fout:
            dbs2pmids = pickle.load(fout)


    outdf = DataFrame()
    outdf.addColumns(["Subset", "Number of PMIDs"])

    allDims = [x for x in dbs2pmids]

    allPowerSets = powerset(sorted(allDims))

    allPMIDs = set()
    for x in dbs2pmids:
        allPMIDs = allPMIDs.union(dbs2pmids[x])


    for pset in allPowerSets:

        if len(pset) == 0:
            continue


        setPMIDs = set([x for x in allPMIDs])
        for dim in pset:
            setPMIDs = setPMIDs.intersection( dbs2pmids[dim] )

        print(pset, len(setPMIDs))

        drdict = {
            "Subset": ", ".join(pset),
            "Number of PMIDs": len(setPMIDs)
        }

        dr = DataRow.fromDict(drdict)
        outdf.addRow(dr)

    print(outdf._makeLatex())







