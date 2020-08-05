import sys, os
from collections import Counter, OrderedDict, defaultdict
from upsetplot import from_contents,plot

import networkx
from natsort import natsorted

from synonymes.GeneOntology import GeneOntology
from synonymes.mirnaID import miRNAPART, miRNA, miRNACOMPARISONLEVEL

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher
from utils.tmutils import normalize_gene_names

import pandas as pd
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt

figidx = 0


cbn2restrict = OrderedDict([(

    "CV-IPN-Endothelial_cell_activation_1", {
        "sentences": "false",
        "cells": [
                {"group": "cells", "name": "endothelial cell",  "termid": "META:52"},
                {"group": "cells", "name": "HUVEC-C",  "termid": "META:02689"}
              ],
        "ncits": [
            {"group": "ncits", "name": "Vascular Cell Adhesion Protein 1", "termid": "NCIT:C48214"},
            {"group": "ncits", "name": "Intercellular Adhesion Molecule 1", "termid": "NCIT:C17304"}
        ]
        })
    ])



def acceptEvidence(ev):

    #return True

    if ev['data_source'] == 'miRTarBase':

        if "Weak" in ev['functional_type']:
            return False

    return True


gene2result = defaultdict(set)



cbn2mirna2evs = defaultdict(lambda: defaultdict(set))

rel2docs = defaultdict(set)

for geneName in ["MEG3", "NEAT1","MMP2", "MMP9", "MMP12"]:


    requestData = {}
    requestData['gene'] = [geneName]
    requestData['sentences'] = "false"

    #print(requestData)

    _,_,_, json = DataBasePlotter.fetchGenes(requestData, gene2name=None, minPMIDEvCount=0, minTgtCount=0, acceptEv=acceptEvidence, MIRNASTRPARTS=[miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

    for x in json["rels"]:

        mirnastr = x["rid"]

        omir = miRNA(mirnastr)
        simpleStr = omir.getStringFromParts(miRNACOMPARISONLEVEL.PRECURSOR.value)

        if simpleStr == "miR-7":
            simpleStr = "let-7"

        #print(x["lid"], x["rid"], simpleStr)
        gene2result[geneName].add(simpleStr)

        for ev in x["evidences"]:
            if "docid" in ev:
                rel2docs[(geneName, simpleStr)].add((ev["docid"], ev["data_source"]))
            else:
                rel2docs[(geneName, simpleStr)].add(ev["data_source"])


for x in gene2result:
    print(x, len(gene2result[x]))

#upIn = from_contents(gene2result)
#plot(upIn, subset_size="auto")
#plt.show()


from plots.DotSetPlot import DotSetPlot

gene2list = defaultdict(lambda: defaultdict(set))
for x in gene2result:

    for mirna in gene2result[x]:
        gene2list[mirna]["Target Genes"].add(x)


filteredList = defaultdict(lambda: defaultdict(set))

for mirna in gene2list:

    if len(gene2list[mirna]["Target Genes"]) > 1:
        filteredList[mirna] = gene2list[mirna]

    if all([x in gene2list[mirna]["Target Genes"] for x in ["MMP2", "MMP9", "MMP12"]]):

        tgenes = gene2list[mirna]["Target Genes"]

        print(mirna, sorted(tgenes))

        for tgene in tgenes:
            print(mirna, tgene, rel2docs[(tgene, mirna)])


DotSetPlot().plot({"Target Genes": [x for x in gene2result]}, filteredList)#, max=30)

plt.show()