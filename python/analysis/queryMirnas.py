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

downMirnas = ["hsa-miR-140",
"hsa-miR-1244",
"hsa-miR-6071",
"hsa-miR-4754",
"hsa-miR-4697",
"hsa-miR-4803",
"hsa-miR-3688",
"hsa-miR-6784",
"hsa-miR-143",
"hsa-miR-624",
"hsa-miR-24",
"hsa-miR-6809",
"hsa-miR-24",
"hsa-miR-605",
"hsa-miR-23b",
"hsa-miR-6853",
"hsa-miR-4504",
"hsa-miR-145",
"hsa-miR-3074",
"hsa-miR-27b",
"hsa-miR-4757",
"hsa-miR-198",
"hsa-miR-6852"]

upMirnas = [

"hsa-miR-4709",
"hsa-miR-3916",
"hsa-miR-3611",
"hsa-miR-5003",
"hsa-miR-7847",
"hsa-miR-6842",
"hsa-miR-4690",
"hsa-miR-4524b",
"hsa-miR-4451",
"hsa-miR-4524a",
"hsa-miR-181b",
"hsa-miR-147b",
"hsa-miR-4420",
"hsa-miR-5001",
"hsa-miR-650"
]

allMirnas = downMirnas#upMirnas

figidx = 0

mirna2genes = defaultdict(set)
gene2mirnas = defaultdict(set)


cbn2mirna2evs = defaultdict(lambda: defaultdict(set))

rel2docs = defaultdict(set)

for miRName in allMirnas:


    requestData = {}
    requestData['mirna'] = [miRName]
    requestData['sentences'] = "false"

    requestData['cells'] = [
                {"group": "cells", "name": "aortic smooth muscle cell",  "termid": "CL:0002539"},
                {"group": "cells", "name": "smooth muscle cell",  "termid": "META:83"}
              ]

    #print(requestData)

    json = DataBasePlotter.fetchSimple(requestData )

    for x in json["rels"]:

        genestr = x["lid"]
        mirnastr = x["rid"]

        try:
            origTarget = miRNA(miRName)
            target = miRNA(mirnastr)

            if not origTarget.accept(target, compLevel = miRNACOMPARISONLEVEL.PRECURSOR):
                #print("Not accepted:", target, " as ", origTarget)
                pass

        except:
            #print("skipping", x)
            continue

        #print(x["lid"], x["rid"], simpleStr)


        for ev in x["evidences"]:

            rdv = ev['rel_direction_verb']
            rcat = ev['rel_category']
            rvb = ev['rel_verb']

            #if genestr == "MEG3":
            #    print(genestr, mirnastr, rdv, rcat)

            acceptEvidence = True
            if rdv == '1V2' and rcat in ["CHANGE", "NEG"]:
                acceptEvidence = True
            elif rdv == '2V1' and rcat in ["CHANGE", "POS"]:
                acceptEvidence = True


            #print(genestr, mirnastr, rdv, rcat, rvb, ev["sentence"])

            if acceptEvidence:
                mirna2genes[miRName].add(genestr)
                gene2mirnas[genestr].add(miRName)

            if "docid" in ev:
                rel2docs[miRName].add((ev["docid"], ev["data_source"]))
            else:
                rel2docs[miRName].add(ev["data_source"])

    #print(miRName, len(mirna2genes[miRName]))


for x in mirna2genes:
    for g in mirna2genes[x]:
        print(x, g, sep="\t")


for gene in gene2mirnas:
    if len(gene2mirnas[gene]) > 1:
        print(gene, ";".join(gene2mirnas[gene]), sep="\t")