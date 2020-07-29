import regex
import spacy
import sys, os
from collections import Counter, OrderedDict, defaultdict
from upsetplot import from_contents,plot

import networkx
from natsort import natsorted

from synonymes.GeneOntology import GeneOntology
from synonymes.mirnaID import miRNAPART, miRNA, miRNACOMPARISONLEVEL
from textmining.MirGeneRelCheck import MirGeneRelCheck, SentenceRelationChecker

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

import regex as re

from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher
from utils.tmutils import normalize_gene_names

import pandas as pd
from pandas.plotting import parallel_coordinates
import matplotlib.pyplot as plt


allMirnas = ["miR-126", "miR-200"]

figidx = 0

mirna2genes = defaultdict(set)
gene2mirnas = defaultdict(set)


cbn2mirna2evs = defaultdict(lambda: defaultdict(set))

rel2docs = defaultdict(set)
#nlp = spacy.load('/mnt/d/spacy/models/en_core_web_lg-2.2.0/en_core_web_lg/en_core_web_lg-2.2.0/')  # create blank Language class #en_core_web_lg

sentNum = 99000
fout = open("/mnt/d/dev/git/BERT-Relation-Extraction/data/SCAImirnas/corpus_train/ADD_CHEM_MIREXPLORE.TXT", 'w')


relChecker = SentenceRelationChecker()



for miRName in allMirnas:


    requestData = {}
    requestData['mirna'] = [miRName]
    requestData['sentences'] = "true"
    #print(requestData)

    json = DataBasePlotter.fetchSimple(requestData )

    allDocEvidences = defaultdict(lambda: defaultdict(list))

    for x in json["rels"]:

        genestr = x["lid"]
        mirnastr = x["rid"]
        accMir = None

        try:
            origTarget = miRNA(miRName)
            target = miRNA(mirnastr)

            accMir = origTarget

            if not origTarget.accept(target, compLevel = miRNACOMPARISONLEVEL.MATUREID):
                print("Not accepted:", target, " as ", origTarget)
                pass

        except:
            #print("skipping", x)
            continue

        #print(x["lid"], x["rid"], simpleStr)

        if accMir == None:
            continue



        for ev in x["evidences"]:

            if not "lontid" in ev:
                continue

            allDocEvidences[ev['lontid']][ev['docid']].append(ev)


    for gene in allDocEvidences:

        if len(allDocEvidences[gene]) > 5:

            for docID in allDocEvidences[gene]:

                printedElems = set()

                for ev in allDocEvidences[gene][docID]:

                    sentence = ev['sentence']

                    if ev["ltype"] == "gene":
                        genePos = list(ev["lpos"])
                        mirPos = list(ev["rpos"])
                    else:
                        genePos = list(ev["rpos"])
                        mirPos = list(ev["lpos"])

                    evTuple = (sentence, tuple(genePos), tuple(mirPos))

                    if evTuple in printedElems:
                        continue

                    printedElems.add(evTuple)

                    relRes = relChecker.check_sentence(sentence
                                                       , {"entity_type": "mirna", "entity_type_token": "e1", "entity_location": genePos}
                                                       , {"entity_type": "gene", "entity_type_token": "e2", "entity_location": mirPos}
                                                       )

                    fullsentence = relRes['full_sentence']
                    acceptInteraction = relRes['accept_relation']


                    assert(fullsentence != None)

                    if not fullsentence.endswith("."):
                        fullsentence += "."

                    assert (all(["<e1>" in fullsentence, "<e2>" in fullsentence, "</e1>" in fullsentence, "</e2>" in fullsentence]))


                    # print(ent1, ent2)
                    print(sentNum, "\"" + fullsentence + "\"", sep="\t", file=fout)

                    if acceptInteraction:
                        print("interaction({},{})".format(relRes["entity1"]["entity_type_token"], relRes["entity2"]["entity_type_token"]), file=fout)
                    else:
                        print("no_interaction", file=fout)

                    if relRes["entity1"]["entity_type_token"] == "e1":
                        print("Comment: {}:::{}".format(relRes["entity1"]["entity_text"], relRes["entity2"]["entity_text"]), file=fout)
                    else:
                        print("Comment: {}:::{}".format(relRes["entity2"]["entity_text"], relRes["entity1"]["entity_text"]), file=fout)

                    print(file=fout)



