import pickle
import tempfile

from collections import defaultdict, Counter
import networkx as nx
import scipy.stats

from networkx.drawing.nx_agraph import graphviz_layout, pygraphviz_layout
import os, sys

from synonymes.GeneOntology import GeneOntology
from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from synonymes.mirnaID import miRNA, miRNAPART, miRNACOMPARISONLEVEL
from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    cellObo = GeneOntology("/mnt/d/owncloud/data/miRExplore/obodir/meta_cells.obo")

    cellTypeName2Terms = {
        "EC": ["META:52"],
        "MC": ["META:148", "META:99"],
        "FC": ["CL:0000891"],
        "SMC": ["META:83"],
    }

    cellType2AccTerms = {}
    for cellT in cellTypeName2Terms:

        cellType2AccTerms[cellT] = set()

        for et in cellTypeName2Terms[cellT]:

            oboT = cellObo.getID(et)

            if oboT != None:
                cellType2AccTerms[cellT].add(et)
                for x in oboT.getAllChildren():
                    cellType2AccTerms[cellT].add(x.termid)

            else:
                print("No such obo term:", et)


    missRates = sorted([2/19, 5/17, 1/4, 5/12, 0/6, 2/8, 3/12, 1/23])

    print("missrate", sum(missRates)/len(missRates))

    print(scipy.stats.describe(missRates))
    print(np.median(missRates))
    #exit(0)

    largeInteractions = {
        'CCL9': ['miR-30d-3p', 'miR-3473c', 'let-7g-5p'],
        'CXCL5': ['miR-204-5p', 'let-7g-5p', 'miR-362-3p', 'miR-155-5p'],
        'CXCL1': ['miR-194-2-3p', 'miR-128-3p', 'miR-194-5p', 'miR-199b-5p', 'miR-467g', 'miR-122-5p'],
        'CXCL13': ['miR-122-5p'],
        'CXCL14': ['miR-301b-3p'],
        'CXCR2': ['let-7g-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7a-5p', 'let-7i-5p', 'miR-98-5p'],
        'CXCL7': ['let-7g-5p'],
        'CCL2': ['let-7a-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7g-5p', 'let-7i-5p', 'miR-181a-5p'],
        'CXCL9': ['miR-1935'],
        'CCL3': ['miR-30a-5p','miR-30b-5p','miR-30c-5p','miR-30d-5p','miR-30e-5p'],
        'CCL7': ['miR-181a-5p', 'miR-322-5p', 'miR-29a-5p', 'miR-29b-1-5p'],
        'CCL22': [  'miR-34a-5p'],
        'CXCL10': ['miR-503-3p', 'miR-186-5p'],
        'CCR5': ['miR-186-5p', 'miR-669j', 'miR-21-5p', 'miR-146a-5p', 'miR-150-5p', 'miR-146b-5p', 'miR-669k-3p', 'miR-142-3p', 'miR-34a-5p'],
        'CCL4': ['miR-27b-3p', 'miR-27a-3p', 'miR-21-3p', 'miR-467f'],
        'CX3CL1': ['miR-15a-5p', 'miR-322-5p', 'miR-706', 'miR-762', 'miR-665-3p', 'miR-758-3p', 'miR-381-3p'],
        'CXCR4': ['miR-381-3p', 'miR-21-3p', 'miR-467a-5p', 'miR-467h', 'miR-218-5p', 'miR-1a-3p', 'miR-181d-5p',
                  'miR-206-3p', 'miR-181b-5p', 'miR-9-5p', 'miR-132-3p', 'miR-25-3p', 'miR-467d-5p', 'miR-669k-3p',
                  'miR-146b-5p', 'miR-467b-5p', 'miR-467e-5p', 'miR-467f', 'miR-146a-5p'],
        'CCR7': ['let-7g-5p', 'miR-23b-3p', 'miR-669p-5p', 'miR-23a-5p', 'let-7e-5p', 'miR-669l-5p', 'miR-15a-5p', 'miR-467e-5p', 'miR-21-5p', 'miR-16-5p', 'let-7d-5p', 'miR-669n', 'miR-98-5p', 'let-7b-5p', 'let-7a-5p', 'let-7i-5p', 'let-7c-5p', 'miR-15b-5p', 'miR-467h'],
        'CXCL12': [
            'miR-532-5p', 'miR-130b-3p', 'miR-222-3p', 'miR-144-3p', 'miR-542-3p', 'miR-149-5p', 'miR-330-3p', 'miR-532-3p', 'miR-3470b', 'miR-125b-5p', 'miR-221-3p', 'miR-19b-3p', 'miR-301b-3p',
            'miR-34b-5p', 'miR-125a-3p', 'miR-126-3p', 'miR-16-1-3p', 'miR-882', 'miR-497-5p', 'miR-26a-5p', 'miR-124-3p', 'miR-26b-5p', 'miR-5620-3p', 'mIR-19a-3p', 'miR-130a-3p', 'miR-690',
            'miR-185-5p', 'miR-31-5p', 'miR-340-5p', 'miR-1843-5p', 'miR-466f-3p', 'miR-301a-3p', 'miR-101a-3p', 'miR-210-3p', 'miR-107-3p', 'miR-706', 'miR-23b-3p', 'miR-146a-5p', 'miR-467f',
            'miR-322-5p', 'miR-15a-5p', 'miR-29b-1-5p', 'let-7e-5p', 'miR-23a-3p', 'miR-338-3p', 'miR-103-3p', 'miR-362-3p', 'let-7g-5p', 'miR-155-5p', 'miR-140-5p', 'miR-122-5p', 'miR-22-3p', 'miR-3470a', 'let-7d-5p'
        ]

    }

    inflammatory_ec = {
        'CCL2': ['miR-495'],
        'PPARA': ['miR-21'],
        'THBS1': ['let-7g'],
        'TGFBR1': ['let-7g'],
        'SMAD2': ['let-7g'],
        'SIRT1': ['let-7g', 'miR-34a'],
        'SOCS5': ['miR-92a'],
        'KLF2': ['miR-92a'],
        'KLF4': ['miR-92a', 'miR-103'],
        'TRC': ['miR-10a'],
        'TAK1': ['miR-10a'],
        'KPNA3': ['miR-181b'],
        'TRAF6': ['miR-146'],
        'IRAK1': ['miR-146'],
        'HUR': ['miR-146'],
        'CXCL12': ['miR-126-3p'],
        'RGS16': ['miR-126-3p'],
        'ETS1': ['miR-155', 'miR-221', 'miR-222'],
        'TNPO2': ['miR-181b']
    }

    ccl2_macrophage = {
        'CCL2': ['miR-124a', 'miR-150'],
        'CHI3L1': ['miR-24'],
        'TLR4': ['miR-146a'],
        'TRAF6': ['miR-146'],
        'IRAK1': ['miR-146'],
        'AKT1': ['miR-342-5p'],
        'BCL6': ['miR-155'],
        'LPL': ['miR-467b']
    }

    normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")

    andreouInts1= {'ICAM1': ['miR-223', 'miR-17-3p'], 'SELE': ['miR-31'], 'VCAM1': ['miR-126', 'miR-126-5p'], 'DLK1': ['miR-126', 'miR-126-5p'], 'JAMA': ['miR-143', 'miR-145'], 'SIRT1': ['miR-217'], 'NOS3': ['miR-155'], 'ETS1': ['miR-155', 'miR-126', 'miR-126-5p'], 'PPARA': ['miR-21'], 'TIMP3': ['miR-712'], 'MAP3K7': ['miR-10'], 'KPNA4': ['miR-181b'], 'TNPO2': ['miR-181b'], 'TRAF6': ['miR-146a', 'miR-146b'], 'IRAK1': ['miR-146a', 'miR-146b'], 'IRAK2': ['miR-146a', 'miR-146b'], 'KLF2': ['miR-143', 'miR-145', 'miR-126', 'miR-126-5p', 'miR-92a'], 'CXCL12': ['miR-126', 'miR-126-5p'], 'KLF4': ['miR-92a', 'miR-663'], 'SOCS5': ['miR-92a']}
    andreouInts2= {'IFNG': ['miR-29', 'miR-155'], 'CD40': ['miR-146a', 'miR-181a'], 'CD40L': ['miR-146a', 'miR-181a'], 'MMP2': ['miR-21', 'miR-21', 'miR-24', 'miR-29'], 'MMP3': ['miR-21', 'miR-24', 'miR-29'], 'MMP8': ['miR-21', 'miR-24', 'miR-29'], 'MMP9': ['miR-21', 'miR-24', 'miR-29'], 'MMP12': ['miR-21', 'miR-24', 'miR-29'], 'MMP13': ['miR-21', 'miR-24', 'miR-29'], 'MMP14': ['miR-21', 'miR-24', 'miR-29']}
    andreouInts3= {'MTTP': ['miR-30c'], 'ABCA1': ['miR-33', 'miR-302a'], 'ABCG1': ['miR-33'], 'CPT1A': ['miR-33'], 'KLF2': ['miR-92a'], 'KLF4': ['miR-92a', 'miR-145'], 'SOCS5': ['miR-92a'], 'DLK1': ['miR-126-5p'], 'RGS16': ['miR-126-3p'], 'BCL6': ['miR-155'], 'SOCS1': ['miR-155'], 'MAP3K10': ['miR-155'], 'KPNA4': ['miR-181b'], 'TNPO2': ['miR-181b'],'AKT1': ['miR-342-5p'], 'TIMP3': ['miR-712']}


    andreouInteractions1 = {}
    andreouInteractions2 = {}
    andreouInteractions3 = {}

    for gene in andreouInts1:
        andreouInteractions1[normGeneSymbols.get(gene, gene)] = andreouInts1[gene]

    for gene in andreouInts2:
        andreouInteractions2[normGeneSymbols.get(gene, gene)] = andreouInts2[gene]

    for gene in andreouInts3:
        andreouInteractions3[normGeneSymbols.get(gene, gene)] = andreouInts3[gene]

    networks = {}
    networks['large_chemokines'] = largeInteractions
    networks['large_chemokines_cv'] = largeInteractions
    networks['large_chemokines_athero'] = largeInteractions


    networks['andreou_fig1'] = andreouInteractions1
    networks['andreou_fig1_cv'] = andreouInteractions1
    networks['andreou_fig1_athero'] = andreouInteractions1

    networks['andreou_fig2'] = andreouInteractions2
    networks['andreou_fig2_cv'] = andreouInteractions2
    networks['andreou_fig2_athero'] = andreouInteractions2

    networks['andreou_table1'] = andreouInteractions3
    networks['andreou_table1_cv'] = andreouInteractions3
    networks['andreou_table1_athero'] = andreouInteractions3

    networks['inflammatory_ec'] = inflammatory_ec
    networks['inflammatory_ec_cv'] = inflammatory_ec
    networks['inflammatory_ec_athero'] = inflammatory_ec

    networks['macrophages'] = ccl2_macrophage
    networks['macrophages_cv'] = ccl2_macrophage
    networks['macrophages_athero'] = ccl2_macrophage
    networks['macrophages_all'] = ccl2_macrophage

    networkGraphs = {}
    makeStory = [
        ('macrophages', 'macrophages_cv', 'macrophages_athero', 'macrophages_all'),
        ('inflammatory_ec','inflammatory_ec_cv', 'inflammatory_ec_athero'),
        ('large_chemokines', 'large_chemokines_cv', 'large_chemokines_athero'),
        ('andreou_fig1', 'andreou_fig1_cv', 'andreou_fig1_athero'),
        ('andreou_fig2', 'andreou_fig2_cv', 'andreou_fig1_athero'),
        ('andreou_table1', 'andreou_table1_cv', 'andreou_table1_athero'),
    ]

    ignoreNetworks = []#['large_chemokines']

    networkRestrictions = {
        'inflammatory_ec':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
                ]
            },
        'inflammatory_ec_cv':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
                ],
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },

        'inflammatory_ec_athero':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
                ],
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
                ]
            },
        'macrophages':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "macrophage", "termid": "META:99"}
                ]
            },
        'macrophages_all':
            {
                'sentences': "false",
            },
        'macrophages_cv':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "macrophage", "termid": "META:99"}
                ],
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },
        'macrophages_athero':
            {
                'sentences': "false",
                "cells": [
                    {"group": "cells", "name": "macrophage", "termid": "META:99"}
                ],
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
                ]
            },
        'large_chemokines': {
            'sentences': "false",
        },
        'large_chemokines_cv': {
            'sentences': "false",
            "disease": [
                {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
            ]        },
        'large_chemokines_athero': {
            'sentences': "false",
            "disease": [
                {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
            ]},
        'andreou_fig1': {
            'sentences': "false",
            "cells": [
                {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
            ],
        },
        'andreou_fig2': {
            'sentences': "false",
        },
        'andreou_table1': {
            'sentences': "false",
        },

        'andreou_fig1_cv':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ],
                "cells": [
                    {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
                ],
            },
        'andreou_fig2_cv':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },
        'andreou_table1_cv':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },

        'andreou_fig1_athero':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'} #{'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ],
                "cells": [
                    {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
                ],
            },
        'andreou_fig2_athero':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
                ]
            },
        'andreou_table1_athero':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}
                ]
            },
    }

    restrictDF = DataFrame()
    restrictDF.addColumns(["Network", "Cells", "Disease", "Other"], "")

    for x in networkRestrictions:

        restricts = networkRestrictions[x]

        networkDRdict = defaultdict(str)
        networkDRdict["Network"] =  x.replace("_", " ")

        diseaseElems = []
        cellElems = []
        otherElems = []

        for restrictType in restricts:

            if restrictType == "sentences":
                continue

            if restrictType in ["disease"]:
                for elem in restricts[restrictType]:
                    diseaseElems.append( elem['name'] + " ("+elem['termid']+")")

            elif restrictType in ["cells"]:
                for elem in restricts[restrictType]:
                    cellElems.append( elem['name'] + " ("+elem['termid']+")")

            else:
                for elem in restricts[restrictType]:
                    otherElems.append( elem['name'] + " ("+elem['termid']+")")


        networkDRdict['Cells'] =  "\makecell[l]{" +"\\\\".join(sorted(cellElems)) + "}"
        networkDRdict['Disease'] = "\makecell[l]{" + "\\\\".join(sorted(diseaseElems)) + "}"
        networkDRdict['Other'] = "\makecell[l]{" + "\\\\".join(sorted(otherElems)) + "}"

        dr = DataRow.fromDict(networkDRdict)
        restrictDF.addRow(dr)

    print(restrictDF._makeLatex())

    mirna2cellOut = open("/mnt/d/yanc_network/important_networks.txt", 'w')


    def acceptEvidence(ev):

        return True


    figidx = 0

    for network in networks:
        figidx+= 1

        networkGraph = nx.Graph()

        if network in ignoreNetworks:
            continue

        #if not network in ['large_chemokines_athero']: #'large_chemokines', 'large_chemokines_cv',
        #    continue

        interactions = networks[network]
        acceptedInteractions = defaultdict(set)
        typeByGene = defaultdict(lambda: Counter())
        elemsByGene = defaultdict(lambda: defaultdict(set))

        allMirna = set()
        miStr2mirna = {}

        normalizedInteractions = defaultdict(set)

        for gene in interactions:

            for mirna in interactions[gene]:

                try:
                    oMirna = miRNA(mirna)
                    mirnaN = oMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])
                    mirna = mirnaN

                except:
                    pass

                allMirna.add(mirna)
                normalizedInteractions[gene].add(mirna)

        allTargetMirna = []

        for x in allMirna:
            try:

                oMirna = miRNA(x)

                allTargetMirna.append( oMirna )
                miStr2mirna[x] = oMirna

            except:
                pass


        requestData = None


        if network in networkRestrictions:
            requestData = networkRestrictions[network]
        else:
            requestData = {
                            'sentences': "false",
                            }

        requestData["gene"]= [x for x in interactions]
        #requestData['disease'] = [{'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'}]


        graph, nodeCounter, edge2datasourceCount, jsonRes = DataBasePlotter.fetchGenes(requestData, minPMIDEvCount=0, minTgtCount=0, acceptEv=acceptEvidence)

        print(len(jsonRes['rels']))

        htmlDF = DataFrame()
        htmlDF.addColumns(['gene rel', 'gene', 'miRNA Group', 'miRNA', 'Original Network', 'PubMed', 'MIRECORD', 'MIRTARBASE', 'DIANA', 'Disease', 'Cells', 'GO'])

        foundGenes = set()

        interactionCountExpDB = 0
        interactionCountPubmed = 0
        interactionCountPubmedExpDB = 0

        interactionsWithPubmed = set()
        interactionsWithExpDB = set()


        interactionCountAccepted = 0
        interactionCountAdditional = 0
        interactionCountMissing = 0

        expDBEdge2Source = defaultdict(set)
        mirnaGeneSourceCounter = Counter()
        mirnaGeneIntSourceCounter = Counter()

        intToSources = defaultdict(set)

        simpleMirna2Genes = defaultdict(set)
        gene2simpleMirna = defaultdict(set)

        distinctMirna2Genes = defaultdict(set)
        gene2distinctMirna = defaultdict(set)

        simpleInteractions = set()
        distinctInteractions = set()

        pmidEvidences = set()

        for rel in jsonRes['rels']:

            orderedEdge = [None, None]

            if rel['ltype'] == "gene":
                orderedEdge[0] = rel['lid']
            elif rel['ltype'] == "mirna":
                orderedEdge[1] = rel['lid']

            if rel['rtype'] == "gene":
                orderedEdge[0] = rel['rid']
            elif rel['rtype'] == "mirna":
                orderedEdge[1] = rel['rid']

            if orderedEdge[0] == "TNPO2":
                orderedEdge[0] = "KPNA4"

            orderedEdges = set()

            if orderedEdge[1].startswith("microRNAS"):
                continue

            wasAccepted = False
            """
            for tMirna in allTargetMirna:

                if tMirna.accept(orderedEdge[1]):

                    wasAccepted = True

                    orderedEdges.add(
                        (orderedEdge[0], str(tMirna))
                    )
            """

            if not wasAccepted:
                orderedEdges.add(tuple(orderedEdge))


            for oEdge in orderedEdges:

                foundGenes.add(oEdge[0])
                origEdge = tuple(oEdge)

                try:
                    miObj = miRNA(oEdge[1])
                    miStr = miObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

                    selfObj = miRNA(miStr)

                    oEdge = list(oEdge)
                    oEdge[1] = miStr
                    oEdge = tuple(oEdge)
                except:
                    print("Ignoring edge for mirna", origEdge, oEdge)
                    continue
                    #pass

                edgeStatus = None

                allGeneMirna = normalizedInteractions[oEdge[0]]
                miAccepted = False

                allAcceptedStr = set()

                if network == "macrophages_athero":
                    print(oEdge)

                for strMirna in allGeneMirna:
                    miObj = miStr2mirna.get(strMirna, None)

                    if miObj == None:
                        continue

                    miObjAccepts = miObj.accept(oEdge[1], compLevel=miRNACOMPARISONLEVEL.PRECURSOR)
                    miAccepted = miAccepted or miObjAccepts

                    if "fig1_athero" in network:
                        if oEdge[0] == "IRAK1":
                            print(oEdge)
                            miObjAccepts = miObj.accept(oEdge[1], compLevel=miRNACOMPARISONLEVEL.PRECURSOR)


                    if miObjAccepts:
                        acceptedInteractions[oEdge[0]].add(strMirna)
                        edgeStatus = "accepted"

                        allAcceptedStr.add(strMirna)

                        #print(oEdge[0], oEdge[1], strMirna)

                if network == "macrophages_athero":
                    print(allAcceptedStr)

                if not miAccepted:
                    edgeStatus = "additional"


                if edgeStatus == 'accepted':
                    interactionCountAccepted += 1
                else:
                    interactionCountAdditional += 1


                # determine kindness of interaction
                edgeResources = edge2datasourceCount.get(origEdge, None)


                edgeLabel = ""
                if edgeResources != None:

                    hasPMID = False
                    hasExpDB = False

                    if edgeResources['pmid'] > 0:
                        interactionCountPubmed += 1
                        interactionsWithPubmed.add(oEdge)
                        hasPMID = True
                        edgeLabel += "p"

                    if edgeResources['DIANA'] > 0 or edgeResources['miRTarBase'] > 0 or edgeResources['miRecords'] > 0:

                        if edgeResources['DIANA'] > 0:
                            edgeLabel += "D"

                        if edgeResources['miRTarBase'] > 0:
                            edgeLabel += "T"

                        if edgeResources['miRecords'] > 0:
                            edgeLabel += "R"

                        if not "D" in edgeLabel and not "T" in edgeLabel and not "R" in edgeLabel:
                            print(edgeResources)

                        interactionCountExpDB += 1
                        interactionsWithExpDB.add(oEdge)
                        hasExpDB = True


                    if hasPMID and hasExpDB:
                        interactionCountPubmedExpDB += 1


                objMirna = miRNA(oEdge[1])
                simpleMirna = objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])
                simpleEdge = tuple([oEdge[0], simpleMirna])


                simpleMirna2Genes[simpleMirna].add(oEdge[0])
                gene2simpleMirna[oEdge[0]].add(simpleMirna)

                distinctMirna2Genes[oEdge[1]].add(oEdge[0])
                gene2distinctMirna[oEdge[0]].add(oEdge[1])

                networkGraph.add_edge(oEdge[0], oEdge[1], edgeResources=edgeResources, label=edgeLabel, color= 'g' if edgeStatus == "accepted" else "b")

                typeByGene[oEdge[0]][edgeStatus] += 1
                elemsByGene[oEdge[0]][edgeStatus].add( oEdge[1] )

                simpleInteractions.add(simpleEdge)
                distinctInteractions.add(oEdge)


                pmidEvs = set()
                mirtarbaseEvs = set()
                mirecordsEvs = set()
                dianaEvs = set()

                interactionEvidences = set()
                docDiseases = []
                docCells = []
                docGOs = []

                for ev in rel['evidences']:

                    if 'docid' in ev:
                        docid = ev['docid']

                        disEvs = jsonRes['pmidinfo'].get('disease', {}).get(docid, {})

                        for disEv in disEvs:
                            did = disEv['termid']
                            dname = disEv['termname']

                            docDiseases.append((dname, did, docid))

                        cellEvs = jsonRes['pmidinfo'].get('cells', {}).get(docid, {})

                        for cellEv in cellEvs:
                            did = cellEv['termid']
                            dname = cellEv['termname']

                            docCells.append((dname, did, docid))

                            for ct in cellType2AccTerms:

                                ctTerms = cellType2AccTerms[ct]

                                if did in ctTerms:
                                    print(network, oEdge[0], oEdge[1], ct, docid, sep="\t", file=mirna2cellOut)

                        goEvs = jsonRes['pmidinfo'].get('go', {}).get(docid, {})

                        for goEv in goEvs:
                            did = goEv['termid']
                            dname = goEv['termname']

                            docGOs.append((dname, did, docid))

                    if ev['data_source'] == "DIANA":
                        dianaEvs.add( (ev['method'], (ev['direction'])))

                        expDBEdge2Source[oEdge].add(('DIANA', ev['method']))

                        mirnaGeneSourceCounter['DIANA'] += 1
                        interactionEvidences.add("DIANA")

                        intToSources[simpleEdge].add("DIANA")

                    elif ev['data_source'] == "miRTarBase":
                        mirtarbaseEvs.add(
                            (ev['data_id'], ",".join(ev['exp_support']), ev['functional_type'], ev['docid'])
                        )
                        expDBEdge2Source[oEdge].add(('MIRTARBASE', (ev['data_id'], ",".join(ev['exp_support']), ev['functional_type'], ev['docid'])))
                        mirnaGeneSourceCounter['MIRTARBASE'] += 1
                        interactionEvidences.add("MIRTARBASE")
                        intToSources[simpleEdge].add("MIRTARBASE")

                        if ev['docid'] != None:
                            pmidEvidences.add(ev['docid'])



                    elif ev['data_source'] == "pmid":
                        pmidEvs.add((ev['docid'],))

                        expDBEdge2Source[oEdge].add(('PMID', (ev['docid'])))
                        mirnaGeneSourceCounter['PMID'] += 1
                        interactionEvidences.add("PMID")
                        intToSources[simpleEdge].add("PMID")

                        pmidEvidences.add(ev['docid'])

                    elif ev['data_source'] == "mirecords":
                        mirecordsEvs.add((ev['docid']))
                        expDBEdge2Source[oEdge].add(('MIRECORDS', (ev['docid'])))
                        mirnaGeneSourceCounter['mirecords'] += 1
                        interactionEvidences.add("mirecords")
                        intToSources[simpleEdge].add("mirecords")

                        if ev['docid'] != None:
                            pmidEvidences.add(ev['docid'])

                    else:
                        print("Unhandled data source", ev['data_source'])

                interactionEvidences = tuple(sorted(interactionEvidences))
                mirnaGeneIntSourceCounter[interactionEvidences] += 1


                dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
                    geneCap=oEdge[0].upper(), geneLow=oEdge[1].upper())


                pmidStr = "<br/>".join(
                    [
                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}\" target=\"_blank\">{pmid}</a>".format(pmid=elem[0]) for elem in pmidEvs
                    ]
                )

                mirtarbaseStr = "<br/>".join(
                    [
                        "<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid={mtbid}\">{mtbid}</a>".format(mtbid=elem[0]) for elem in mirtarbaseEvs
                    ]
                )

                mirecordStr = "<br/>".join(
                    [
                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}\" target=\"_blank\">{pmid}</a>".format(pmid=elem[0]) for elem in mirecordsEvs
                    ]
                )

                dianaStr = "<br/>".join(
                    [
                        "{method} {direction}".format(method=elem[0], direction=elem[1]) for elem in dianaEvs
                    ]
                )

                goStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for elem in docGOs
                    ]
                )

                cellStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for elem in docCells
                    ]
                )

                diseaseStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for elem in docDiseases
                    ]
                )

                addRow = {
                    'gene rel': oEdge[0] + "<br/>" + oEdge[1],
                    'gene': oEdge[0],
                    'miRNA Group': objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                    'miRNA': "<br/>".join(allAcceptedStr),
                    'Original Network': "{edgestate}</br>".format(edgestate=edgeStatus) +
                                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term={miRes}+{miShort}\">Search PUBMED</a>".format(miRes=oEdge[1], miShort=objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID]))+
                                        "</br><a href=\"{dianaLink}\">Search DIANA</a>".format(dianaLink=dianaLink)
                    ,
                    'PubMed': pmidStr,
                    'MIRECORD': mirecordStr,
                    'MIRTARBASE': mirtarbaseStr,
                    'DIANA': dianaStr,
                    "Cells": cellStr,
                    "Disease": diseaseStr,
                    "GO": goStr
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)

        for gene in normalizedInteractions:

            for mirna in normalizedInteractions[gene]:

                edgeWasFound = mirna in acceptedInteractions[gene]

                if edgeWasFound:
                    continue

                edgeStatus = "missed"

                interactionCountMissing += 1

                networkGraph.add_edge(gene, mirna, color='r')


                typeByGene[gene][edgeStatus] += 1
                elemsByGene[gene][edgeStatus].add(mirna)

                objMirna = miRNA(mirna)
                dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
                    geneCap=gene.upper(), geneLow=gene.capitalize())

                cellStr = ""
                diseaseStr = ""
                goStr = ""

                addRow = {
                    'gene rel': gene,
                    'gene': gene,
                    'miRNA Group': objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                    'miRNA': mirna,
                    'Original Network': "{edgestate}</br>".format(edgestate=edgeStatus) +
                                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term={gene} {miRes}+{miShort}\">Search PUBMED</a>".format(gene=gene,miRes=mirna, miShort=objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID]))+
                                        "</br><a href=\"{dianaLink}\">Search DIANA</a>".format(dianaLink=dianaLink)
                    ,
                    'PubMed': "",
                    'MIRECORD': "",
                    'MIRTARBASE': "",
                    'DIANA': "",
                    "Cells": cellStr,
                    "Disease": diseaseStr,
                    "GO": goStr
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)


        allDegrees = networkGraph.degree()


        singleMirnaGenes = set()
        multiMirnaGenes = set()

        singleGeneMirnas = set()
        multiGeneMirnas = set()



        """
        THIS ANALYSIS ALSO HAS ADDED EDGES FROM MISSING!!!!!!
        """
        for node, nodedegree in allDegrees:

            if node.startswith("miR-") or node.startswith("let-"):

                if nodedegree == 1:
                    singleGeneMirnas.add(node)
                elif nodedegree >= 1:
                    multiGeneMirnas.add(node)

            else:

                if nodedegree == 1:
                    singleMirnaGenes.add(node)
                elif nodedegree >= 1:
                    multiMirnaGenes.add(node)

        print(network)


        print("Simple Interactions", len(simpleInteractions))
        print("Distinct Interactions", len(distinctInteractions))

        print("Number of simple mirnas", len(simpleMirna2Genes))
        print("Number of simple mirnas with 1 gene", sum([1 for x in simpleMirna2Genes if len(simpleMirna2Genes[x]) == 1]))
        print("Number of simple mirnas > 1 gene",
              sum([1 for x in simpleMirna2Genes if len(simpleMirna2Genes[x]) > 1]))

        print("Number of genes", len(gene2simpleMirna))
        print("Number of genes with 1 simple mirna",
              sum([1 for x in gene2simpleMirna if len(gene2simpleMirna[x]) == 1]))
        print("Number of genes with >1 simple mirna",
              sum([1 for x in gene2simpleMirna if len(gene2simpleMirna[x]) > 1]))
        # anzahl gene mit nur einer miRNA
        print("Number of genes with only one miRNA", len(singleMirnaGenes))


        print("Number of distinct mirnas", len(distinctMirna2Genes))
        print("Number of distinct mirnas with 1 gene", sum([1 for x in distinctMirna2Genes if len(distinctMirna2Genes[x]) == 1]))
        print("Number of distinct mirnas > 1 gene",
              sum([1 for x in distinctMirna2Genes if len(distinctMirna2Genes[x]) > 1]))

        print("Number of genes", len(gene2distinctMirna))
        print("Number of genes with 1 distinct mirna",
              sum([1 for x in gene2distinctMirna if len(gene2distinctMirna[x]) == 1]))
        print("Number of genes with >1 distinct mirna",
              sum([1 for x in gene2distinctMirna if len(gene2distinctMirna[x]) > 1]))
        # anzahl gene mit nur einer miRNA
        print("Number of genes with only one miRNA", len(singleMirnaGenes))

        # anzahl gene mit mehr als einer miRNA
        print("Number of genes with more than one miRNA", len(multiMirnaGenes))

        print("Number of miRNAs with only one genes", len(singleGeneMirnas))
        print("Number of miRNAs with more than one genes", len(multiGeneMirnas))
        print("miRNA Gene Interaction source counter", sum([mirnaGeneSourceCounter[x] for x in mirnaGeneSourceCounter]), mirnaGeneSourceCounter)
        print("miRNA Gene Interactions detailed", sum([mirnaGeneIntSourceCounter[x] for x in mirnaGeneIntSourceCounter]), mirnaGeneIntSourceCounter)

        intToSourcesCount = Counter()
        for ed in intToSources:
            uniedge = tuple(sorted(intToSources[ed]))
            intToSourcesCount[uniedge] += 1

        print("miRNA(simple) Gene Interactions detailed", sum([intToSourcesCount[x] for x in intToSourcesCount]), intToSourcesCount)


        totalInteractionCount = interactionCountAccepted+ interactionCountAdditional+ interactionCountMissing
        totalGraphInteractionCount = interactionCountAccepted  + interactionCountAdditional
        # anzahl kanten in graph
        print("Number of interactions (acc, add, graph, mis, total)", interactionCountAccepted, interactionCountAdditional, totalGraphInteractionCount, interactionCountMissing, totalInteractionCount)

        # anzahl kanten nur durch experimental db
        print("Number of interactions supported by experimental db", len(interactionsWithExpDB))

        printSet = interactionsWithExpDB.difference(interactionsWithPubmed)
        if len(printSet) > 50:
            printSet = {}
        print("Number of interactions supported only by experimental db", len(interactionsWithExpDB.difference(interactionsWithPubmed)), printSet)

        for edge in printSet:
            print(edge, expDBEdge2Source.get(edge, None))


        # anzahl kanten nur durch pubmed
        print("Number of interactions supported by pubmed", len(interactionsWithPubmed))
        print("Number of interactions supported only by pubmed", len(interactionsWithPubmed.difference(interactionsWithExpDB)))


        # anzahl kanten mit pubmed UND experimental db
        printSet = interactionsWithExpDB.intersection(interactionsWithPubmed)
        if len(printSet) > 50:
            printSet = {}
        print("Number of interactions supported by pubmed+expdb", len(interactionsWithExpDB.intersection(interactionsWithPubmed)), printSet)
        for edge in printSet:
            print(edge, expDBEdge2Source.get(edge, None))
        # anzahl kanten nur mit exp hinweisen aber zu weak

        # anzahl kanten nur durch pubmed aber zu weak
        print("Number of PMIDs", len(pmidEvidences))

        print()
        print()
        print()
        print()
        for gene in sorted([x for x in typeByGene]):
            if "athero" in network:
                print(gene, typeByGene[gene], elemsByGene[gene]['additional'], elemsByGene[gene]['missing'])
            else:
                print(gene, typeByGene[gene], elemsByGene[gene]['missing'])

        print()
        print()

        networkDF = DataFrame()
        networkDF.addColumns(["gene", "accepted", "additional", "missing", "missing miRNAs"])


        print(network)
        for gene in sorted([x for x in typeByGene]):

            infoDict = {
                "gene": gene,
                "accepted": typeByGene[gene]["accepted"],
                "additional": typeByGene[gene]["additional"],
                "missing": typeByGene[gene]["missing"],
                "missing miRNAs": ", ".join(elemsByGene[gene]['missing'])
            }

            dr = DataRow.fromDict(infoDict)
            networkDF.addRow(dr)

            print("Gene:", gene, "Status: ", ", ".join([": ".join([x, str(typeByGene[gene][x])]) for x in typeByGene[gene]]), "Missing miRNAs: "+",".join(elemsByGene[gene]['missing']))

        print()
        print()

        print(networkDF._makeLatex())

        print()
        print()


        networkGraphs[network] = networkGraph

        htmlDF.export("/mnt/d/yanc_network/" + network.replace(" ", "_") + ".html", ExportTYPE.HTML)
        htmlDF.export("/mnt/d/yanc_network/" + network.replace(" ", "_") + ".tsv", ExportTYPE.TSV)


    print("Network covered genes", network, foundGenes)


    figidx = 0
    for stages in makeStory:

        mergedGraph = networkGraphs[stages[0]]

        for i in range(1, len(stages)):
            mergedGraph = nx.compose(mergedGraph, networkGraphs[stages[i]])


        hasLargeStage = any(['large' in stage for stage in stages])

        pos = nx.spring_layout(mergedGraph)

        for stage in stages:


            networkGraph = networkGraphs[stage]

            edges = networkGraph.edges()
            colors = [networkGraph[u][v]['color'] for u, v in edges]

            d = nx.degree(networkGraph)
            nodes = networkGraph.nodes()
            nodeColors = []
            nodeSizes = []

            allSizes = [x[1] for x in d]
            minSize = min(allSizes)
            maxSize = max(allSizes)
            diffSize = maxSize-minSize

            fontSize = 16
            minNodeSize = 1200
            figSize = (20, 14)
            edgeWidth = 3

            if hasLargeStage:
                fontSize = 8
                minNodeSize = 100
                figSize = (20,30)
                edgeWidth = 0.75


            plt.figure(figidx, figsize=figSize)
            figidx += 1

            maxNodeSize = 3000
            diffNodeSize = maxNodeSize-minNodeSize
            nodeList = []






            for x in nodes:
                if any([x.lower().startswith(y) for y in ['mir', 'let']]):
                    nodeColors.append('blue')
                else:
                    nodeColors.append('green')

                nodeDegree = d[x]

                nodeDegree -= minSize
                nodeDegree = nodeDegree / diffSize

                nodeSize = minNodeSize + diffNodeSize * nodeDegree

                nodeSizes.append( nodeSize )
                nodeList.append(x)

            nx.draw(networkGraph, pos, font_size=fontSize, with_labels=False, node_color=nodeColors, edges=edges, edge_color=colors, nodelist=nodeList, node_size=nodeSizes, width=edgeWidth, font_weight='bold', dpi=1000)


            for p in pos:  # raise text positions

                clist = list(pos[p])

                if p in nodeList:
                    if nodeSizes[nodeList.index(p)] < 1000:
                        clist[1] = clist[1] + 0.005
                    else:
                        clist[1] = clist[1] + 0.02
                pos[p] = tuple(clist)

            nx.draw_networkx_labels(networkGraph, pos, font_weight='bold', font_size=fontSize)

            plt.suptitle(stage)

            plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".png")
            plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".pdf")

            CytoscapeGrapher.showGraph(networkGraph, location="/mnt/d/yanc_network/", name="cyjs_" + stage.replace(" ", "_"), title=stage)
            CytoscapeGrapher.exportGraphML(networkGraph, location="/mnt/d/yanc_network/", name="cyjs_" + stage.replace(" ", "_"))

            CytoscapeGrapher.export_d3Json(networkGraph, location="/mnt/c/Users/mjopp/Desktop/d3test/",name= "cyjs_" + stage.replace(" ", "_"), singletons=True )


            if stage == "large_chemokines":

                with open("/mnt/c/Users/mjopp/Desktop/gergely.tsv", "w") as fgout:

                    existingEdges = set()

                    for gene in largeInteractions:

                        for mirna in largeInteractions[gene]:
                            try:
                                miObj = miRNA(mirna)
                                miStr = miObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])
                                selfObj = miRNA(miStr)

                                existingEdges.add((gene, miStr))
                            except:
                                continue


                    for node in exportGraph.nodes():
                        nodename = str(node).lower()

                        if nodename.startswith("mir") or nodename.startswith("let"):
                            fmtNode = "<node y=\"20\" x=\"360\" shape=\"ellipse\" color=\"FFFFFF\" label=\"{elemid}\" id=\"{elemid}\" bgfixed=\"-1\"/>".format(elemid=node)
                            #print(fmtNode)

                    for edge in exportGraph.edges():
                        fmtEdge = "<edge dir=\"false\" bold=\"false\" color=\"000000\" width=\"null\" inhibit=\"false\" rev=\"false\" arrowcolor=\"000000\" to=\"{src}\" from=\"{tgt}\"/>".format(src=edge[0], tgt=edge[1])
                        #print(fmtEdge)

                        if edge[0].lower().startswith("mir") or edge[0].lower().startswith("let"):
                            mirna = edge[0]
                            gene = edge[1]

                        else:

                            mirna = edge[1]
                            gene = edge[0]

                        if (gene, mirna) in existingEdges:
                            continue

                        print(gene, mirna, sep="\t", file=fgout)


    #plt.show()


