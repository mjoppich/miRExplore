import pickle
import tempfile

from collections import defaultdict, Counter
import networkx as nx

from networkx.drawing.nx_agraph import graphviz_layout, pygraphviz_layout
import os, sys

from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from synonymes.mirnaID import miRNA, miRNAPART
from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher

import matplotlib.pyplot as plt

if __name__ == '__main__':

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
        'CXCR4': ['miR-381-3p', 'miR-21-3p', 'miR-467a-5p', 'miR-467h', 'miR-218-5p', 'miR-1a-3p', 'miR-181d-5p', 'miR-206-3p', 'miR-181b-5p', 'miR-9-5p', 'miR-132-3p', 'miR-25-3p', 'miR-467d-5p', 'miR-669k-3p', 'miR-146b-5p', 'miR-467b-5p', 'miR-467e-5p', 'miR-467f', 'miR-146a-5p'],
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
        'ETS1': ['miR-155', 'miR-221', 'miR-222']
    }

    ccl2_macrophage = {
        'CCL2': ['miR-124a', 'miR-150'],
        'CHI3L1': ['miR-24'],
        'TLR4': ['miR-146'],
        'TRAF6': ['miR-146'],
        'IRAK1': ['miR-146'],
        'AKT1': ['miR-342-5p'],
        'BCL6': ['miR-155'],
        'LPL': ['miR-467a']
    }

    normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")

    andreouInts1= {'ICAM1': ['miR-223', 'miR-17-3p'], 'SELE': ['miR-31'], 'VCAM1': ['miR-126', 'miR-126-5p'], 'DLK1': ['miR-126', 'miR-126-5p'], 'JAMA': ['miR-143', 'miR-145'], 'SIRT1': ['miR-217'], 'NOS3': ['miR-155'], 'ETS1': ['miR-155', 'miR-126', 'miR-126-5p'], 'PPARA': ['miR-21'], 'TIMP3': ['miR-712'], 'MAP3K7': ['miR-10'], 'KPNA4': ['miR-181b'], 'TRAF6': ['miR-146a', 'miR-146b'], 'IRAK1': ['miR-146a', 'miR-146b'], 'IRAK2': ['miR-146a', 'miR-146b'], 'KLF2': ['miR-143', 'miR-145', 'miR-126', 'miR-126-5p', 'miR-92a'], 'CXCL12': ['miR-126', 'miR-126-5p'], 'KLF4': ['miR-92a', 'miR-663'], 'SOCS5': ['miR-92a']}
    andreouInts2= {'IFNG': ['miR-29', 'miR-155'], 'CD40': ['miR-146a', 'miR-181a'], 'CD40L': ['miR-146a', 'miR-181a'], 'MMP2': ['miR-21', 'miR-21', 'miR-24', 'miR-29'], 'MMP3': ['miR-21', 'miR-24', 'miR-29'], 'MMP8': ['miR-21', 'miR-24', 'miR-29'], 'MMP9': ['miR-21', 'miR-24', 'miR-29'], 'MMP12': ['miR-21', 'miR-24', 'miR-29'], 'MMP13': ['miR-21', 'miR-24', 'miR-29'], 'MMP14': ['miR-21', 'miR-24', 'miR-29']}
    andreouInts3= {'MTTP': ['miR-30c'], 'ABCA1': ['miR-33', 'miR-302a'], 'ABCG1': ['miR-33'], 'CPT1A': ['miR-33'], 'KLF2': ['miR-92a'], 'KLF4': ['miR-92a', 'miR-145'], 'SOCS5': ['miR-92a'], 'DLK1': ['miR-126-5p'], 'RGS16': ['miR-126-3p'], 'BCL6': ['miR-155'], 'SOCS1': ['miR-155'], 'MAP3K10': ['miR-155'], 'KPNA4': ['miR-181b'], 'AKT1': ['miR-342-5p'], 'TIMP3': ['miR-712']}


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

    networks['macrophages'] = ccl2_macrophage
    networks['macrophages_cv'] = ccl2_macrophage
    networks['macrophages_all'] = ccl2_macrophage

    networkGraphs = {}
    makeStory = [
        ('macrophages', 'macrophages_cv', 'macrophages_all'),
        ('inflammatory_ec','inflammatory_ec_cv'),
        ('large_chemokines', 'large_chemokines_cv'),
        ('andreou_fig1', 'andreou_fig1_cv'),
        ('andreou_fig2', 'andreou_fig2_cv'),
        ('andreou_table1', 'andreou_table1_cv'),
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
        'large_chemokines': {
            'sentences': "false",
        },
        'large_chemokines_cv': {
            'sentences': "false",
            "disease": [
                {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
            ]        },

        'andreou_fig1': {
            'sentences': "false",
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
                ]
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
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },
        'andreou_fig2_athero':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },
        'andreou_table1_athero':
            {
                'sentences': "false",
                "disease": [
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]
            },
    }

    figidx = 0

    for network in networks:
        figidx+= 1

        networkGraph = nx.Graph()

        if network in ignoreNetworks:
            continue

        interactions = networks[network]
        acceptedInteractions = defaultdict(set)
        typeByGene = defaultdict(lambda: Counter())
        elemsByGene = defaultdict(lambda: defaultdict(set))

        allMirna = set()

        miStr2mirna = {}

        for gene in interactions:

            for mirna in interactions[gene]:

                allMirna.add(mirna)

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


        graph, nodeCounter, edge2datasourceCount, jsonRes = DataBasePlotter.fetchGenes(requestData)

        print(len(jsonRes['rels']))

        htmlDF = DataFrame()
        htmlDF.addColumns(['gene rel', 'gene', 'miRNA Group', 'miRNA', 'Original Network', 'PubMed', 'MIRECORD', 'MIRTARBASE', 'DIANA'])

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

                edgeStatus = None

                allGeneMirna = interactions[oEdge[0]]
                miAccepted = False

                allAcceptedStr = set()

                for strMirna in allGeneMirna:
                    miObj = miStr2mirna.get(strMirna, None)

                    if miObj == None:
                        continue

                    miObjAccepts = miObj.accept(oEdge[1])
                    miAccepted = miAccepted or miObjAccepts


                    if miObjAccepts:
                        acceptedInteractions[oEdge[0]].add(strMirna)
                        edgeStatus = "accepted"

                        allAcceptedStr.add(strMirna)

                        #print(oEdge[0], oEdge[1], strMirna)

                if not miAccepted:
                    edgeStatus = "additional"

                networkGraph.add_edge(oEdge[0], oEdge[1], color= 'g' if edgeStatus == "accepted" else "b")

                typeByGene[oEdge[0]][edgeStatus] += 1
                elemsByGene[oEdge[0]][edgeStatus].add(oEdge[1])

                objMirna = miRNA(oEdge[1])

                pmidEvs = set()
                mirtarbaseEvs = set()
                mirecordsEvs = set()
                dianaEvs = set()

                for ev in rel['evidences']:

                    if ev['data_source'] == "DIANA":
                        dianaEvs.add( (ev['method'], ev['direction']) )

                    elif ev['data_source'] == "miRTarBase":
                        mirtarbaseEvs.add(
                            (ev['data_id'], ",".join(ev['exp_support']), ev['functional_type'], ev['docid'])
                        )
                    elif ev['data_source'] == "pmid":
                        pmidEvs.add((ev['docid'],))
                    elif ev['data_source'] == "mirecords":
                        mirecordsEvs.add((ev['docid']))
                    else:
                        print("Unhandled data source", ev['data_source'])


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
                    'DIANA': dianaStr
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)

        for gene in interactions:

            for mirna in interactions[gene]:

                edgeWasFound = mirna in acceptedInteractions[gene]

                if edgeWasFound:
                    continue

                edgeStatus = "missing"

                networkGraph.add_edge(gene, mirna, color='r')


                typeByGene[gene][edgeStatus] += 1
                elemsByGene[gene][edgeStatus].add(mirna)

                objMirna = miRNA(mirna)
                dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
                    geneCap=gene.upper(), geneLow=gene.upper())

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
                    'DIANA': ""
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)



        print(network)
        for gene in sorted([x for x in typeByGene]):
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

    #plt.show()


