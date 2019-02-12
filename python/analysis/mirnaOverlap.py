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

    networks = {}

    # endothelial cell activation

    targetMirsECA = [
        'miR-21',
        'miR-92a',
        'miR-217',
        'miR-663',
        'miR-712',
        'miR-7g',
        'let-7g',
        'miR-10a',
        'miR-17-3p',
        'miR-31',
        'miR-124a',
        'miR-125',
        'miR-126',
        'miR-126-5p',
        'miR-143',
        'miR-145',
        'miR-146a',
        'miR-146b',
        'miR-155',
        'miR-181b',
        'miR-221',
        'miR-222']

    networks['targetMirsECA'] = targetMirsECA

    # monocyte

    targetMirsMonocyte = [
        'miR-222',
        'miR-323',
        'miR-503',
        'miR-125b',
        'miR-155',
        'miR-342-5p',
        'miR-17',
        'miR-20a',
        'miR-106a',
        'miR-9',
        'miR-21',
        'miR-124',
        'miR-125a-5p',
        'miR-146a',
        'miR-146b',
        'miR-147',
        'miR-223']

    networks['targetMirsMonocyte'] = targetMirsMonocyte

    # foam cell formation
    targetMirsFCF = [
        'miR-9',
        'miR-125a-5p',
        'miR-146a-5p',
        'miR-155'
    ]

    networks['targetMirsFCF'] = targetMirsFCF


    # Angiogenesis
    targetMirsAngio = [
        'let-7f',
        'miR-7f',
        'miR-23',
        'miR-24',
        'miR-27',
        'miR-126',
        'miR-130a',
        'miR-132',
        'miR-150',
        'miR-210',
        'miR-218',
        'miR-378',
        'miR-15b',
        'miR-16',
        'miR-20a',
        'miR-21',
        'miR-26a',
        'miR-17',
        'miR-92',
        'miR-100',
        'miR-200',
        'miR-221',
        'miR-222',
        'miR-223']

    networks['targetMirsAngio'] = targetMirsAngio


    # Vascular remodeling
    targetMirsVasRemod = [
        'miR-21',
        'miR-155',
        'miR-222',
        'miR-126',
        'miR-143',
        'miR-145']

    networks['targetMirsVasRemod'] = targetMirsVasRemod


    # T - cell differentiation and activation

    targetMirsTCell = [
        'miR-17',
        'miR-92',
        'miR-146a',
        'miR-155',
        'miR-182',
        'miR-326',
        'miR-125b',
        'miR-181a']

    networks['targetMirsTCell'] = targetMirsTCell


    # Cholestrol efflux

    targetMirsCholEfflux = [
        'miR-10b',
        'miR-26',
        'miR-27',
        'miR-33a',
        'miR-106b',
        'miR-144',
        'miR-145',
        'miR-155',
        'miR-302a',
        'miR-758',
        'miR-223',
        'miR-378']


    networks['targetMirsCholEfflux'] = targetMirsCholEfflux


    # SMC proliferation / migration
    targetMirsSMCProlif = [
        'miR-24',
        'miR-26a',
        'miR-31',
        'miR-146a',
        'miR-155',
        'miR-208',
        'miR-221',
        'miR-222',
        'miR-7d',
        'let-7d',
        'miR-1',
        'miR-10a',
        'miR-21',
        'miR-29',
        'miR-100',
        'miR-132',
        'miR-133',
        'miR-143',
        'miR-145',
        'miR-195',
        'miR-204',
        'miR-424',
        'miR-638',
        'miR-663']

    networks['targetMirsSMCProlif'] = targetMirsSMCProlif

    summaryDF = DataFrame()
    summaryDF.addColumns(["Network", "Accepted miRNAs", "Missing miRNAs"])

    networkGraphs = {}
    makeStory = [

    ]

    allNetworks = [x for x in networks]
    print(allNetworks)

    #exit()

    ignoreNetworks = []

    networkRestrictions = {
        'targetMirsECA': {
            "cells": [
                {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
            ]
        },
        'targetMirsMonocyte': {
            "cells": [
                {"group": "cells", "name": "monocyte", "termid": "META:148"},
                {"group": "cells", "name": "macrophage", "termid": "META:99"}
            ]
        },
        'targetMirsFCF': {
            "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"}]
        },
        'targetMirsAngio': {
            #"cells": [{"group": "cells", "name": "blood vessel", "termid": "UBERON:0001981"}, {"group": "cells", "name": "blood vessel elastic tissue", "termid": "UBERON:0003614"} , {"group": "cells", "name": "arterial blood vessel", "termid": "UBERON:0003509"}],
            "go": [{"group": "go", "name": "angiogenesis", "termid": "GO:0001525"}]
        },
        'targetMirsVasRemod': {
            "cells": [{"group": "cells", "name": "vascular disease", "termid": "DOID:178"}]
        }
        ,'targetMirsTCell': {
            "cells": [{"group": "cells", "name": "T cell", "termid": "META:44"}]
        },
        'targetMirsCholEfflux': {
            "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"}]
        },
        'targetMirsSMCProlif': {
            "cells": [{"group": "cells", "name": "smooth muscle cell", "termid": "META:83"}]
        }

    }

    networkToTitle = {
        "targetMirsECA": "\makecell[l]{Endothelial cell activation\\\\and inflammation}",
        "targetMirsMonocyte": "\makecell[l]{Monocyte differentiation\\\\Macrophage activation}",
        "targetMirsFCF": "\makecell[l]{Foam cell formation}",
        "targetMirsAngio": "\makecell[l]{Angiogenesis}",
        "targetMirsVasRemod": "\makecell[l]{Vascular remodeling}",
        "targetMirsTCell": "\makecell[l]{T-cell differentiation\\\\and activation}",
        "targetMirsCholEfflux": "\makecell[l]{Cholesterol efflux}",
        "targetMirsSMCProlif": "\makecell[l]{SMC proliferation\\\\SMC migration}"

    }

    figidx = 0

    for network in networks:
        figidx+= 1

        networkGraph = nx.Graph()

        if network in ignoreNetworks:
            continue


        interactions = defaultdict(set)

        acceptedInteractions = defaultdict(set)
        typeByGene = defaultdict(lambda: Counter())
        elemsByGene = defaultdict(lambda: defaultdict(set))

        allMirna = set(networks[network])

        miStr2mirna = {}
        allTargetMirna = []
        mirnaObj2str = {}

        for x in allMirna:
            try:

                oMirna = miRNA(x)

                allTargetMirna.append( oMirna )
                miStr2mirna[x] = oMirna
                mirnaObj2str[oMirna] = x

            except:
                pass

        #allMirna = set([str(x) for x in allTargetMirna])
        requestData = None


        if network in networkRestrictions:
            requestData = networkRestrictions[network]
        else:
            requestData = {
                            'sentences': "false",
                            }

        requestData['sentences'] = "false"


        requestData["mirna"]= list(allMirna)
        print(allMirna)
        requestData['disease'] = [
                    {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
                    {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
                ]


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

                origEdge = tuple(oEdge)
                edgeStatus = None
                oEdge = list(oEdge)

                for miObj in mirnaObj2str:

                    if miObj.accept(oEdge[1]):
                        oEdge[1] = mirnaObj2str[miObj]
                        break

                allGeneMirna = interactions[oEdge[0]]
                miAccepted = False

                allAcceptedStr = set()

                if "758" in oEdge[1] or "758" in oEdge[0]:
                    print(oEdge)

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



        elemsByMirna = defaultdict(set)

        for gene in elemsByGene:

            for expMirna in allMirna:

                for category in elemsByGene[gene]:

                    if category == 'missing':
                        continue

                    for foundMirna in elemsByGene[gene][category]:

                        elemsByMirna[foundMirna].add(gene)


        foundMirnas = set([x for x in elemsByMirna])

        print(network)
        print("Found Mirnas", len(foundMirnas), list(foundMirnas))
        print("Expected Mirnas", len(allMirna), list(allMirna))
        print("Intersected Mirnas", len(foundMirnas.intersection(allMirna)), list(foundMirnas.intersection(allMirna)))
        print("Missing Mirnas", len(allMirna.difference(foundMirnas)), allMirna.difference(foundMirnas))

        rowDict = {}
        rowDict['Network'] = networkToTitle[network]
        rowDict['Accepted miRNAs'] = "\makecell[l]{" +"\\\\".join(sorted(foundMirnas.intersection(allMirna))) + "}"
        rowDict['Missing miRNAs'] = "\makecell[l]{"  + "\\\\".join(sorted(allMirna.difference(foundMirnas))) + "}"

        newRow = DataRow.fromDict(rowDict)

        #["Network", "Accepted miRNAs", "Missing miRNAs"]
        summaryDF.addRow( newRow )

        if False:

            print(network)
            for gene in sorted([x for x in typeByGene]):
                print(gene, typeByGene[gene], elemsByGene[gene]['missing'])

            print()
            print()

            print(network)
            for gene in sorted([x for x in typeByGene]):
                print("Gene:", gene, "Status: ", ", ".join([": ".join([x, str(typeByGene[gene][x])]) for x in typeByGene[gene]]), "Missing miRNAs: "+",".join(elemsByGene[gene]['missing']))

            print()
            print()
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

    print(summaryDF._makeLatex())



