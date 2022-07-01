import sys, os
from collections import Counter, OrderedDict, defaultdict

import networkx
from natsort import natsorted

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))


from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher
from utils.tmutils import normalize_gene_names

import pandas as pd
from pandas.tools.plotting import parallel_coordinates
import matplotlib.pyplot as plt

figidx = 0



normGeneSymbols = normalize_gene_names(path="/mnt/d/owncloud/data/miRExplore/obodir/" + "/hgnc_no_withdrawn.syn")


cells2restrict = OrderedDict([(

    "endothelial_cells", {
        "sentences": "false",
        "cells": [
                {"group": "cells", "name": "endothelial cell",  "termid": "META:52"},
                {"group": "cells", "name": "HUVEC-C",  "termid": "META:02689"}
              ]
        }),
    ("monocytes", {
        "sentences": "false",
        "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"},
            {"group": "cells", "name": "monocyte", "termid": "META:148"},
            {"group": "cells", "name": "macrophage", "termid": "META:99"},
            {"group": "cells", "name": "T cell", "termid": "META:44"}
            ]
    }),
    ("foam_cells", {
        "sentences": "false",
        "cells": [
                  {"group": "cells", "name": "foam cell", "termid": "CL:0000891"},
                  ],
    }),
    ("smooth_muscle_Cell", {
        "sentences": "false",
        "cells": [
            {"group": "cells", "name": "smooth muscle cell", "termid": "META:83"},
        ]
    }),
    ("platelet", {
            'sentences': "false",
            "cells": [
                { "group": "cells", "name": "thromboblast", "termid": "CL:0000828" }
            ],
        }),
    ])

network2nicename = {

}


def acceptEvidence(ev):

    return True

    if ev['data_source'] == 'miRTarBase':

        if "Weak" in ev['functional_type']:
            return False

    return True




cbn2graph = {}

cells2edges = {}

for x in cells2restrict:
    cells2edges[x] = defaultdict(set)


mirnaNodes = set()

stageMir2Cells = defaultdict(lambda: defaultdict(set))
stageMir2CellEvidence = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

for cbnNW in cells2edges:

    allNWGenes = set()
    newEdges = set()

    for edge in cells2edges[cbnNW]:

        src = edge[0].upper()
        tgt = edge[1].upper()

        if src in normGeneSymbols:
            src = normGeneSymbols[src]

        if tgt in normGeneSymbols:
            tgt = normGeneSymbols[tgt]

        cells2edges[cbnNW] = newEdges

    #print(cbnNW, allNWGenes)

    requestData = cells2restrict[cbnNW]

    requestData['disease'] = [
        {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
        {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
    ]

    requestData['gene'] = []

    print(cbnNW)
    print(requestData)

    graph, nodeCounts, evcount, json = DataBasePlotter.fetchGenes(requestData, gene2name=None, minPMIDEvCount=0, minTgtCount=0, acceptEv=acceptEvidence)

    newNodes = []

    for (node, nodeAttr) in graph.nodes(data=True):

        if node in nodeCounts:
            nodeAttr['size'] = 20 + nodeCounts[node]

        if nodeAttr['color'] == 'blue':
            mirnaNodes.add(node)

        newNodes.append((node, nodeAttr))

    removeNodes = [node for node, degree in graph.degree() if degree >= 2 and node in mirnaNodes]
    #graph.remove_nodes_from(removeNodes)


    removeEdges = set()

    singleEvidenceCounter = Counter()
    edgeEvidenceCounter = Counter()

    for (edgeSrc, edgeTgt, edgeAttr) in graph.edges(data=True):

        edge = (edgeSrc, edgeTgt)
        edgeR = (edgeTgt, edgeSrc)

        evC = None

        if edge in evcount:
            evC = evcount[edge]
        elif edgeR in evcount:
            evC = evcount[edgeR]

        if evC != None:

            if len(evC) == 1:

                for x in evC:
                    singleEvidenceCounter[x] += 1

                edgeAttr['color'] = '#00ff00'

            elabel = ""
            for x in sorted([x for x in evC]):
                elabel += x[0]

            edgeAttr["label"] = elabel

            edgeEvidenceCounter[elabel] += 1

            miRName = None

            if edge[0].startswith("miR") or edge[0].startswith("let-"):
                miRName = edge[0]
            elif edge[1].startswith("miR") or edge[1].startswith("let-"):
                miRName = edge[1]


            if "celldata" in edgeAttr and miRName != None :

                for x in edgeAttr["celldata"]:
                    stageMir2Cells[cbnNW][miRName].add(x)

                for x in edgeAttr['cellEvidence']:
                    for pmid in edgeAttr['cellEvidence'][x]:
                        stageMir2CellEvidence[cbnNW][miRName][x].add(pmid)

            elif miRName == None:
                print("INVALID MIRNA:", edge, file=sys.stderr)

        else:
            print(edge, "not in evcount")

    print(singleEvidenceCounter)
    print(edgeEvidenceCounter)

    for edge in removeEdges:
        graph.remove_edge(edge[0], edge[1])

    graph.remove_nodes_from(networkx.isolates(graph))

    for edge in cells2edges[cbnNW]:
        if edge[0] == edge[1]:
            continue

        graph.add_edge(edge[0], edge[1], {'color': '#aa0000'})

    cbn2graph[cbnNW] = graph


### which cell types are most common per stage?
print()
print("cells per stage")
print()
print()
print()

importantCellTypesPerMirna = defaultdict(lambda: Counter())

for stage in stageMir2Cells:

    cellCounter = Counter()

    stageMirnaCellPairs = Counter()
    mirnaCellPairs = defaultdict(lambda: Counter())
    stageCellCount = Counter()

    for mirna in stageMir2Cells[stage]:


        for cell in stageMir2Cells[stage][mirna]:
            cellCounter[cell] += 1

        mirnaStageCells = [x for x in stageMir2Cells[stage][mirna]]

        for i in range(0, len(mirnaStageCells)):
            for j in range(i+1, len(mirnaStageCells)):

                if mirnaStageCells[i][0].startswith("CVCL") or mirnaStageCells[j][0].startswith("CVCL"):
                    continue

                cellpair = tuple(sorted((mirnaStageCells[i],mirnaStageCells[j])))

                mirnaCellPairs[mirna][cellpair] += 1
                stageMirnaCellPairs[cellpair] += 1

                stageCellCount[cellpair[0]] += 1
                stageCellCount[cellpair[1]] += 1

                importantCellTypesPerMirna[mirna][cellpair] += 1


    mostCommonCellPairs = []
    for (mpair, count) in stageMirnaCellPairs.most_common(): #20
        mostCommonCellPairs.append(mpair)


    edge2support = defaultdict(set)

    for mirna in mirnaCellPairs:
        for cellpair in mirnaCellPairs[mirna]:
            if cellpair in mostCommonCellPairs:

                edge2support[cellpair].add(mirna)

                print(stage, mirna, cellpair[0], cellpair[1], mirnaCellPairs[mirna][cellpair], stageMirnaCellPairs[cellpair], stageMir2CellEvidence[stage][mirna].get(cellpair[0]),stageMir2CellEvidence[stage][mirna].get(cellpair[1]) )

    cellgraph = networkx.Graph()

    allnodes = set()
    for edge in edge2support:
        allnodes.add(edge[0])
        allnodes.add(edge[1])

    for node in allnodes:
        cellgraph.add_node(node[1] + " ("+node[0]+")", size=20 + stageCellCount[node])


    cellCommunicatorDF = DataFrame()
    cellCommunicatorDF.addColumns(["miRNA", "cells"])

    mirna2cells = defaultdict(set)

    for edge in edge2support:
        cellgraph.add_edge(
            edge[0][1] + " (" + edge[0][0] + ")",
            edge[1][1] + " (" + edge[1][0] + ")",
            label=", ".join(edge2support.get(edge, [])))

        mirnas = edge2support.get(edge, [])

        for mirna in mirnas:
            mirna2cells[mirna].add(edge[0][1] + " (" + edge[0][0] + ")")
            mirna2cells[mirna].add(edge[1][1] + " (" + edge[1][0] + ")")


    cells2mirnas = defaultdict(set)
    for mirna in mirna2cells:
        cells = tuple(sorted(mirna2cells[mirna]))

        cells2mirnas[cells].add(mirna)


    for cells in cells2mirnas:

        ncells = []
        for cell in cells:
            ncells.append( ""+cell+"}" )

        rowdict = {
            'miRNA': "\\makecell[l]{" + "\\\\".join(sorted(cells2mirnas[cells])) + "}",
            'cells': "\\makecell[l]{" + "\\\\".join(cells) + "}"
        }
        cellCommunicatorDF.addRow(DataRow.fromDict(rowdict))


    cellCommunicatorDF.export("/mnt/d/yanc_network/stage_cells_cliques"+stage+".latex", ExportTYPE.LATEX)
    cellCommunicatorDF.export("/mnt/d/yanc_network/stage_cells_cliques"+stage+".tsv", ExportTYPE.TSV)


    CytoscapeGrapher.showGraph(cellgraph, '/mnt/d/yanc_network/', name="stage_cells_" + stage)
    figidx = CytoscapeGrapher.plotNXGraph(cellgraph, stage, ["/mnt/d/yanc_network/stage_cells_" + stage.replace(" ", "_") + ".png", "/mnt/d/yanc_network/stage_cells_" + stage.replace(" ", "_") + ".pdf"], figidx)

    print()
    print()
    print()
    print()
    print()

    plotKeys = []
    plotValues = []
    for (cell, count) in cellCounter.most_common(20):
        print(stage, cell[0], cell[1], count)
        plotKeys.append(cell[1] + " ("+cell[0]+")")
        plotValues.append(count)


    plt.figure(figidx)
    figidx+= 1
    plt.bar(plotKeys, plotValues)
    plt.xticks(rotation=90)
    plt.tight_layout()

    plt.savefig("/mnt/d/yanc_network/cells_per_stage_" + stage + ".png",bbox_inches='tight')
    plt.savefig("/mnt/d/yanc_network/cells_per_stage_" + stage + ".pdf",bbox_inches='tight')

    plt.show()
    print()
    print()
    print()
    print()
    print()

for mirna in importantCellTypesPerMirna:

    for (cp, count) in importantCellTypesPerMirna[mirna].most_common(10):
        print(mirna, cp[0], cp[1], count)

print()
print()
print()
print()
print()

### which miRNAs play a role in all networks?

### are there miRNAs which orchestrate the transition from one state to the other? which ones are unique to each stage? only for miRNAs connected to genes present in both stages

makeStory = [
    ['CV-IPN-Endothelial_cell_activation_1','CV-IPN-Endothelial_cell-monocyte_interaction_1','CV-IPN-Foam_cell_formation_1','CV-IPN-Smooth_muscle_cell_activation_1','CV-IPN-Platelet_activation_1',
     'CV-IPN-Plaque_destabilization_1'
     ]
]

storyTransitions = []


for storyline in makeStory:
    for i in range(0, len(storyline)-1):
        storyTransitions.append( (storyline[i], storyline[i+1]) )


def getMirsForGene(gene, graph):

    alledges = graph.edges(gene)

    targetMirs = set()

    for u,v in alledges:

        if u == gene:
            if v.startswith("miR"):
               targetMirs.add(v)

        elif v == gene:
            if u.startswith("miR"):
                targetMirs.add(u)

    return targetMirs

outfile = open("/mnt/d/yanc_network/" + "mirs_per_stage.tsv", 'w')

print("Condition1", "Condition2", "Gene", "miRNAs Only Before", "miRNAs Only After", "Common Mirs", sep="\t")
print("Condition1", "Condition2", "Gene", "miRNAs Only Before", "miRNAs Only After", "Common Mirs", sep="\t", file=outfile)
for storyTransition in storyTransitions:

    graphFrom = cbn2graph[storyTransition[0]]
    graphTo = cbn2graph[storyTransition[1]]

    fromGenes = set()
    toGenes = set()

    for (node, nodeAttr) in graphFrom.nodes(data=True):
        if not node.startswith("miR"):
            fromGenes.add(node)

    for (node, nodeAttr) in graphTo.nodes(data=True):
        if not node.startswith("miR"):
            toGenes.add(node)

    commonGenes = fromGenes.intersection(toGenes)

    genesWithDiff = {}

    for cgene in commonGenes:

        fromMirs = getMirsForGene(cgene, graphFrom)
        toMirs = getMirsForGene(cgene, graphTo)

        diffMirs = set()
        toDiffMirs = set()
        fromDiffMirs = set()
        commonMirs = set()
        for x in fromMirs.union(toMirs):
            if not x in fromMirs and x in toMirs:
                toDiffMirs.add(x)
            elif x in fromMirs and not x in toMirs:
                fromDiffMirs.add(x)
            elif x in fromMirs and x in toMirs:
                commonMirs.add(x)


        if len(toDiffMirs) > 0 or len(fromDiffMirs) > 0:
            genesWithDiff[cgene] = (fromDiffMirs, toDiffMirs, commonMirs)


    for cgene in genesWithDiff:
        mirset = genesWithDiff[cgene]
        print(storyTransition[0], storyTransition[1], cgene, ",".join(sorted(mirset[0])), ",".join(sorted(mirset[1])),",".join(sorted(mirset[2])), sep="\t")
        print(storyTransition[0], storyTransition[1], cgene, ",".join(sorted(mirset[0])), ",".join(sorted(mirset[1])),",".join(sorted(mirset[2])), sep="\t", file=outfile)

outfile.close()


latexStageMirsFile = open("/mnt/d/yanc_network/cbn_stages_most_reg.txt", 'w')

mirSByCBN = defaultdict(set)

mir2networks = defaultdict(set)
interactionsPerNetworkAndMirna = defaultdict(lambda: dict())

overallMostRegulatingMIRs = Counter()

for cbnNW in cells2restrict:

    graph = cbn2graph[cbnNW]

    for (node, nodeAttr) in graph.nodes(data=True):

        if 'color' in nodeAttr and nodeAttr['color'] == 'blue':

            otherGraphHasNode = False

            for oNW in cbn2graph:
                if oNW == cbnNW:
                    continue

                ograph = cbn2graph[oNW]

                if ograph.has_node(node):
                    otherGraphHasNode = True
                    break

            if not otherGraphHasNode:
                nodeAttr['color'] = '#d942f4'


    intNodes = [(node,degree) for node, degree in graph.degree() if node in mirnaNodes]
    intNodes = sorted(intNodes, key=lambda x: x[1], reverse=True)


    for node in intNodes:
        if node[0].startswith("miR"):
            mirSByCBN[cbnNW].add(node[0])
            mir2networks[node[0]].add(cbnNW)

            overallMostRegulatingMIRs[node[0]] += node[1]

            interactionsPerNetworkAndMirna[cbnNW][node[0]] = node[1]


    print(cbnNW)
    cbnNWMirs = []
    for i in range(0, min(10, len(intNodes))):
        print(i, intNodes[i])
        cbnNWMirs.append(  intNodes[i][0] + " ("+str(intNodes[i][1])+")"  )

    latexStageMirsFile.write(cbnNW + "&" + ", ".join(cbnNWMirs) + "\n")

    print()



    print()

    mygraph = CytoscapeGrapher.showGraph(graph, location='/mnt/d/yanc_network/', name=cbnNW)

latexStageMirsFile.close()

print("Overall most regulating")
print(overallMostRegulatingMIRs.most_common(10))
print(", ".join([x[0] + " (" + str(x[1]) + ")" for x in overallMostRegulatingMIRs.most_common(10)]))
print()
print()


mirna2nwcount = Counter()
for mirna in mir2networks:
    mirna2nwcount[mirna] = len(mir2networks[mirna])

print("most common mirnas")
mcMirnas = set()
for x in mirna2nwcount.most_common(10):
    print(x)
    mcMirnas.add(x[0])

print(mcMirnas)

mirna2stagecount = defaultdict(list)
stages = [x for x in cells2restrict]
for mirna in natsorted(mcMirnas):

    mirnarow = [mirna]
    for stage in stages:
        icount = interactionsPerNetworkAndMirna[stage].get(mirna, 0)
        mirnarow.append(icount)

        print(stage, mirna, icount, sep="\t")


    mirna2stagecount[mirna] = mirnarow




stageDF = pd.DataFrame.from_dict(mirna2stagecount, orient='index', columns=["miRNA"] + [network2nicename[x] for x in stages])

ax = parallel_coordinates(stageDF, 'miRNA', colormap=plt.get_cmap("Set2"))
plt.xticks(rotation=90)

lgd = plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=5)

plt.tight_layout()

plt.savefig("/mnt/d/yanc_network/mirna_stage_parallel.png",bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.savefig("/mnt/d/yanc_network/mirna_stage_parallel.pdf",bbox_extra_artists=(lgd,), bbox_inches='tight')

plt.show()


print()
print()


allSet = None

for nw in mirSByCBN:

    if allSet == None:
        allSet = mirSByCBN[nw]
    else:
        allSet = allSet.intersection(mirSByCBN[nw])

print("miRNAs in all nws")
print(allSet)


import matplotlib.pyplot as plt

for stages in makeStory:

    mergedGraph = cbn2graph[stages[0]]

    for i in range(1, len(stages)):
        mergedGraph = networkx.compose(mergedGraph, cbn2graph[stages[i]])


    pos = networkx.spring_layout(mergedGraph)

    for stage in stages:
        plt.figure(figidx, figsize=(50,30))
        figidx += 1

        networkGraph = cbn2graph[stage]

        edges = networkGraph.edges()

        colors = []
        for u,v in edges:
            elem = networkGraph[u][v]

            if not 'color' in elem:
                elem['color'] = "#0000FF"

            colors.append(elem['color'])


        nodes = networkGraph.nodes()
        nodeColors = []

        mirNodes = 0
        geneNodes = 0

        for x in nodes:
            if any([x.lower().startswith(y) for y in ['mir', 'let']]):
                nodeColors.append('blue')
                mirNodes += 1
            else:
                nodeColors.append('green')
                geneNodes += 1


        lineWidth = 0.5

        if len(nodes) < 200 and len(edges) < 200:
            print("change linewidth", stage)
            lineWidth = 3

        else:
            print("linewidth unchanged", stage)
            lineWidth=0.75

        networkx.draw(networkGraph, pos, font_size=25, with_labels=False, node_color=nodeColors, edges=edges, edge_color=colors, node_size=1200, width=lineWidth, font_weight='bold', dpi=1000)
        for p in pos:  # raise text positions
            clist = list(pos[p])
            clist[1] = clist[1] + 0.02
            pos[p] = tuple(clist)

        networkx.draw_networkx_labels(networkGraph, pos, font_weight='bold', font_size=25)

        plt.suptitle("{title}: {mirc} miRs, {genec} genes, {intc} interactions".format(title=stage, mirc=mirNodes, genec=geneNodes, intc=len(colors)))

        plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".png")
        plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".pdf")



