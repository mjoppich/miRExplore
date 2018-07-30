from collections import Counter

import networkx
import requests
import json

from utils.cytoscape_grapher import CytoscapeGrapher

class DataBasePlotter:

    @classmethod
    def makePlotForGenes(cls, path, name, gene2name, add=None, cv=False):


        allGenes = set()

        for x in gene2name:

            allGenes.add(x)
            allGenes.add(gene2name[x])

        allGenes = list(allGenes)



        if cv:
            requestData = {
            'disease': [{'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'}],
            'gene': allGenes, 'sentences': "false"}
        else:
            requestData = {'gene': allGenes, 'sentences': "false"}


        if add != None:

            for x in add:
                requestData[x] = add[x]

        print(requestData)

        serverAddress = "https://turingwww.bio.ifi.lmu.de"
        serverPort = None
        serverPath = "yancDB"

        def makeServerAddress(address, port, path):

            ret = address

            if port != None:
                ret += ":" + str(port)

            if path != None:
                ret += "/" + path + "/"

            return ret

        r = requests.post(makeServerAddress(serverAddress, serverPort, serverPath) + "/find_interactions",
                          data=json.dumps(requestData))

        print(r)

        jsonRes = r.json()

        graph = networkx.Graph()

        nodeCounter = Counter()

        print(len(jsonRes['rels']))

        for rel in jsonRes['rels']:

            source = rel['lid']
            target = rel['rid']

            if source.upper() in gene2name:
                source = gene2name[source]

            if target.upper() in gene2name:
                target = gene2name[target]

            graph.add_node(source, {'color': 'red'})
            graph.add_node(target, {'color': 'blue'})

            graph.add_edge(source, target)

            nodeCounter[source] += 1
            nodeCounter[target] += 1

        newNodes = []

        for (node, nodeAttr) in graph.nodes(data=True):

            if node in nodeCounter:
                nodeAttr['size'] = 20 + nodeCounter[node]

            newNodes.append((node, nodeAttr))

        seenNodes = set()

        for (node, nodeAttr) in newNodes:

            if nodeCounter[node] > 0:
                seenNodes.add(node)
            graph.add_node(node, nodeAttr)

        print(set(allGenes).difference(seenNodes))

        mygraph = CytoscapeGrapher.showGraph(graph, location=path,
                                             name=name)



if __name__ == '__main__':

    interactions = {
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
            'miR-532-5p', 'miR-130b-3p', 'miR-222-3p', 'miR144-3p', 'miR-542-3p', 'miR-149-5p', 'miR-330-3p', 'miR-532-3p', 'miR-3470b', 'miR-125b-5p', 'miR-221-3p', 'miR-19b-3p', 'miR-301b-3p',
            'miR-34b-5p', 'miR-125a-3p', 'miR-126-3p', 'miR-16-1-3p', 'miR-882', 'miR-497-5p', 'miR-26a-5p', 'miR-124-3p', 'miR-26b-5p', 'miR-5620-3p', 'mIR-19a-3p', 'miR-130a-3p', 'miR-690',
            'miR-185-5p', 'miR-31-5p', 'miR-340-5p', 'miR-1843-5p', 'miR-466f-3p', 'miR-301a-3p', 'miR-101a-3p', 'miR-210-3p', 'miR-107-3p', 'miR-706', 'miR-23b-3p', 'miR-146a-5p', 'miR-467f',
            'miR-322-5p', 'miR-15a-5p', 'miR-29b-1-5p', 'let-7e-5p', 'miR-23a-3p', 'miR-338-3p', 'miR-103-3p', 'miR-362-3p', 'let-7g-5p', 'miR-155-5p', 'miR-140-5p', 'miR-122-5p', 'miR-22-3p', 'miR-3470a', 'let-7d-5p'
        ]

    }

    genes = [x for x in interactions]
    #genes = ['CCL2', 'CCL3']

    genes2name = {'CXCL12': {'CXCL12'}, 'CXCL13': {'CXCL13'}, 'CXCL14': {'CXCL14'}, 'PPBP': {'PPBP'}, 'XCL1': {'XCL1'},
     'CCL2': {'CCL2'}, 'PF4': {'PF4'}, 'CXCL10': {'CXCL10'}, 'CXCL5': {'CXCL5'}, 'CCL3': {'CCL3'},
     'XCL2': {'XCL1', 'XCL2'}, 'CXCL1': { 'CXCL1'}, 'CCL13': {'CCL13'}, 'CXCL11': {'CXCL11'},
     'CXCL8': {'CXCL8'}, 'CCL14': {'CCL14'}, 'CCL3': {'CCL3', 'CCL3L3'}, 'CCL5': {'CCL5'},
     'CXCL2': { 'CXCL2'}, 'CCL15': {'CCL15'}, 'CXCL3': {'CXCL3'},
     'CCL21': {'CCL21C', 'CCL21A', 'CCL21', 'GM10591', 'GM13304', 'GM21541', 'CCL21B'}, 'CCL17': {'CCL17'},
     'CXCL6': {'CXCL6'}, 'CCL11': {'CCL11'}, 'CCL7': {'CCL7'}, 'CCL4': {'CCL4'}, 'CCL1': {'CCL1'},
     'CXCL16': {'CXCL16'}, 'CCL18': {'CCL18'}, 'CCL19': {'GM2564', 'CCL19'}, 'CXCL9': {'CXCL9'},
     'CCL8': {'CCL8', 'CCL12'}, 'CCL20': {'CCL20'}, 'C5': {'C5', 'HC'}, 'CCL22': {'CCL22'}, 'CCL24': {'CCL24'},
     'CX3CL1': {'CX3CL1'}, 'CCL25': {'CCL25'}, 'CCL28': {'CCL28'}, 'CCL23': {'CCL23'},
     'CCL26': {'CCL26'}, 'CCL27': {'CCL27A', 'CCL27', 'CCL27B', 'GM13306'}, 'CCL6': {'CCL6'}, 'CCL9':{'CCL9'}, 'CCL3': {'CCL3'}}

    gene2name = {}


    allGenes = set()

    for x in genes2name:
        allGenes.add(x)

        for g in genes2name[x]:
            allGenes.add(g)

            gene2name[g] = x


    for x in interactions:

        if x not in allGenes:
            gene2name[x] = x
            allGenes.add(x)

            print("Manual add", x)


    allGenes = list(allGenes)


    print(len(allGenes), allGenes)

    DataBasePlotter.makePlotForGenes('/mnt/c/Users/mjopp/Desktop/yanc_network/', 'all_chemokines', gene2name)






    def subsetGene2Name(newgenes):
        sgene2name = {}

        for x in gene2name:

            if x in newgenes or gene2name[x] in newgenes:
                sgene2name[x] = gene2name[x]

                if x in newgenes:
                    newgenes.remove(x)
                else:
                    newgenes.remove(gene2name[x])

        for gene in newgenes:
            sgene2name[gene] = gene

        return sgene2name



    subPlot = ['CCL2', 'CXCL1', 'CXCL12']
    subPlot += ['CXCR4', 'RGS16', 'HUR', 'ETS1', 'IRAK1', 'TRAF6']
    subPlot += ['SOCS5', 'KLF2', 'KLF4', 'TAK1', 'SIRT1', 'THBS1', 'TGFBR1', 'SMAD2', 'JUN']
    subPlot += ['KPNA4', 'BTRC', 'PPARA']

    sgene2name = subsetGene2Name(subPlot)

    addRestrict = {
        'cells': [{ "group": "cells", "name": "endothelial cell", "termid": "CL:0000115" }]
    }

    DataBasePlotter.makePlotForGenes('/mnt/c/Users/mjopp/Desktop/yanc_network/', 'chemokines_sp1', sgene2name, add=addRestrict, cv=True)

    subPlot = ['CCL2', 'CXCL1', 'CXCL12']
    subPlot += ['CXCR4', 'RGS16', 'HUR', 'ETS1', 'IRAK1', 'TRAF6']
    subPlot += ['SOCS5', 'KLF2', 'KLF4', 'TAK1', 'SIRT1', 'THBS1', 'TGFBR1', 'SMAD2', 'JUN']
    subPlot += ['KPNA4', 'BTRC', 'PPARA']

    sgene2name = subsetGene2Name(subPlot)

    addRestrict = {
        'cells': [{ "group": "cells", "name": "endothelial cell", "termid": "CL:0000115" }]
    }

    DataBasePlotter.makePlotForGenes('/mnt/c/Users/mjopp/Desktop/yanc_network/', 'chemokines_sp1_all_cv', sgene2name, add=None, cv=True)


    subPlot = ['KLF2', 'CHI3L1', 'TLR4', 'TRAF6', 'IRAK1', 'BMPR2', 'AKT1', 'BCL6', 'LPL', 'CCL2']

    sgene2name = subsetGene2Name(subPlot)

    addRestrict = {
        'cells': [{ 'group': "cells", 'name': "monocyte", 'termid': "CL:0000576" }]
    }

    DataBasePlotter.makePlotForGenes('/mnt/c/Users/mjopp/Desktop/yanc_network/', 'chemokines_sp2', sgene2name, add=addRestrict, cv=True)

    subPlot = ['KLF2', 'CHI3L1', 'TLR4', 'TRAF6', 'IRAK1', 'BMPR2', 'AKT1', 'BCL6', 'LPL', 'CCL2']

    sgene2name = subsetGene2Name(subPlot)

    addRestrict = {
        'cells': [{ 'group': "cells", 'name': "monocyte", 'termid': "CL:0000576" }]
    }

    DataBasePlotter.makePlotForGenes('/mnt/c/Users/mjopp/Desktop/yanc_network/', 'chemokines_sp2_all_cv', sgene2name, add=None, cv=True)