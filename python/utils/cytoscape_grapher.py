import json
import os
from shutil import copyfile

import networkx
from jinja2 import Environment, PackageLoader, select_autoescape, Template, FileSystemLoader

import matplotlib.pyplot as plt

class CytoscapeGrapher:

    @classmethod
    def export_d3Json(cls, graph, location, name, singletons=False):



        if singletons:
            ng = networkx.Graph()

            for node in graph.nodes():

                nlow = node.lower()
                if nlow.startswith("mir") or nlow.startswith("let"):
                    continue


                metaNode = str(node) + "_meta"

                allEdges = [e for e in graph.edges(node)]

                ng.add_edge(node, metaNode)

                for edge in allEdges:
                    orderedEdge = None
                    if edge[0] == node:
                        orderedEdge = [edge[0], edge[1]]
                    else:
                        orderedEdge = [edge[1], edge[0]]

                    if orderedEdge != None:
                        ng.add_edge(metaNode, orderedEdge[1] + "_" +orderedEdge[0])


            graph = ng

        data1 = networkx.json_graph.node_link_data(graph)

        for x in data1["nodes"]:
            x["name"] = x["id"]

        outpath = os.path.join(location, name)

        if not outpath.endswith(".json"):
            outpath += ".json"

        with open(outpath, 'w') as fp:
            json.dump(data1, fp, indent=4)

    @classmethod
    def jinjaRender(cls,tpl_path, context):
        path, filename = os.path.split(tpl_path)
        return Environment( loader=FileSystemLoader(path or './') ).get_template(filename).render(context)

    @classmethod
    def showGraph(cls, graph, location='/tmp/', name='graph', title=None, nodeLabel=lambda x: x, edgeLabel=lambda x: ''):

        if title == None:
            title=name

        this_dir, this_filename = os.path.split(__file__)
        copyfile(this_dir + '/cytoscape.js', location + "/cytoscape.js")


        graphNodes = []

        for (node, nodeAttr) in graph.nodes(data=True):

            nodeName = nodeLabel(node)
            nodeColor = '#555555'
            nodeId = str(node)

            plotAttr = {'id': nodeId,'name': nodeName, 'color': nodeColor, "size": 20}

            if 'color' in nodeAttr:
                plotAttr['color'] = nodeAttr['color']

                if nodeAttr['color'] in ['g', 'b', 'r']:
                    repMap = {'g': 'green', 'b': 'blue', 'r': 'red'}
                    plotAttr['color'] = repMap[nodeAttr['color']]

            if 'size' in nodeAttr:
                plotAttr['size'] = nodeAttr['size']


            graphNodes.append(plotAttr)

        graphEdges = []
        for edge in graph.edges():

            edata = graph.get_edge_data(edge[0], edge[1])

            eSource = str(edge[0])
            eTarget = str(edge[1])
            eColor = edata['color'] if 'color' in edata else '#9dbaea'
            eLabel = edata['label'] if 'label' in edata else edgeLabel(edge)

            if eColor in ['g', 'b', 'r']:
                repMap = {'g': 'green', 'b': 'blue', 'r': 'red'}
                eColor = repMap[eColor]

            graphEdges.append( {'source': eSource, 'target': eTarget, 'color': eColor, 'label': eLabel} )



        content = {
            'title': title,
            'nodes': graphNodes,
            'edges': graphEdges
        }


        outfile = location + "/" + name + ".html"
        with open(outfile, 'w') as outfile:

            output = cls.jinjaRender(this_dir + '/interactiveCYjsGraph.html', content)

            outfile.write(output)


        return outfile

    @classmethod
    def plotNXGraph(cls, ngraph, title, savePrefix, figidx, figsize=(25,15)):
        plt.figure(figidx, figsize=figsize)

        networkGraph = ngraph

        edges = networkGraph.edges()

        colors = []
        for u, v in edges:
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
            print("change linewidth", title)
            lineWidth = 3

        else:
            print("linewidth unchanged", title)
            lineWidth = 0.75

        pos = networkx.spring_layout(networkGraph)

        networkx.draw(networkGraph, pos, font_size=25, with_labels=False, node_color=nodeColors, edges=edges,
                      edge_color=colors, node_size=1200, width=lineWidth, font_weight='bold', dpi=1000)
        for p in pos:  # raise text positions
            clist = list(pos[p])
            clist[1] = clist[1] + 0.01
            pos[p] = tuple(clist)

        networkx.draw_networkx_labels(networkGraph, pos, font_weight='bold', font_size=25)

        plt.suptitle(title)
        plt.tight_layout()

        if isinstance(savePrefix, list):
            for sp in savePrefix:
                plt.savefig(sp)

        else:
            plt.savefig(savePrefix)



        return figidx+1

    @classmethod
    def exportGraphMLColored(cls, graph, location, name):

        exportGraph = networkx.Graph()

        for edge in graph.edges():
            exportGraph.get_edge_data()
            exportGraph.add_edge(edge[0], edge[1])

        nodeAttribs = {}
        for node in exportGraph.nodes():

            nodename = str(node).lower()

            if nodename.startswith("mir") or nodename.startswith("let"):
                nodeType = "mirna"
            else:
                nodeType = "gene"

            nodeAttribs[node] = {
                "type": nodeType
            }
        networkx.set_node_attributes(exportGraph, nodeAttribs)

        outpath = os.path.join(location, name)

        if not outpath.endswith(".graphml"):
            outpath += ".graphml"


        networkx.write_graphml(exportGraph,  outpath,
                         prettyprint=True)

    @classmethod
    def exportGraphML(cls, graph, location, name):

        exportGraph = networkx.Graph()

        for edge in graph.edges():

            label = None

            edata = graph.get_edge_data(edge[0], edge[1], None)

            if edata != None:
                if 'label' in edata:
                    label = edata['label']
                    print(edge, label)

            exportGraph.add_edge(edge[0], edge[1], label=label)

        nodeAttribs = {}
        for node in exportGraph.nodes():

            nodename = str(node).lower()

            if nodename.startswith("mir") or nodename.startswith("let"):
                nodeType = "mirna"
            else:
                nodeType = "gene"

            nodeAttribs[node] = {
                "type": nodeType
            }
        networkx.set_node_attributes(exportGraph, nodeAttribs)

        outpath = os.path.join(location, name)

        if not outpath.endswith(".graphml"):
            outpath += ".graphml"


        networkx.write_graphml(exportGraph,  outpath,
                         prettyprint=True)

if __name__ == '__main__':


    graph = networkx.Graph()

    graph.add_edge("a", "b")

    mygraph = CytoscapeGrapher.showGraph(graph)

