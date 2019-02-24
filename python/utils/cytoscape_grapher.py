import os
from shutil import copyfile

import networkx
from jinja2 import Environment, PackageLoader, select_autoescape, Template, FileSystemLoader

import matplotlib.pyplot as plt

class CytoscapeGrapher:


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

if __name__ == '__main__':


    graph = networkx.Graph()

    graph.add_edge("a", "b")

    mygraph = CytoscapeGrapher.showGraph(graph)

