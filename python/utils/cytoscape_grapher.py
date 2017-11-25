import os
from shutil import copyfile

import networkx
from jinja2 import Environment, PackageLoader, select_autoescape, Template, FileSystemLoader


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

        for node in graph.nodes_iter():

            nodeName = nodeLabel(node)
            nodeColor = '#555555'
            nodeId = str(node)

            graphNodes.append({'id': nodeId,'name': nodeName,'color': nodeColor})

        graphEdges = []
        for edge in graph.edges():

            edata = graph.get_edge_data(edge[0], edge[1])

            eSource = str(edge[0])
            eTarget = str(edge[1])
            eColor = edata['color'] if 'color' in edata else '#9dbaea'
            eLabel = edata['label'] if 'label' in edata else edgeLabel(edge)

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

if __name__ == '__main__':


    graph = networkx.Graph()

    graph.add_edge("a", "b")

    mygraph = CytoscapeGrapher.showGraph(graph)
