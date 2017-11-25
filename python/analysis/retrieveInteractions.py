from database.Neo4JInterface import neo4jInterface
import networkx as nx

from utils.cytoscape_grapher import CytoscapeGrapher


class InteractionRetriever:

    def __init__(self, chemokines = list(), db=neo4jInterface(simulate=False, printQueries=False)):

        self.chemokines = list(set(chemokines))
        self.db = db

        self.all_nodes = dict()
        self.all_relation_ids = set()
        self.all_relations = dict()

        self.simple_connections = set()

        if self.chemokines != None and len(self.chemokines) > 0:
            for chemokine in self.chemokines:
                self._queryDB(chemokine)

    def _addNode(self, node):
        self.all_nodes[node.id] = node

    def _addRelationship(self, rel, test=True):

        if test==True:
            hasStart = rel.start in self.all_nodes
            hasEnd = rel.end in self.all_nodes

            if not (hasStart and hasEnd):
                raise ValueError("relationship nodes not yet included " +  str(rel))

        self.all_relation_ids.add(rel.id)
        self.all_relations[rel.id] = rel


    def _queryDB(self, chemokine):

        print("Query: ", chemokine)
        queryStr = "match p=(g:GENE)-[h]-(s:EVIDENCE)-[r*2]-(t:MIRNA) where g.id='{geneID}' and t.name starts with '{orgPrefix}' return p, LENGTH(p) as length".format(geneID=chemokine, orgPrefix='hsa-')

        res = self.db.runInDatabase(queryStr)

        resVals = [x for x in res]

        print("Found Results: ", len(resVals))

        for hit in resVals:

            path = hit['p']

            pathNodes = path.nodes
            pathRelations = path.relationships

            for node in pathNodes:
                self._addNode(node)

            for rel in pathRelations:
                self._addRelationship(rel)

            simpleRelation = ( path.start.properties['id'], path.end.properties['name'] )
            self.simple_connections.add(simpleRelation)


if __name__ == '__main__':

    #test = InteractionRetriever(chemokines=['CXCL9, CXCL10, CCL22, CCL4, CXCL13, CXCR2, CXCL7, CXCL5'])
    test = InteractionRetriever(chemokines=['CCL17', 'CXCL9', 'CXCL10', 'CCL4', 'CCL22', 'CXCL13', 'CXCR2', 'CXCL7', 'CXCL5'])

    graph = nx.Graph()

    for x in test.simple_connections:
        print(x)
        graph.add_edge(x[0], x[1], attr_dict={'color': "#FF0000"})

    CytoscapeGrapher.showGraph(graph, name='chemokines', title='Chemokines', nodeLabel=lambda x: x, edgeLabel=lambda x: '')
