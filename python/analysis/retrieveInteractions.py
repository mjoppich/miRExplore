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
    test = InteractionRetriever(chemokines=['CXCR2','CCL9', 'CXCL5', 'CXCL1', 'CXCL13', 'CXCL7', 'CXCL12', 'CXCL14', 'CCR7', 'CXCL9', 'CCL2', 'CCR7', 'CCL7', 'CCL3', 'CXCL10', 'CCL22', 'CCL22', 'CCL4', 'CXCR4', 'CX3CL1'])


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
        'CCL22': ['miR-34a-5p'],
        'CXCL10': ['miR-503-3p', 'miR-186-5p'],
        'CCR5': ['miR-186-5p', 'miR-669j', 'miR-21-5p', 'miR-146a-5p', 'miR-150-5p', 'miR-136b-5p', 'miR-669k-3p', 'miR-142-3p', 'miR-34a-5p'],
        'CCL4': ['miR-27b-3p', 'miR-27a-3p', 'miR-21-3p', 'miR-467f'],
        'CX3CL1': ['miR-15a-5p', 'miR-322-5p', 'miR-706', 'miR-762', 'miR-665-3p', 'miR-758-3p', 'miR-381-3p'],
        'CXCR4': ['miR-381-3p', 'miR-21-3p', 'miR-467a-5p', 'miR-467h', 'miR-218-5p', 'miR-1a-3p', 'miR-181d-5p', 'miR-206-3p', 'miR-181b-5p', 'miR-9-5p', 'miR-132-3p', 'miR-25-3p', 'miR-467d-5p', 'miR-669k-3p', 'miR-146b-5p', 'miR467b-5p', 'miR-467e-5p', 'miR-467f', 'miR-146a-5p'],
        'CCR-7': ['let-7g-5p', 'miR-23b-3p', 'miR-669p-5p', 'miR-23a-5p', 'let-7e-5p', 'miR-669l-5p', 'miR-15a-5p', 'miR-467e-5p', 'miR-21-5p', 'miR-16-5p', 'let-7d-5p', 'miR-669n', 'miR-98-5p', 'let-7b-5p', 'let-7a-5p', 'let-7i-5p', 'let-7c-5p', 'miR-15b-5p', 'miR-467h'],
        'CXCL12': [
            'mIR-532-5p', 'miR-130b-3p', 'miR-222-3p', 'miR144-3p', 'miR-542-3p', 'miR-149-5p', 'miR-330-3p', 'miR-532-3p', 'miR-3470b', 'miR-125b-5p', 'miR-221-3p', 'miR-19b-3p', 'miR-301b-3p',
            'miR-34b-5p', 'miR-125a-3p', 'miR-126-3p', 'miR-16-1-3p', 'miR-882', 'miR-497-5p', 'miR-26a-5p', 'miR-124-3p', 'miR-26b-5p', 'miR-5620-3p', 'mIR-19a-3p', 'miR-130a-3p', 'miR-690',
            'miR-185-5p', 'miR-31-5p', 'miR-340-5p', 'miR-1843-5p', 'miR-466f-3p', 'miR-301a-3p', 'miR-101a-3p', 'miR-210-3p', 'miR-107-3p', 'miR-706', 'miR-23b-3p', 'miR-146a-5p', 'miR-467f',
            'miR-322-5p', 'miR-15a-5p', 'miR-29b-1-5p', 'let-7e-5p', 'miR-23a-3p', 'miR-338-3p', 'miR-103-3p', 'miR-362-3p', 'let-7g-5p', 'miR-155-5p', 'miR-140-5p', 'miR-122-5p', 'miR-22-3p', 'miR-3470a', 'let-7d-5p'
        ]

    }

    graph = nx.Graph()

    for x in test.simple_connections:
        print(x)
        graph.add_edge(x[0], x[1], attr_dict={'color': "#FF0000"})

    CytoscapeGrapher.showGraph(graph, name='chemokines', title='Chemokines', nodeLabel=lambda x: x, edgeLabel=lambda x: '')
