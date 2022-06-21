import igraph

G = igraph.Graph()

vertex2vid = dict()

with open("/tmp/blood_top", 'r') as fin:

    cnt = 0
    for line in fin:
        line = line.strip()

        aline = line.split('\t')

        fromNode = int(aline[0])
        toNode = int(aline[1])
        weight = float(aline[2])

        if not fromNode in vertex2vid:
            newid = G.add_vertex(fromNode)
            vertex2vid[fromNode] = newid

        if not toNode in vertex2vid:
            newid = G.add_vertex(toNode)
            vertex2vid[toNode] = newid

        fromNodeID = vertex2vid[fromNode]
        toNodeID = vertex2vid[toNode]

        G.add_edge(fromNodeID, toNodeID, weight=weight)

        cnt += 1

        if cnt % 100000 == 0:
            print("Added: " + str(cnt))
            print("Edges: ", len(G.get_edgelist()))
            print("Nodes: ", G.vcount())
            print()
