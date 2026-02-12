import networkx as nx

def MD_ordering(graph):
    j = 0
    H = graph.copy()
    deg_list = sorted(list(graph.degree), key=lambda node: node[1])
    MD_list = []
    while j < len(graph.nodes):
        MD_list.append(deg_list[j][0])
        H = H.subgraph(v for v in list(graph) if v not in MD_list)
        j += 1

    return MD_list


G = nx.path_graph(4)
print(MD_ordering(G))
