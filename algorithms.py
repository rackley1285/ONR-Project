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

def degeneracy(graph, MD = None):
    if MD == None:
        MD = MD_ordering(graph)
    
    D = []
    for v in MD:
        Nv = list(graph.neighbors(v))
        rnbrhd = [n for n in MD[MD.index(v):] if n in Nv]
        D.append(len(rnbrhd))
    
    return max(D)
    

G = nx.path_graph(4)
print(degeneracy(G))

