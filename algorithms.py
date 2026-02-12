import networkx as nx

def MD_ordering(graph):
    MD_list = []
    j = 0
    H = graph.copy()
    deg_list = sorted(list(graph.degree), key=lambda x: x[1])
    while j < len(graph.nodes):
        MD_list.append(deg_list[j][0])
        H = H.subgraph(v for v in list(graph) if v not in MD_list)
        j += 1
    
    return MD_list


def degeneracy(graph):
    L = set()
    deg_list = dict(graph.degree)
    max_deg = max(deg_list.values())

    D = [[] for i in range(max_deg + 1)]  #Buckets for possible vertex degrees
    for v, deg in deg_list.items():  #Filling buckets
        D[deg].append(v)
    d = 0  #Degeneracy of the graph
    
    for x in graph.nodes:
        i = 0
        # Scan buckets until a nonempty bucket is found
        while i < max_deg and len(D[i]) == 0:
            i += 1
        v = D[i].pop()
        L.add(v)
        d = i if i > d else d #Update degeneracy with degree of found vertex
            
        
        #Update buckets for neighbors of v
        for u in graph.neighbors(v):
            if u in L:
                continue
            for i in range(len(D)):
                if u in D[i]:
                    D[i].remove(u)
                    D[i-1].append(u)
                    break
        
    return d
    

G = nx.path_graph(4)
print(MD_ordering(G), degeneracy(G))

