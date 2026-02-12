import networkx as nx


#Create a minimum-degree ordered list
#--------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------


# Calculate the degeneracy of a graph
#--------------------------------------------------------------------------------
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
            u_deg = deg_list[u]
            deg_list[u] -= 1
            D[u_deg].remove(u)
            D[u_deg - 1].append(u)
        
    return d
#--------------------------------------------------------------------------------
 

#Testing grounds
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    G = nx.path_graph(4)
    print(MD_ordering(G), degeneracy(G))

