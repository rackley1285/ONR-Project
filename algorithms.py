import networkx as nx
import igraph as ig

# Create a minimum-degree ordered list
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

# Walteros-Buchanan maximum clique algorithm
#--------------------------------------------------------------------------------
def WB_max_clq(graph, p):
    V = MD_ordering(graph)
    indices = range(len(V))
    d = degeneracy(graph)
    n = len(V)

    D = []
    D_rns = []
    for i in indices:
        r_neighborhood = [V[j] for j in indices[i:] if graph.has_edge(V[i], V[j])]
        rdeg = len(r_neighborhood)
        if i < n - d and rdeg >= d - p:
            D.append(V[i])
            D_rns.append(r_neighborhood)
    
    for v_rn in D_rns:
        G_bar = graph.subgraph()

#--------------------------------------------------------------------------------


# Find random clique of size q
#--------------------------------------------------------------------------------
def find_clique(graph, q):
    v1 = 0
    for i in range(graph.vcount()):
        if graph.degree(i) == q:
            v1 = i
            break
    
    C = [v1]
    C_nbrs = [(v1, neighbor) for neighbor in graph.neighborhood(v1)]
    for v1, neighbor in C_nbrs:
        if graph.degree(neighbor) >= q and neighbor not in C:
            C.append(neighbor)
    
    i = 1
    while len(C) > q and i < len(C):
        for j in range(1, len(C)):
            if not graph.are_adjacent(C[i], C[j]) and C[i] != C[j]:
                C.pop(j)
        i += 1
    
    if len(C) < q:
        return []
    else:
        return sorted(C)
#--------------------------------------------------------------------------------


# Read in DIMACS10 file format
#--------------------------------------------------------------------------------
def rd(pathname,filename):
    """Read graph instance in DIMACS clustering challenge format."""
    print("DIMACS10 Instance:",filename)
    with open(pathname+filename,'r') as infile:
        line = infile.readline()
        vertices = int(line.split()[0])
        edges = int(line.split()[1])
        fmt = int(line.split()[2])
        G = ig.Graph()
        if fmt==0:
            print("#Vertices",vertices)
            print("#Edges",edges)
            for i in range(vertices):
                G.add_vertex(name = str(i+1))
            
            u = 0
            while (u <= vertices):
                line = infile.readline()
                if (not line.strip()):
                    u = u + 1
#                    print("Vertex", u, " is isolated.")
                elif (line[0] == '%'): 
                    print("Skipping comment: ",line)
                else:
                    u = u + 1
                    for word in line.split():
                        G.add_edge(u-1,(int(word)-1))
        else:
            print("DIMACS10 weighted graphs need a new reader!")
    return G
#--------------------------------------------------------------------------------


# Testing grounds
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    G = rd("/workspaces/ONR-Project/testbed/", "CIP_example.graph").simplify()
    print(find_clique(G, 4))
    #WB_max_clq(G, 0)

