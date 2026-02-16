import networkx as nx
import igraph as ig
import sys, time

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


# Read in DIMACS10 file format
#--------------------------------------------------------------------------------
def rd(pathname,filename):
    """Read graph instance in DIMACS clustering challenge format."""
    print("DIMACS10 Instance:",filename)
    start = time.time()
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
            edges = []
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
                        if int(word) > u:
                            edges.append((u-1,(int(word)-1)))
            G.add_edges(edges)
        else:
            print("DIMACS10 weighted graphs need a new reader!")
    print(f"Graph read in {time.time() - start:.2f} seconds")
    return G
#--------------------------------------------------------------------------------


# Testing grounds
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    G = rd(r"C:\Users\rackl\ONR-Project\testbed\\", r"cond-mat.graph").simplify()
    print(find_clique(G, 2))
    #WB_max_clq(G, 0)

