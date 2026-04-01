import igraph as ig
import sys, time
import os

# Read in DIMACS10 file format
#--------------------------------------------------------------------------------
def rd(pathname,filename, printsense=True):
    """Read graph instance in DIMACS clustering challenge format."""
    if printsense:
        print("DIMACS10 Instance:",filename)
    start = time.time()
    with open(pathname+filename,'r') as infile:
        line = infile.readline()
        vertices = int(line.split()[0])
        edge_count = int(line.split()[1])
        fmt = int(line.split()[2])
        G = ig.Graph()
        if fmt==0:
            if printsense:
                print("#Vertices",vertices)
                print("#Edges",edge_count)
            G.add_vertices(range(vertices))
            
            u = 0
            edges = []
            while (u <= vertices):
                line = infile.readline()
                if (not line.strip()):
                    u = u + 1
#                    print("Vertex", u, " is isolated.")
                elif (line[0] == '%'): 
                    if printsense:
                        print("Skipping comment: ",line)
                else:
                    u = u + 1
                    for word in line.split():
                        if int(word) > u:
                            edges.append((u-1,(int(word)-1)))
            G.add_edges(edges)
            G.vs["name"] = range(1, vertices + 1)
        elif fmt == 1:
            if printsense:
                print("#Vertices",vertices)
                print("#Edges",edge_count)
            G.add_vertices(range(vertices))

            u = 0
            edges = []
            weights = []
            add_weight = False
            while (u <= vertices):
                line = infile.readline()
                if (not line.strip()):
                    u = u + 1
#                    print("Vertex", u, " is isolated.")
                elif (line[0] == '%'): 
                    if printsense:
                        print("Skipping comment: ",line)
                else:
                    u = u + 1
                    for i, word in enumerate(line.split()):
                        if int(word) > u and i % 2 == 0:
                            edges.append((u-1,(int(word)-1)))
                            add_weight = True
                        elif i % 2 == 1 and add_weight:
                            weights.append(int(word))
                            add_weight = False
            G.add_edges(edges, attributes={"weight": weights})
            G.vs["name"] = range(1, vertices + 1)
        else:
           print("DIMACS10 weighted graphs need a new reader!")
    if printsense:
        print(f"Graph read in {time.time() - start:.2f} seconds")
    return G
#--------------------------------------------------------------------------------


# Greedy heuristic function for the lower bound of a maximum clique
#--------------------------------------------------------------------------------
def greedy_lb_max_clq(graph):
    lb_clq = []
    for v in sorted(graph.vs, key=lambda x: graph.degree(x), reverse=True):
        clique = [v["name"]]
        v_nbrs = set(graph.neighbors(v))
        while v_nbrs:
            u = max(v_nbrs, key=lambda x: len(v_nbrs & set(graph.neighbors(x))))
            clique.append(u)
            v_nbrs &= set(graph.neighbors(u))
        if len(clique) > len(lb_clq):
            lb_clq = clique
    return len(lb_clq)
#--------------------------------------------------------------------------------


# Graph peeling by vertex coreness (number for max k-core that includes vertex)
#--------------------------------------------------------------------------------
def core_peel(graph, k):
    core = graph.coreness()
    graph.delete_vertices([v.index for v in graph.vs if core[v.index] < k])
    return graph
#--------------------------------------------------------------------------------


# Testing grounds
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    # for file in os.listdir(r"C:\Users\rackl\ONR-Project\troublesome_graphs\\"):
    #     G = rd(r"C:\Users\rackl\ONR-Project\troublesome_graphs\\", file, printsense=False)
    #     print(file)
    
    G = rd(r"C:\Users\rackl\ONR-Project\troublesome_graphs\\", "lesmis.graph", )
    print(list(G.es["weight"]))