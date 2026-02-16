import networkx as nx
import igraph as ig
import sys, time


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
            G.add_vertices(range(vertices))
            
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
            G.vs["name"] = range(1, vertices + 1)
        else:
            print("DIMACS10 weighted graphs need a new reader!")
    print(f"Graph read in {time.time() - start:.2f} seconds")
    return G
#--------------------------------------------------------------------------------


# Testing grounds
#--------------------------------------------------------------------------------
if __name__ == "__main__":
    G = rd(r"C:\Users\rackl\ONR-Project\testbed\\", r"cond-mat.graph")
    print(G)

