import igraph as ig
ig.config['plotting.backend'] = 'matplotlib'
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import algorithms as a

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


# Callback class
#--------------------------------------------------------------------------------
class CIPCallback:
    def __init__(self, z, theta, graph):
        self.z = z
        self.theta = theta
        self.G = graph
        self.kappa = graph.largest_cliques()

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.restrict_clique(model)
            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        theta = model.cbGetSolution(self.theta)
        z = model.cbGetSolution(self.z)
        V_bar = [v for v in range(self.G.vcount()) if model.cbGetSolution(self.z[v]) < 0.5]
        G_int = self.G.induced_subgraph(V_bar)
        max_clique = G_int.largest_cliques()[0]
        
        if theta + sum(z[i] for i in max_clique) < len(max_clique):
            model.cbLazy(self.theta + gp.quicksum(self.z[i] for i in max_clique) >= len(max_clique))
        # for K in self.kappa: #FIX THIS: Find a maximal clique to add a lazy constraint for
        #     if theta + sum(z[i] for i in K) < len(K):
        #         model.cbLazy(self.theta + gp.quicksum(self.z[i] for i in K) >= len(K))
                
#--------------------------------------------------------------------------------


# Solver function
#--------------------------------------------------------------------------------
def solve_clq_int(graph, budget):
    with gp.Env() as env, gp.Model(env=env) as m:
        # Variable definition
        theta = m.addVar(vtype=GRB.INTEGER, name='theta')
        z = m.addVars(range(graph.vcount()), vtype=GRB.BINARY, name='z')

        # Constraint definition
        m.addConstr(gp.quicksum(z.values()) <= budget, name='budget')
        m.addConstr(theta >= 0, name='theta_lb')

        # Optimization
        m.Params.LazyConstraints = 1
        cb = CIPCallback(z, theta, graph)
        m.setObjective(theta, GRB.MINIMIZE)
        m.optimize(cb)

        u_int = [i for i in z if z[i].getAttr(GRB.Attr.X) > 0.5]
        return u_int, m.ObjVal
#--------------------------------------------------------------------------------


# Define test graphs
#--------------------------------------------------------------------------------
G = rd("/workspaces/ONR-Project/testbed/", "netscience.graph").simplify()

# Solve problem
int_nodes, max_clq_size = solve_clq_int(G, 10)
V2 = [i for i in range(G.vcount()) if i not in int_nodes]

'''# Visualize graphs
layout = G.layout(layout="graphopt")
ig.plot(G, target="pre-graph.png", layout=layout, vertex_label = range(G.vcount()))
plt.close()

G2 = G.induced_subgraph(V2)
ig.plot(G2, target="post-graph.png", layout=ig.Layout([layout[i] for i in V2]), vertex_label = V2)
plt.close()'''