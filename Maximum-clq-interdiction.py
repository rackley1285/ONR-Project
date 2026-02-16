import igraph as ig
ig.config['plotting.backend'] = 'matplotlib'
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import algorithms as a

# Callback class
#--------------------------------------------------------------------------------
class CIPCallback:
    def __init__(self, z, theta, graph):
        self.z = z
        self.theta = theta
        self.G = graph

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.restrict_clique(model)
            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        theta_hat = model.cbGetSolution(self.theta)
        z_hat = model.cbGetSolution(self.z)
        V_bar = [v for v in self.G.vs["name"] if z_hat[v] < 0.5]
        V_bar = self.G.vs.select(name_in=V_bar).indices
        G_int = self.G.induced_subgraph(V_bar)
        cliques = G_int.maximal_cliques()
        max_clique = max(cliques, key=len)
        
        if theta_hat < len(max_clique):
            model.cbLazy(self.theta + gp.quicksum(self.z[G_int.vs[i]["name"]] for i in max_clique) >= len(max_clique))
            
#--------------------------------------------------------------------------------


# Solver function
#--------------------------------------------------------------------------------
def solve_clq_int(graph, budget):
    with gp.Env() as env, gp.Model(env=env) as m:
        # Variable definition
        theta = m.addVar(vtype=GRB.INTEGER, name='theta')
        z = m.addVars(graph.vs["name"], vtype=GRB.BINARY, name='z')

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
#G = rd("/workspaces/ONR-Project/testbed/", "power.graph")
G = a.rd(r"C:\Users\rackl\ONR-Project\testbed\\", r"CIP_example.graph")


# Solve problem
int_nodes, max_clq_size = solve_clq_int(G, 2)
V2 = G.vs.select(name_notin=int_nodes)

# Visualize graphs
layout = G.layout(layout="graphopt")
name_layout = {G.vs[i]["name"]: layout.coords[i] for i in G.vs.indices}


ig.plot(G, target="pre-graph.png", layout=layout, vertex_label = G.vs["name"])
plt.close()

G2 = G.induced_subgraph(V2)
ig.plot(G2, target="post-graph.png", layout=ig.Layout([name_layout[v["name"]] for v in V2]), vertex_label = V2["name"])
plt.close()