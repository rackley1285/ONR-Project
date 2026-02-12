import networkx as nx
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
        self.kappa = list(nx.find_cliques(self.G))

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.restrict_clique(model)
            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        theta = model.cbGetSolution(self.theta)
        z = model.cbGetSolution(self.z)
        G_int = self.G.subgraph([v for v in self.G if model.cbGetSolution(self.z[v]) < 0.5])
        max_clique = nx.max_weight_clique(G_int, weight=None)
        
        if theta + sum(z[i] for i in max_clique[0]) < max_clique[1]:
            model.cbLazy(self.theta + gp.quicksum(self.z[i] for i in max_clique[0]) >= max_clique[1])
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
        z = m.addVars(graph.nodes, vtype=GRB.BINARY, name='z')

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


# Maximal cliques function
#--------------------------------------------------------------------------------
def maximal_clique(graph, p):
    MD = a.MD_ordering(graph)
    d = a.degeneracy(graph)
#--------------------------------------------------------------------------------


# Define test graph
#--------------------------------------------------------------------------------
edges = [(0, i + 1) for i in range(4)
    ] + [(1, i) for i in [2, 3]
    ] + [(2, i) for i in [3, 4]
    ] + [(3, 4)
    ] + [(4, i) for i in [5, 6, 11]
    ] + [(5, i) for i in [6, 7, 8, 10, 11]
    ] + [(7, i) for i in [8, 9, 10]
    ] + [(8, i) for i in [9, 10]
    ] + [(9, 10)]
G = nx.Graph(edges)
nodes = sorted(list(G.nodes))

# Solve problem
int_nodes, max_clq_size = solve_clq_int(G, 2)
#--------------------------------------------------------------------------------


# Visualize graphs
#--------------------------------------------------------------------------------
pos = {
    0: (1 - 1/6, 0.5 - 1/6), 1: (1, 0.5 - 1/6), 2: (1, 0.5 + 1/6), 3: (1 - 1/6, 0.5 + 1/6), 
    4: (1 - 1/3, 0.5), 5: (1/3, 0.5), 6: (0.5, 0.5 + 1/6), 7: (0, 0.5 + 1/6), 
    8: (1/6, 0.5 - 1/6), 9: (0, 0.5 - 1/6), 10: (1/6, 0.5 + 1/6), 11: (0.5, 0.5 - 1/6)
}
plt.figure()
nx.draw_networkx(G, pos=pos, with_labels=True, node_color="#9EC3E2")
plt.savefig("pre-graph.png")
plt.close()

G2 = G.copy()
G2.remove_nodes_from(int_nodes)

plt.figure()
nx.draw_networkx(G2, pos=pos, with_labels=True, node_color="#9EC3E2")
plt.savefig("post-graph.png")
plt.close()