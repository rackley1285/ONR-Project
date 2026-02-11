import networkx as nx
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt

# Callback class
#--------------------------------------------------------------------------------
class CIPCallback:
    def __init__(self, z, theta, cliques):
        self.z = z
        self.theta = theta
        self.kappa = cliques

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.restrict_clique(model)
            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        theta = model.cbGetSolution(self.theta)
        z = model.cbGetSolution(self.z)

        for K in self.kappa:
            if theta + sum(z[i] for i in K) < len(K):
                model.cbLazy(self.theta + gp.quicksum(self.z[i] for i in K) >= len(K))
                break
#--------------------------------------------------------------------------------

# Solver function
#--------------------------------------------------------------------------------
def solve_clq_int(nodes, cliques, budget):
    with gp.Env() as env, gp.Model(env=env) as m:
        # Variable definition
        theta = m.addVar(vtype=GRB.INTEGER, name='theta')
        z = m.addVars(nodes, vtype=GRB.BINARY, name='z')

        # Constraint definition
        # for K in cliques:
        #    m.addConstr(theta + gp.quicksum(z[i] for i in K) >= len(K), name='clq_max')
        m.addConstr(gp.quicksum(z.values()) <= budget, name='budget')
        m.addConstr(theta >= 0, name='theta_lb')

        # Optimization
        m.Params.LazyConstraints = 1
        cb = CIPCallback(z, theta, cliques)
        m.setObjective(theta, GRB.MINIMIZE)
        m.optimize(cb)

        u_int = [i for i in z if z[i].getAttr(GRB.Attr.X) > 0.5]
        return u_int, m.ObjVal
#--------------------------------------------------------------------------------

# Maximal cliques function
#--------------------------------------------------------------------------------
def maximal_clique(graph):
    pass
#--------------------------------------------------------------------------------

# Test graph
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

# Enumerating maximal cliques (BK algo) 
    ### REPLACE WITH INTEGER PROGRAM
kappa = list(nx.find_cliques(G))

# Solve problem
int_nodes, max_clq_size = solve_clq_int(nodes, kappa, 2)

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