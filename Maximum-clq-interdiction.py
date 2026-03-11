import igraph as ig
ig.config['plotting.backend'] = 'matplotlib'
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import algorithms as a

# Clique interdiction callback class
#--------------------------------------------------------------------------------
class CIPCallback:
    def __init__(self, x, theta, graph):
        self.x = x
        self.theta = theta
        self.G = graph

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.restrict_clique(model)
            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        # Callback solution
        theta_hat = model.cbGetSolution(self.theta)
        x_hat = model.cbGetSolution(self.x)
        
        # Creating interdicted subgraph
        V_bar = [v for v in self.G.vs["name"] if x_hat[v] < 0.5] #Selects names of non-interdicted vertices
        V_bar = self.G.vs.select(name_in=V_bar) #Selects vertex objects for interdicted graph by name attribute
        G_int = self.G.induced_subgraph(V_bar)
        
        # Finding a maximum clique for lazy constraint
        max_clique = solve_max_clq(G_int, theta_hat)
        if theta_hat < len(max_clique):
            model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[i] for i in max_clique))
        
        # cliques = G_int.maximal_cliques()
        # max_clique = max(cliques, key=len)
        # if theta_hat < len(max_clique):
            # model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[G_int.vs[i]["name"]] for i in max_clique))           
#--------------------------------------------------------------------------------


# Clique interdiction solver function
#--------------------------------------------------------------------------------
def solve_clq_int(graph, budget):
    with gp.Env() as env, gp.Model(env=env) as m:
        n = graph.vcount()
        adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}

        # Variable definitions
        theta = m.addVar(vtype=GRB.INTEGER, name='theta') #Maximum size of remaining cliques
        z = m.addVars(graph.vs["name"], vtype=GRB.BINARY, name='z') #Interdicted vertices
        x = m.addVars(graph.vs["name"], vtype=GRB.BINARY, name='x') #Deleted vertices

        # Constraint definitions
        m.addConstr(gp.quicksum(z) <= budget, name='budget')
        m.addConstr(theta >= 0, name='theta_lb')
        for u, nbrs in adj.items():
            for v in nbrs:
                m.addConstr(x[v] >= z[u], name='clq_int')
            m.addConstr(gp.quicksum(z[v] for v in nbrs) >= x[u], name='clq_nonint')
            m.addConstr(x[u] >= z[u], name='del_int')

        # Optimization
        m.Params.LazyConstraints = 1
        cb = CIPCallback(x, theta, graph)
        m.setObjective(theta, GRB.MINIMIZE)
        m.optimize(cb)

        # Returning vertex names (as stored in z) and objective value
        u_int = [i for i in z if z[i].getAttr(GRB.Attr.X) > 0.5]
        u_del = [i for i in x if x[i].getAttr(GRB.Attr.X) > 0.5]
        return u_int, u_del, m.ObjVal
#--------------------------------------------------------------------------------


# Maximum clique callback class
#--------------------------------------------------------------------------------
class MCCallback:
    def __init__(self, x, graph):
        self.G = graph
        self.x = x
        self.adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}
    
    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.enforce_adjacency(model)
            except Exception:
                model.terminate()
    
    def enforce_adjacency(self, model):
        pass
#--------------------------------------------------------------------------------


# Maximum clique solver function
#--------------------------------------------------------------------------------
def solve_max_clq(graph, theta):
    with gp.Env() as env, gp.Model(env=env) as m:
        x = m.addVars(graph.vs["name"], vtype=GRB.BINARY, name="x") # Vertices in the clique
        
        adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}
        for v, nbrs in adj.items():
            for i in range(graph.vcount()):
                u = graph.vs[i]["name"]
                if u not in nbrs and u != v:
                    m.addConstr(x[v] + x[u] <= 1, name="clq_adj")

        # m.Params.LazyConstraints = 1
        # cb = MCCallback(x, graph)
        m.Params.BestObjStop = theta
        m.setObjective(gp.quicksum(x), GRB.MAXIMIZE)
        m.optimize()

        clique = [i for i in x if x[i].getAttr(GRB.Attr.X) > 0.5]
        return clique
#--------------------------------------------------------------------------------


# Define test graphs
#--------------------------------------------------------------------------------
#G = rd("/workspaces/ONR-Project/testbed/", "power.graph")
G = a.rd(r"C:\Users\rackl\ONR-Project\testbed\\", r"CIP_example.graph")

# c = solve_max_clq(G)
# print(c)
# print(max(G.maximal_cliques(), key=len))

# Solve problem
int_nodes, del_nodes, max_clq_size = solve_clq_int(G, 2)
V2 = G.vs.select(name_notin=del_nodes) #Selects vertex objects needed for interdicted graph by name attribute
G2 = G.induced_subgraph(V2)
print(int_nodes)
'''
# Visualize graphs

layout = G.layout(layout="graphopt")
name_layout = {G.vs[i]["name"]: layout.coords[i] for i in G.vs.indices}


ig.plot(G, target="pre-graph.png", layout=layout, vertex_label = G.vs["name"])
plt.close()

ig.plot(G2, target="post-graph.png", layout=ig.Layout([name_layout[v["name"]] for v in V2]), vertex_label = V2["name"])
plt.close()
'''