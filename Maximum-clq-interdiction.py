import igraph as ig
ig.config['plotting.backend'] = 'matplotlib'
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import algorithms as a
import time
import os
import pandas as pd

# Clique interdiction callback class
#--------------------------------------------------------------------------------
class CIPCallback:
    def __init__(self, x, theta, graph, solver):
        self.x = x
        self.theta = theta
        self.G = graph
        self.solver = solver
        self.CB = 0
        self.LC = 0
        self.CB_time = 0

    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                t0 = time.time()
                self.restrict_clique(model)
                
                # Update statistics
                self.CB += 1
                self.CB_time += time.time() - t0

            except Exception:
                model.terminate()
    
    def restrict_clique(self, model):
        # Callback solution
        theta_hat = model.cbGetSolution(self.theta)
        x_hat = model.cbGetSolution(self.x)
        
        if self.solver == "enum":
            # Creating interdicted subgraph
            V_bar = [v for v in self.G.vs["name"] if x_hat[v] < 0.5] #Selects names of non-interdicted vertices
            V_bar = self.G.vs.select(name_in=V_bar) #Selects vertex objects for interdicted graph by name attribute
            G_int = self.G.induced_subgraph(V_bar)

            # Enumerate maximal cliques
            cliques = G_int.maximal_cliques()
            max_clique = max(cliques, key=len)
            if theta_hat < len(max_clique):
                model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[G_int.vs[i]["name"]] for i in max_clique))
                self.LC += 1
        
        else:
            # Finding a maximum clique for lazy constraint
            V_del = [v for v in self.G.vs["name"] if x_hat[v] > 0.5] # List of deleted vertices
            max_clique = self.solver.solve(V_del, theta_hat) # Max clique in remaining graph
            if theta_hat < len(max_clique):
                model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[i] for i in max_clique))
                self.LC += 1
#--------------------------------------------------------------------------------


# Clique interdiction solver function
#--------------------------------------------------------------------------------
def solve_clq_int(graph, budget, separation):
    with gp.Env() as env, gp.Model(env=env) as m:
        # Python variables
        sep_procs = ["enum", "MIP"]
        t0 = time.time()
        solver = MCSolver(graph) if separation == 1 else sep_procs[separation]
        adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}

        # Gurobi variables
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
        m.Params.OutputFlag = 0
        m.Params.LazyConstraints = 1
        m.Params.TimeLimit = 3600
        cb = CIPCallback(x, theta, graph, solver)
        m.setObjective(theta, GRB.MINIMIZE)
        m.optimize(cb)
        total_t = time.time() - t0 if m.Status != GRB.TIME_LIMIT else "TL"

        # Returning solution statistics and data
        u_int = [i for i in z if z[i].getAttr(GRB.Attr.X) > 0.5]
        u_del = [i for i in x if x[i].getAttr(GRB.Attr.X) > 0.5]
        return len(u_int), m.ObjVal, m.NodeCount, cb.CB, cb.LC, total_t, cb.CB_time
#--------------------------------------------------------------------------------


# Maximum clique solver class
#--------------------------------------------------------------------------------
class MCSolver:
    def __init__(self, graph):
        self.G = graph
        self.env = gp.Env()
        self.m = gp.Model(env=self.env)
        self.m.Params.OutputFlag = 0

        self.x = self.m.addVars(graph.vs["name"], vtype=GRB.BINARY, name="x") # Vertices included in the clique

        # Find connected components first, then add constraints for each component (when non-neighbors AND in same component)

        adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}
        vertices = graph.vs["name"]
        for v, nbrs in adj.items():
            for u in vertices:
                if u not in nbrs and u != v:
                    self.m.addConstr(self.x[v] + self.x[u] <= 1, name="clq_adj")
        
        self.m.setObjective(gp.quicksum(self.x), GRB.MAXIMIZE)
        self.m.update()

    def solve(self, interdicted, theta):
        # Update variable upper bounds
        interdicted_set = set(interdicted)
        for v in self.x:
            self.x[v].ub = 0 if v in interdicted_set else 1

        # Optimize model and return solution
        self.m.Params.BestObjStop = theta + 0.5
        self.m.optimize()
        #print(f"Status: {self.m.Status}, SolCount: {self.m.SolCount}")
        
        # Check status to see if it reached optimality or target

        clique = [i for i in self.x if self.x[i].X > 0.5]
        return clique

    def __del__(self):
        self.m.dispose()
        self.env.dispose()
#--------------------------------------------------------------------------------


# Maximum clique solver function
#--------------------------------------------------------------------------------
def solve_max_clq(graph, theta, interdicted):
    with gp.Env() as env, gp.Model(env=env) as m:
        x = m.addVars(graph.vs["name"], vtype=GRB.BINARY, name="x") # Vertices in the clique
        for v in interdicted:
            x[v].ub = 0
        
        adj = {v['name']: set(graph.vs[u]["name"] for u in graph.neighbors(v.index)) for v in graph.vs}
        vertices = graph.vs["name"]
        for v, nbrs in adj.items():
            for u in vertices:
                if u not in nbrs and u != v:
                    m.addConstr(x[v] + x[u] <= 1, name="clq_adj")

        m.Params.BestObjStop = theta
        m.setObjective(gp.quicksum(x), GRB.MAXIMIZE)
        m.optimize()

        clique = [i for i in x if x[i].getAttr(GRB.Attr.X) > 0.5]
        return clique
#--------------------------------------------------------------------------------


# Function to run solver on a graph file
#--------------------------------------------------------------------------------
def max_clq_int(path, file, separation):
    G = a.rd(path, file, printsense=False)
    filename = file.split(".")[0]
    return [filename] + list(solve_clq_int(G, 2, separation))
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
if __name__ == "__main__":
    ex_col = ["Graph G", "z(V)", "theta", "#BC", "#CB", "#LC", "Total time (s)", "CB time (s)"]
    
    data = []
    for file in os.listdir(r"C:\Users\rackl\ONR-Project\testbed\\"):
        data.append(max_clq_int(r"C:\Users\rackl\ONR-Project\testbed\\", file, 1))
    
    df = pd.DataFrame(data, columns = ex_col)
    df.to_excel(r"C:\Users\rackl\ONR-Project\MIP_statistics.xlsx", index=False)


    data = []
    for file in os.listdir(r"C:\Users\rackl\ONR-Project\testbed\\"):
        data.append(max_clq_int(r"C:\Users\rackl\ONR-Project\testbed\\", file, 0))
    
    df = pd.DataFrame(data, columns = ex_col)
    df.to_excel(r"C:\Users\rackl\ONR-Project\Enum_statistics.xlsx", index=False)

#G = rd("/workspaces/ONR-Project/testbed/", "power.graph")


'''
# Visualize graphs

V2 = G.vs.select(name_notin=del_nodes) #Selects vertex objects needed for interdicted graph by name attribute
G2 = G.induced_subgraph(V2)

layout = G.layout(layout="graphopt")
name_layout = {G.vs[i]["name"]: layout.coords[i] for i in G.vs.indices}


ig.plot(G, target="pre-graph.png", layout=layout, vertex_label = G.vs["name"])
plt.close()

ig.plot(G2, target="post-graph.png", layout=ig.Layout([name_layout[v["name"]] for v in V2]), vertex_label = V2["name"])
plt.close()
'''