import igraph as ig
ig.config['plotting.backend'] = 'matplotlib'
import gurobipy as gp
from gurobipy import GRB
import matplotlib.pyplot as plt
import algorithms as a
import time, os, math
import pandas as pd

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
        del solver

        # Returning solution statistics and data
        total_t = time.time() - t0 if m.Status != GRB.TIME_LIMIT and m.SolCount > 0 else "TL"
        u_int = [i for i in z if z[i].getAttr(GRB.Attr.X) > 0.5] if m.SolCount > 0 else []
        #u_del = [i for i in x if x[i].getAttr(GRB.Attr.X) > 0.5]
        obj = m.ObjVal if m.SolCount > 0 else -1
        return len(u_int), obj, m.NodeCount, cb.CB, cb.LC, total_t, cb.CB_time, "Optimal" if m.Status == 2 else m.Status
#--------------------------------------------------------------------------------


# Maximum clique solver class
#--------------------------------------------------------------------------------
class MCSolver:
    def __init__(self, graph):
        self.G = graph
        self.env = gp.Env()
        self.m = gp.Model(env=self.env)
        self.m.Params.OutputFlag = 0
        self.m.Params.MIPFocus = 1
        
        self.components = graph.components()
        self.x = self.m.addVars(graph.vs["name"], vtype=GRB.BINARY, name="x") # Vertices included in the clique
        self.f = self.m.addVars(range(len(self.components)), vtype=GRB.BINARY, name="f")

        # Find connected components first, then add constraints for each component (when non-neighbors AND in same component)
        for i, component in enumerate(self.components):
            comp_adj = {graph.vs[v]["name"]: set(graph.vs[u]["name"] for u in graph.neighbors(v)) for v in component}
            for v, nbrs in comp_adj.items():
                self.m.addConstr(self.x[v] <= self.f[i], name="compnt_memb")
                for u in component:
                    u = graph.vs[u]["name"]
                    if u not in nbrs and u != v:
                        self.m.addConstr(self.x[v] + self.x[u] <= 1, name="clq_adj")
        
        self.m.addConstr(gp.quicksum(self.f) == 1, name="1_compt")

        self.m.setObjective(gp.quicksum(self.x), GRB.MAXIMIZE)
        self.m.update()

    def solve(self, deleted, theta):
        # Update variable upper bounds
        peeled_set = {v["name"] for v in self.G.vs if self.G.coreness()[v.index] < theta - 1}
        deleted_set = set(deleted)
        for v in self.x:
            self.x[v].ub = 0 if v in deleted_set or v in peeled_set else 1

        # Optimize model and return solution
        self.m.Params.BestObjStop = theta + 0.5
        self.m.optimize()
        #print(f"Status: {self.m.Status}, SolCount: {self.m.SolCount}")
        
        # Check status to see if it reached optimality or target

        clique = [i for i in self.x if self.x[i].X > 0.5] if self.m.Status != GRB.TIME_LIMIT and self.m.SolCount > 0 else "TL"
        return clique

    def __del__(self):
        self.m.dispose()
        self.env.dispose()
#--------------------------------------------------------------------------------


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
        self.start = time.time()

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
                import traceback
                print(f"Error: {Exception}")
                print(traceback.format_exc())
    
    def restrict_clique(self, model):
        # Callback solution
        theta_hat = math.ceil(model.cbGetSolution(self.theta))
        x_hat = model.cbGetSolution(self.x)
        
        if self.solver == "enum":
            # Creating interdicted subgraph
            V_bar = [v for v in self.G.vs["name"] if x_hat[v] < 0.5] #Selects names of non-deleted vertices
            V_bar = self.G.vs.select(name_in=V_bar) #Selects vertex objects for interdicted graph by name attribute
            G_int = self.G.induced_subgraph(V_bar)

            # Enumerate maximal cliques
            # cliques = G_int.maximal_cliques()
            # max_clique = max(cliques, key=len) if len(cliques) > 0 else []

            # clique = G_int.cliques(min = theta_hat + 1, max_results = 1)
            # max_clique = clique[0] if len(clique) > 0 else []

            max_clique = G_int.largest_cliques()[0]

            if theta_hat < len(max_clique):
                model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[G_int.vs[i]["name"]] for i in max_clique))
                self.LC += 1
        
        else:
            # Finding a maximum clique for lazy constraint
            V_del = [v for v in self.G.vs["name"] if x_hat[v] > 0.5] # List of deleted vertices
            max_clique = self.solver.solve(V_del, theta_hat) # Max clique in remaining graph
            if max_clique == "TL":
                return
            elif theta_hat < len(max_clique):
                model.cbLazy(self.theta >= len(max_clique) - gp.quicksum(self.x[i] for i in max_clique))
                self.LC += 1
            self.solver.m.Params.TimeLimit = 3600 - (time.time() - self.start)
#--------------------------------------------------------------------------------


# Function to run solver on a graph file
#--------------------------------------------------------------------------------
def max_clq_int(path, file, budget, separation):
    G = a.rd(path, file, printsense=False)
    filename = file.split(".")[0]
    budget = budget if budget >= 1 else math.ceil(budget * G.vcount())
    return [filename] + list(solve_clq_int(G, budget, separation))
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
if __name__ == "__main__":
    ex_col = ["Graph G", "z(V)", "theta", "#BC", "#CB", "#LC", "Total time (s)", "CB time (s)", "Status"]
    sheets = []
    budgets = [2, 0.05, 0.1]

    # print(max_clq_int(r"C:\Users\rackl\ONR-Project\testbed\\", "football.graph", 0.1, 0))

    for b in budgets:
        dataMIP = []
        dataIG = []
        for file in os.listdir(r"C:\Users\rackl\ONR-Project\testbed\\"):
            print(file)
            dataMIP.append(max_clq_int(r"C:\Users\rackl\ONR-Project\testbed\\", file, b, 1))
            dataIG.append(max_clq_int(r"C:\Users\rackl\ONR-Project\testbed\\", file, b, 0))
        
        sheets.append((pd.DataFrame(dataMIP, columns = ex_col), b, "MIP"))
        sheets.append((pd.DataFrame(dataIG, columns = ex_col), b, "IG"))

    with pd.ExcelWriter(r"C:\Users\rackl\ONR-Project\clq_int_statistics.xlsx") as writer:
        for df, b, solver in sheets:
            df.to_excel(writer, sheet_name = solver + ", b = " + str(b), index=False)


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