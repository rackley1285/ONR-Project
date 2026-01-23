# Initial Python file for ONR Project
import networkx as nx
import gurobipy as gp
from gurobipy import GRB


# Test data (based on graph from Fig 4b)
A = [
    10,11,
    20,21,
    45,46,47,48,
    54,55,56,57,
    65,66
]
h = [i for i in A]

E_A = []
Ebar_A = []

checked_v = []
for u in A:
    for v in A:
        if u == v or v in checked_v:
            continue

        if abs(u - v) == 1 or (abs(u - v) >= 9 and abs(u - v) <= 11):
            E_A.append((u, v))
        else:
            Ebar_A.append((u, v))
    checked_v.append(u)

G = nx.Graph()
G.add_edges_from(E_A)


# Enumerating maximal cliques
Cscr = list(nx.find_cliques(G))


# Inner optimization problem
inner_model = gp.Model()

x = inner_model.addVars(A, vtype = GRB.BINARY, name = 'vertex_active')
fscr = inner_model.addVars(range(len(Cscr)), vtype=GRB.BINARY, name='clique_active')
z = inner_model.addVars(A, vtype=GRB.BINARY, name='vertex_interdicted')

inner_model.addConstrs((x[u] + x[v] <= 1 for (u, v) in Ebar_A), name="clique_con")
inner_model.addConstrs(
    (gp.quicksum(x[v] for v in set(A) - set(nx.neighbors(G, u))) >= 1 - x[u] for u in A),
    name="maximal_clique_con"
)
inner_model.addConstr(gp.quicksum(fscr) == 1, name='terminal_clique_pick')
inner_model.addConstrs((fscr[C] <= 1 - z[u] for C in range(len(Cscr)) for u in Cscr[C]), name='interdict_clique_con')
inner_model.addConstrs((fscr[C] >= x[u] for C in range(len(Cscr)) for u in Cscr[C]), name='interdict_vertex_con')

inner_model.setObjective(gp.quicksum(h[A.index(u)] * x[u] for u in A), GRB.MAXIMIZE)
inner_model.optimize()

for u in A:
    print(x[u])

print(z)