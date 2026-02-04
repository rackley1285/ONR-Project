#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edge formulation for the maximum clique problem.

Implemented and solved using Gurobi Optimization solver.
Reads DIMACS 2nd clique challenge format instances and 
unweighted DIMACS 10th clustering challenge format (METIS format) instances.
@author: baski363
"""
from gurobipy import *
import sys, time
import networkx as nx


# SET INPUT PARAMETERS HERE
pathname = "/Users/baski363/Library/CloudStorage/OneDrive-OklahomaAandMSystem/CODES/Data/"
filename = "karate.graph"
fileformat = "dimacs10"
timelimit = 10

# ############################################################################

def read_dimacs2(pathname,filename):
    """Read graph instance in DIMACS clique challenge format."""
    print("DIMACS2 Instance:",filename)
    with open(pathname+filename,'r') as infile:
        for line in infile:
            if(not line.strip()):
                continue               
            elif(line.split()[0] == "c"):
                continue
            elif(line.split()[0] == "p"):
                vertices = int(line.split()[2])
                edges = int(line.split()[3])
                print("#Vertices",vertices)
                for i in range(vertices):
                   G.add_node(i+1)
                print("#Edges",edges)
            elif(line.split()[0] == "e"):
                u = int(line.split()[1])
                v = int(line.split()[2])
                G.add_edge(u,v)


def read_dimacs10(pathname,filename):
    """Read graph instance in DIMACS clustering challenge format."""
    print("DIMACS10 Instance:",filename)
    with open(pathname+filename,'r') as infile:
        line = infile.readline()
        vertices = int(line.split()[0])
        edges = int(line.split()[1])
        fmt = int(line.split()[2])
        if fmt==0:
            print("#Vertices",vertices)
            print("#Edges",edges)
            for i in range(vertices):
                G.add_node(i+1)
                
            u = 0;
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
                        G.add_edge(u,(int(word)))
        else:
            print("DIMACS10 weighted graphs need a new reader!")

# FUNCTION DEFINITIONS ARE ABOVE THIS LINE
# ############################################################################

start = time.time()
G = nx.Graph()
if fileformat == "dimacs2": 
    read_dimacs2(pathname,filename)
elif fileformat == "dimacs10":  
    read_dimacs10(pathname,filename)
else:
    sys.exit("File format unknown.")

try:  
    # Model & Solve
    model = Model('MaxClqEdgeFormulation')
    n = G.order()
    # Create variables
    x = model.addVars(G.nodes(), vtype=GRB.BINARY, obj=1, name='x')
    model.ModelSense = GRB.MAXIMIZE

    # Add non-edge constraints
    model.addConstrs(x[u]+x[v] <= 1 for u,v in nx.non_edges(G))

    # Set time limit    
    model.setParam('TimeLimit', float(timelimit)) 
    
    model.optimize()
    
    # Write the model in a .lp format (turn-off for large models)
    #model.write(filename+'_maxclq.lp')
    
    # Retrieve and print result
    if model.Status == GRB.Status.OPTIMAL or GRB.Status.SUBOPTIMAL:
        termstat = {2 : 'OPTIMAL', 9 : 'SUBOPTIMAL'}
        print('\nTermination was ' + termstat[model.Status])
        print('\nMax clique:')
        solution = model.getAttr('X', x)
        sol = ''
        for i in G.nodes():
            if solution[i] > 0.5:
                sol += str(i)+" "
        print(sol)
    end = time.time()
    print('Elapsed time: {t:.2f} seconds'.format(t= end - start))
except GurobiError:
    print('\nError reported')