import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# H2O Cycles

nodes = np.loadtxt("SysInfo.txt")[1].astype(int)
nodes = [i for i in range(0, nodes)]

edges = np.genfromtxt(r"Water/Edges.csv", delimiter=",", dtype=int)

G = nx.Graph()
G.add_nodes_from(nodes, nodetype=int)
G.add_edges_from(edges, edgetype=int)

cycles = nx.cycle_basis(G)
sizes = [len(x) for x in cycles]
sizes = np.array(sizes)
sizes = -np.sort(-sizes)

with open("Water/Cycles.txt", "a") as f:
    print(*sizes, file=f)

prob = np.zeros(11)

for i in sizes:
    if i < 11:
        prob[i] += 1

prob = prob / len(sizes)
prob = prob[3:]

with open("Water/CyclesProbs.txt", "a") as f:
    print(*np.round(prob, 4), file=f)


components = nx.connected_components(G)
largest_component = max(components, key=len)
LCC = G.subgraph(largest_component)


# Combined
numW = np.loadtxt("SysInfo.txt")[1].astype(int)
nodes = [i for i in range(0, numW)]

e = np.genfromtxt(r"Combined/EdgesT.csv", delimiter=",", dtype=int)

e_no_tt = e[e[:, 0] < numW]


G = nx.Graph()
G.add_nodes_from(nodes, nodetype=int)
G.add_edges_from(e_no_tt, edgetype=int)

cycles = nx.cycle_basis(G)
sizes = [len(x) for x in cycles]
sizes = np.array(sizes)
sizes = -np.sort(-sizes)

with open("Combined/Cycles.txt", "a") as f:
    print(*sizes, file=f)

prob = np.zeros(11)

for i in sizes:
    if i < 11:
        prob[i] += 1

prob = prob / len(sizes)
prob = prob[3:]

with open("Combined/CyclesProbs.txt", "a") as f:
    print(*np.round(prob, 4), file=f)


components = nx.connected_components(G)
largest_component = max(components, key=len)
LCC = G.subgraph(largest_component)


with open("Combined/Largest.txt", "a") as f:
    print(len(LCC), file=f)
