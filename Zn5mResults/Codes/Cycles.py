import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


nodes = np.loadtxt("SysInfo.txt")[1].astype(int)
nodes = [i for i in range(0, nodes)]

edges = np.genfromtxt(r"Water/Edges.csv", delimiter=",", dtype=int)

G = nx.Graph()
G.add_nodes_from(nodes, nodetype=int)
G.add_edges_from(edges, edgetype=int)

degree = [v for k, v in nx.degree(G)]

# Clusters
components = list(nx.connected_components(G))
components = list(sorted(components, key=lambda x: len(x), reverse=True))

comp_sizes = [len(comp) for comp in components]

prob = np.zeros(10)

for i in comp_sizes:
    if i == 1:
        prob[0] += i
    elif i == 2:
        prob[1] += i
    elif i == 3:
        prob[2] += i
    elif i == 4:
        prob[3] += i
    elif i >= 5 and i <= 20:
        prob[4] += i
    elif i >= 21 and i <= 50:
        prob[5] += i
    elif i >= 51 and i <= 100:
        prob[6] += i
    elif i >= 101 and i <= 300:
        prob[7] += i
    elif i >= 301 and i <= 500:
        prob[8] += i
    elif i >= 501:
        prob[9] += i

prob = prob * 100 / nx.number_of_nodes(G)

heights, bins, plot = plt.hist(
    comp_sizes, bins=[*range(1, 5), 20, 50, 100, 300, 500, 5000]
)
with open("Water/clusterFreq.txt", "a") as f:
    print(*heights, file=f)

with open("Water/clusterPerc.txt", "a") as f:
    print(*np.round(prob, 4), file=f)

with open("Water/clusterSizes.txt", "a") as f:
    print(*comp_sizes, file=f)


# Largest conncted component
with open("Water/lcc.txt", "a") as f:
    print(max(comp_sizes), file=f)

# Cluster fraction
with open("Water/clusterFraction.txt", "a") as f:
    print(1 - (degree.count(1) / len(nodes)), file=f)

# Normalised number of clusters
large_components = [comp for comp in components if len(comp) >= 2]
fraction_large_components = 2 * len(large_components) / len(nodes)

with open("Water/normClusters.txt", "a") as f:
    print(fraction_large_components, file=f)


# H2O Cycles


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
