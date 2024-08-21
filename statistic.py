import networkx as nx
import os
file_path = "./dataset/railway/network.dat"

#print average degree, maximum degree

G = nx.read_edgelist(file_path)
G.remove_edges_from(nx.selfloop_edges(G))
degree = G.degree()
V = G.number_of_nodes()
E = G.number_of_edges()

print("number of nodes", V)
print("number of edges", E)
degree = [d for n, d in degree]
print("average degree", sum(degree)/len(degree), (2*E)/V)
print("maximum degree", max(degree))
print("connceted", nx.is_connected(G))