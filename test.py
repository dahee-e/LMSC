import networkx as nx
file_path = "./dataset/dblp"
G = nx.read_edgelist(file_path+"/network.dat")
conncected = nx.is_connected(G)
print(conncected)