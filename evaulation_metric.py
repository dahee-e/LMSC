#edge density
#coefficient clustering
#modularity
import networkx as nx
import os

def read_network_file(file_path):
    with open(file_path+"network.dat", 'r') as f:
        G = nx.Graph()
        for line in f:
            u, v = list(map(int, line.strip().split()))
            G.add_edge(u, v)
    return G

def read_community_file(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        line = lines[1]
        nodes_str = line.split("\t")[1].strip()
        nodes = list(map(int, nodes_str.split()))
    return nodes

def degreeSum(G, nodes):
    sum = 0
    for u in nodes:
        sum = sum + G.degree(u)
    return sum

def get_modularity(G, C):
    E = G.number_of_edges()

    lc = C.number_of_edges()
    var1 = lc / E
    var2 = (degreeSum(G, C.nodes) / (2.0 * E)) ** 2
    mod = var1 - var2
    return mod


if __name__ == "__main__":
    file_path = "./dataset/amazon/"
    G = read_network_file(file_path)
    file = os.listdir(file_path)
    file = [f for f in file if f.startswith("SMA") or f.startswith("IGA")]

    for file_name in file:
        C = read_community_file(file_path+file_name)
        C = G.subgraph(C)

        community_size = len(C.nodes())
        edge_density = nx.density(C)
        clustering_coefficient = nx.average_clustering(C)
        modularity = get_modularity(G, C)

        print("File name: ", file_path+file_name)
        print("Community size: ", community_size)
        print("Edge density: ", edge_density)
        print("Clustering coefficient: ", clustering_coefficient)
        print("Modularity: ", modularity)
        print("--------------------------------------------------------")

