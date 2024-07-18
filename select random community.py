import networkx as nx
import os
import random



def read_community_file(file_path):
    with open(file_path+"community.dat", 'r') as f:
        communities = []
        for line in f:
            communities.append(list(map(int, line.strip().split())))
    return communities

def read_network_file(file_path):
    with open(file_path+"network.dat", 'r') as f:
        G = nx.Graph()
        for line in f:
            u, v = list(map(int, line.strip().split()))
            G.add_edge(u, v)
    return G


def select_query_nodes(G, community, num_queries):
    degrees = G.degree(community)
    avg_degree = (2*G.number_of_edges())/G.number_of_nodes()
    query_nodes = [node for node, degree in degrees if degree >= avg_degree]
    if len(query_nodes) < num_queries:
        return query_nodes
    else:
        return random.sample(query_nodes, k=num_queries)

if __name__ == "__main__":
    file_path = "./dataset/youtube/"
    GT_communities = read_community_file(file_path)
    G = read_network_file(file_path)

    t = 1
    l = 10
    h = 20
    with open(file_path + "query.txt", 'w') as f:
        for i in range(1,11):
            while True:
                cmty_id = random.randrange(0, len(GT_communities))
                community = GT_communities[cmty_id]
                if len(community) > l:
                    break

            num_queries = random.randint(1, 5)
            # query selection , query는 average degree이상인 node로 선택
            query_nodes = select_query_nodes(G, community, num_queries)
            query_nodes = " ".join(map(str, query_nodes))

            f.write("iteration"+str(i)+"\n")
            f.write("Query\t"+query_nodes+"\n")
            f.write("select cmty id\t"+str(cmty_id)+"\n")
            f.write(
                f"python run.py --algorithm MMA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write(
                f"python run.py --algorithm NGA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write("---------------------------------------------\n")
    f.close()


