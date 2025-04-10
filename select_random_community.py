import networkx as nx
import os
import random
import argparse
import networkx as nx


def read_community_file(file_path):
    with open(file_path+"community.dat", 'r') as f:
        communities = []
        for line in f:
            communities.append(list(map(int, line.strip().split())))
    return communities

def read_node_community_file(file_path):
    with open(file_path + "community.dat", 'r') as f:
        node_cmty = {}
        for i, line in enumerate(f):
            for node in line.strip().split():
                node = int(node)
                if node not in node_cmty:
                    node_cmty[node] = []
                node_cmty[node].append(i)
    return node_cmty


def read_network_file(file_path):
    with open(file_path+"network.dat", 'r') as f:
        G = nx.Graph()
        for line in f:
            u, v = list(map(int, line.strip().split()))
            G.add_edge(u, v)
    return G


def select_query_nodes(community, num_queries):
    query_nodes = community.copy()

    if len(query_nodes) < num_queries:
        return query_nodes
    else:
        return random.sample(query_nodes, k=num_queries)

def select_single_query(G, num_queries):
    query_nodes = nx.k_truss(G, k=3).nodes()
    return random.sample(query_nodes, k=num_queries)


if __name__ == "__main__":

    t = 1
    l = 10
    h = 20
    parser = argparse.ArgumentParser()
    parser.add_argument('--l', type=int, default=l,
                        help='lower bound l')
    parser.add_argument('--h', type=int, default=h,
                        help='upper bound h')
    parser.add_argument('--t', type=float, default=t,
                        help='threshold tau for local sketch modularity`')
    parser.add_argument('--network', default="./dataset/karate/",)
    args = parser.parse_args()
    overlap = []
    # overlap = ['amazon', 'dblp', 'youtube',"livejournal","orkut"]
    l = args.l
    h = args.h
    t = args.t

    file_path = args.network
    GT_communities = read_community_file(file_path)
    node_cmty = read_node_community_file(file_path)
    G = read_network_file(file_path)

    with open(file_path + "query.txt", 'w') as f:
        for i in range(10):


            num_queries = 1
            query_nodes = select_single_query(G, num_queries)
            query_nodes = " ".join(map(str, query_nodes))
            #community 는 query 가 포함된 GT_communities에서 찾음
            cmty_id = random.choice(node_cmty[int(query_nodes)])
            community = GT_communities[cmty_id]


            f.write("iteration "+str(i)+"\n")
            f.write("Query\t"+query_nodes+"\n")
            f.write("Community size\t"+str(len(community))+"\n")
            f.write("Community id\t"+str(cmty_id)+"\n")
            f.write("-------------Given a community of size------[C/2, 3C/2]------- \n")
            l = len(community)//2
            h = min(3*len(community)//2,G.number_of_nodes())
            f.write("l\t"+str(l)+"\n")
            f.write("h\t"+str(h)+"\n")
            f.write("t\t"+str(t)+"\n")
            f.write(
                f"python run.py --algorithm SMA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write(
                f"python run.py --algorithm IGA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write("---------------------------------------------\n")
            f.write("---------------------------------------------\n")
    f.close()


