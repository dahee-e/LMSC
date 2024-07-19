import networkx as nx
import os
import random
import argparse


def read_community_file(file_path,l):
    with open(file_path+"community.dat", 'r') as f:
        communities = []
        for line in f:
            community = list(map(int, line.strip().split()))
            if len(community) > l-1:
                communities.append(community)
    return communities

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

    l = args.l
    h = args.h
    t = args.t

    file_path = args.network
    GT_communities = read_community_file(file_path,l)
    G = read_network_file(file_path)

    with open(file_path + "query.txt", 'w') as f:
        for i in range(1,11):
            cmty_id = random.randrange(0, len(GT_communities))
            community = GT_communities[cmty_id]


            num_queries = random.randint(1, 5)

            query_nodes = select_query_nodes(community, num_queries)
            query_nodes = " ".join(map(str, query_nodes))

            f.write("iteration"+str(i)+"\n")
            f.write("Query\t"+query_nodes+"\n")
            f.write("Community\t"+str(len(community))+"\n")
            f.write("-------------Given a community of size------------- \n ")
            f.write(
                f"python run.py --algorithm MMA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write(
                f"python run.py --algorithm NGA --q {query_nodes} --l {str(l)} --h {str(h)} --t {str(t)} --network {file_path}network.dat\n")
            f.write("---------------------------------------------\n")
            f.write("---------------------------------------------\n")
    f.close()


