# Evaluate NMI between detected communities and ground truth
import networkx as nx
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics import f1_score
from sklearn.metrics.cluster import adjusted_rand_score
import argparse
import glob
import os


# Reading the Ground-Truth Community Data
def read_GT_community(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        community = []
        for line in lines:
            community.append(list(map(int, line.strip().split())))
    return community

def read_community_results(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        line = lines[1]
        nodes_str = line.split("\t")[1].strip()
        nodes = list(map(int, nodes_str.split()))
        seconds = lines[0].split("\t")[1].strip()
    return nodes,seconds



# Calculating NMI Score
def calculate_nmi(true_communities, detected_communities,G):
    # Create a set of all nodes
    all_nodes = list(map(int,G.nodes()))

    # Create node-to-community mapping
    true_communities_dict = {node: 1 for node in true_communities}
    detected_communities_dict = {node: 1 for node in detected_communities}

    # Map nodes to their community labels
    labels_true = [true_communities_dict[node] if node in true_communities_dict else 0 for node in all_nodes]
    labels_detected = [detected_communities_dict[node] if node in detected_communities_dict else 0 for node in
                       all_nodes]

    # Calculate NMI and ARI
    return normalized_mutual_info_score(labels_true, labels_detected)

def calculate_ari(true_communities, detected_communities,G):
    # Create a set of all nodes
    all_nodes = list(map(int, G.nodes()))

    # Create node-to-community mapping
    true_communities_dict = {node: 1 for node in true_communities}
    detected_communities_dict = {node: 1 for node in detected_communities}

    # Map nodes to their community labels
    labels_true = [true_communities_dict[node] if node in true_communities_dict else 0 for node in all_nodes]
    labels_detected = [detected_communities_dict[node] if node in detected_communities_dict else 0 for node in
                       all_nodes]
    # Calculate NMI and ARI
    return adjusted_rand_score(labels_true, labels_detected)
def calculate_F1score(true_communities, detected_communities,G):
    # Flatten the lists and create label vectors
    all_nodes = list(map(int, G.nodes()))

    # Create node-to-community mapping
    true_communities_dict = {node: 1 for node in true_communities}
    detected_communities_dict = {node: 1 for node in detected_communities}

    # Map nodes to their community labels
    labels_true = [true_communities_dict[node] if node in true_communities_dict else 0 for node in all_nodes]
    labels_detected = [detected_communities_dict[node] if node in detected_communities_dict else 0 for node in
                       all_nodes]

    # Calculate NMI and ARI
    return f1_score(labels_true, labels_detected,  average='macro')

def eval(file_path,f,cmty_id    ):
    # Replace or append file extensions as necessary to construct paths
    community_file_path = file_path+f

    network_file_path = file_path+"network.dat"
    ground_truth_file_path = file_path+"community.dat"

    G = nx.read_edgelist(network_file_path)
    detected_communities,seconds  = read_community_results(community_file_path)
    true_communities = read_GT_community(ground_truth_file_path)
    true_communities = true_communities[cmty_id]
    nmi_score = calculate_nmi(true_communities, detected_communities,G)
    fscore = calculate_F1score(true_communities, detected_communities,G)
    ari_score = calculate_ari(true_communities, detected_communities,G)
    return nmi_score, fscore, ari_score, seconds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate community search algorithm')
    parser.add_argument('--network', default="./dataset/amazon/",
                        help='a folder name containing network.dat')

    args = parser.parse_args()
    print("network ", args.network)

    # file_path = os.listdir(args.network)


    file_path = ["MMA_q_['456995']_l_10_h_20_t_1.0.dat"]
    file_path.append(file_path[0].replace("MMA", "NGA"))

    cmty_id = 49481


    for f in file_path:
        if f.startswith("NGA") or f.startswith("MMA"):
            print("file name: ",f)
            nmi_score, fscore, ari_score,seconds = eval(args.network,f,cmty_id)
            print("NMI Score: ", str(nmi_score))
            print("ARI Score: ", str(ari_score))
            print("F1 Score: ", str(fscore))

            with open(args.network+"accuracy.txt", 'a') as f1:
                f1.write("file name: "+f+"\n")
                f1.write("seconds: "+str(seconds)+"\n")
                f1.write("NMI Score: "+str(nmi_score)+"\n")
                f1.write("ARI Score: "+str(ari_score)+"\n")
                f1.write("F1 Score: " + str(fscore) + "\n")
                f1.write("\n")
            f1.close()

