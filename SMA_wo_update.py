import networkx as nx
from networkx.algorithms import approximation as ax
import matplotlib.pyplot as plt
import common_utility as cu


class SubgraphData:
    def __init__(self, graph, lsm_value, size,sequence):
        self.graph = graph
        self.lsm_value = lsm_value
        self.size = size
        self.sequence = sequence




def getChains(G, Core,  h1):
    P = set(cu.get_neighbour(G, Core.graph.nodes()))
    Z = []
    C_nodes = set(Core.graph.nodes())
    for u in P:
        chain = [u]
        clsm = cu.CLSM(G, C_nodes, chain,None)
        while len(chain) < h1:
            neighbors = set(cu.get_neighbour(G, set(chain))) - C_nodes
            if len(neighbors) == 0:
                break
            v = max(neighbors, key=lambda w: (
            cu.CLSM(G, Core.graph.nodes(), chain, w), -nx.shortest_path_length(G, chain[0], w), -int(w)))

            clsm_append = cu.CLSM(G, Core.graph.nodes(), chain, v)
            if clsm_append < clsm:
                break
            chain.append(v)
            clsm = clsm_append
        Z.append(chain)
    return Z






def findBestSubchain(G, C, Z, t, h_prime, q_nodes):# find best sequence
    C = C.graph.nodes()
    S_max = []
    lsm_max = -1

    for index, T in enumerate(Z):
        Subchain = []
        Core = C
        for v in T:
            C2 = G.subgraph(set(Core).union({v}))
            lsm_current,_,_ = cu.LSM(G, C2, t)
            Subchain.append(v)
            if lsm_max == lsm_current and S_max != []:
                cand = [S_max, Subchain]
                S_max = min(cand, key=lambda x: (
                min(nx.shortest_path_length(G, q_node, x[0]) for q_node in q_nodes), len(x), int(x[0]))).copy()
            elif lsm_max < lsm_current:
                lsm_max = lsm_current
                S_max = Subchain.copy()
            if len(Subchain) == h_prime:
                break

    return S_max



def run(G, q, l, h, t, weak):
    C = []
    best_lsm = 0
    best_graph = None


    #STEP 0 : preprocessing
    if len(q) == 1:
        initial_graph = G.subgraph(q)
    else:
        initial_graph = G.subgraph(cu.steiner_tree(G, q).nodes())
    print("seed V=", initial_graph.number_of_nodes(), "\tE=", initial_graph.number_of_edges())
    initial_sequence = {0: list(initial_graph.nodes())}
    lsm_current,_,_ = cu.LSM(G, initial_graph, t)
    C = SubgraphData(initial_graph.copy(), lsm_current, len(initial_graph), initial_sequence)
    i = 0

    # STEP 1 : Greedy expansion procedure
    while C.size < h:

        neighbours = cu.get_neighbour(G, C.graph.nodes())

        v = max(neighbours, key=lambda x: (
        cu.LSM(G, C.graph, t, x)[0], -min(nx.shortest_path_length(G, q_node, x) for q_node in q), -int(x)))

        lsm_append, in_degree, out_degree = cu.LSM(G, C.graph.copy(), t, v)
        current_C = list(C.graph.nodes())
        current_C.append(v)
        if lsm_append > best_lsm:
            C.sequence[i + 1] = v

            community = G.subgraph(current_C).copy()

            C = SubgraphData(community, lsm_append, len(community), C.sequence)
            if C.size >= l:
                if best_lsm < C.lsm_value:
                    best_lsm = C.lsm_value
                    best_graph = C.graph.copy()
            i += 1
        else:
            break

    Z = []
    S_max = None
    P = None
    while C.size < h:

        lsm_max = C.lsm_value
        # STEP 2 : Chain identification procedure
        Z = getChains(G, C, h - C.size)
        with open("dataset/chain_exact.txt", 'a') as f:
            f.write("------Chain Identification Procedure------\n")
            f.write(str(len(Z)) + '\n')
            for chain in sorted(Z, key=lambda x: (len(x),x), reverse=True):
                f.write(str(chain) + '\n')
            f.close()

        # STEP 3 : Subchain merge procedure
        S_max = findBestSubchain(G, C, Z, t, h - C.size,q)

        with open("dataset/chain_exact.txt", 'a') as f:
            f.write("------Chain Identification Procedure------\n")
            f.write(str(S_max) + '\n')

        current_C = list(C.graph.nodes())
        current_C.extend(S_max)
        C.sequence[i+1] = S_max
        community = G.subgraph(current_C).copy()
        lsm_current,in_degree,out_degree = cu.LSM(G, community, t)
        C = SubgraphData(community, lsm_current, len(community), C.sequence)

        if C.size >= l:
            if best_lsm < C.lsm_value:
                best_lsm = C.lsm_value
                best_graph = C.graph.copy()

        i += 1



    return best_graph,best_lsm,C