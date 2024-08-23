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

#
# def LM(G,C,T,v):
#     core = set(C)
#     chain = set(T)
#     chain.add(v)
#     G1 = core.union(chain)
#     G1 = G.subgraph(G1)
#     in_degree = G1.number_of_edges()
#     out_degree = sum(G.degree(u) for u in G1.nodes()) - 2 * in_degree
#     LM = (in_degree) / (out_degree) if (out_degree) > 0 else 0
#     return LM


def LM(G,C,T,v): #core 고려 안함
    core = set(C)
    chain = set(T)
    chain.add(v)
    G1 = core.union(chain)
    core_subgraph = G.subgraph(core).copy()
    chain_subgraph = G.subgraph(chain).copy()
    G1 = G.subgraph(G1)
    internal_edge = G1.number_of_edges() - core_subgraph.number_of_edges() - chain_subgraph.number_of_edges()
    in_degree = chain_subgraph.number_of_edges() + internal_edge
    out_degree = sum(G.degree(u) for u in chain_subgraph.nodes()) - internal_edge
    if out_degree == 0:
        return float('inf')
    LM = (in_degree) / (out_degree)
    return LM

def getChains(G, Core, SS, h1, P):

    P = cu.get_neighbour(G, Core.graph.nodes())
    Z = []
    P = set(P) - set(Core.graph.nodes())
    for u in P:
        chain = [u]
        lm = 0
        while len(chain) < h1:
            neighbors = cu.get_neighbour(G.subgraph(G.nodes()-Core.graph.nodes()), set(chain))
            if len(neighbors) == 0:
                break
            v = max(neighbors, key=lambda w: (LM(G, Core.graph.nodes(),chain, w), -nx.shortest_path_length(G,chain[0], w), -int(w)))

            lm_append = LM(G, Core.graph.nodes(),chain, v)
            if lm_append >= lm:
                chain.append(v)
                lm = lm_append
            else:
                break
        Z.append(chain)
    return Z,P






def findBestSubchain(G, C, Z, t, h_prime, q_nodes):# find best sequence
    C = C.graph
    S_max = []
    lsm_max = -1

    for index, T in enumerate(Z):
        Subchain = []
        Core = C.copy()
        for v in T:
            Core = G.subgraph(set(Core.nodes()).union({v}))
            lsm_current,_,_ = cu.LSM(G, Core, t)
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

        v = max(neighbours, key=lambda x: (cu.LSM(G, C.graph, t, x)[0], -min(nx.shortest_path_length(G, q_node, x) for q_node in q), -int(x)))

        lsm_append,in_degree,out_degree = cu.LSM(G, C.graph.copy(), t, v)
        current_C = list(C.graph.nodes())
        current_C.append(v)
        if lsm_append > best_lsm:
            C.sequence[i + 1] = v

            community = G.subgraph(current_C).copy()

            C = SubgraphData(community, lsm_append, len(community), C.sequence)
            if C.size >= l:
                if weak == True:
                    if in_degree > out_degree:
                        if best_lsm < C.lsm_value:
                            best_lsm = C.lsm_value
                            best_graph = C.graph.copy()
                else:
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
        Z, P = getChains(G, C, S_max, h - C.size, P)
        with open("dataset/karate/chain_exact.txt", 'a') as f:
            f.write("------Chain Identification Procedure------\n")
            f.write(str(len(Z)) + '\n')
            for chain in sorted(Z, key=lambda x: (len(x),x), reverse=True):
                f.write(str(chain) + '\n')
            f.close()

        # STEP 3 : Subchain merge procedure
        S_max = findBestSubchain(G, C, Z, t, h - C.size,q)

        with open("dataset/karate/chain_exact.txt", 'a') as f:
            f.write("------Chain Identification Procedure------\n")
            f.write(str(S_max) + '\n')

        current_C = list(C.graph.nodes())
        current_C.extend(S_max)
        C.sequence[i+1] = S_max
        community = G.subgraph(current_C).copy()
        lsm_current,in_degree,out_degree = cu.LSM(G, community, t)
        C = SubgraphData(community, lsm_current, len(community), C.sequence)



        if C.size >= l:
            if weak == True:
                if in_degree > out_degree:
                    if best_lsm < C.lsm_value:
                        best_lsm = C.lsm_value
                        best_graph = C.graph.copy()
            else:
                if best_lsm < C.lsm_value:
                    best_lsm = C.lsm_value
                    best_graph = C.graph.copy()

        i += 1



    return best_graph,best_lsm,C