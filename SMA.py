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


def LM(G,C,T,v):
    core = set(C)
    chain = set(T)
    chain.add(v)
    G1 = core.union(chain)
    G1 = G.subgraph(G1)
    in_degree = G1.number_of_edges()
    out_degree = sum(G.degree(u) for u in G1.nodes()) - 2 * in_degree
    LM = (in_degree) / (out_degree) if (out_degree) > 0 else 0
    return LM


def getChains(G, Core, SS, h1, Z, P):
    if len(Z) == 0:
        P = cu.get_neighbour(G, Core.graph.nodes())
    else:
        P = set(cu.get_neighbour(G, SS)) - set(cu.get_neighbour(G,P))

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



def updateChains(G, C, SS, h1, Z):
    core = C.graph.nodes()
    one_hop = set(cu.get_neighbour(G, SS))
    j = 0
    while j < len(Z):
        chain = Z[j]
        if chain[0] in SS:
            Z.remove(chain)
            continue
        R = set(one_hop)  - set(core)
        two_hop = cu.get_neighbour_include_set(G, R)
        if set(two_hop).intersection(set(chain)):
            T_new = []
            candidates = {chain[0]}
            while len(T_new) < h1:
                u = max(candidates,key=lambda x: (LM(G, core, T_new, x), -nx.shortest_path_length(G,chain[0], x), -int(x)))
                T_new = T_new + [u]
                candidates = candidates - set([u])
                if any(u == item for item in two_hop):
                    candidates = candidates.union(R.intersection(set(cu.get_neighbour(G, [u]))))
                if any(u == item for item in chain):
                    next_node = int(chain.index(u))+1
                    next_node = set([chain[next_node]])
                    next_node -= set(core)
                    candidates = candidates.union(set(next_node))

                candidates -= set(T_new)
            chain = T_new.copy()
            Z[j] = chain
        j += 1
    return Z



def findBestSubchain(G, C, Z, t, h_prime):# find best sequence
    C = C.graph
    S_max = []
    lsm_max = 0

    for index, T in enumerate(Z):
        Subchain = []
        Core = C.copy()
        for v in T:
            Core = G.subgraph(set(Core.nodes()).union({v}))
            lsm_current,_,_ = cu.LSM(G, Core, t)
            Subchain.append(v)
            if lsm_max < lsm_current:
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
        Z, P = getChains(G, C, S_max, h - C.size, Z, P)

        with open("./dataset/karate/chain.txt", 'a') as f:
            f.write("------Chain Identification Procedure------\n")
            f.write(str(len(Z)) + '\n')
            for chain in Z:
                f.write(str(chain) + '\n')
            f.close()

        # STEP 3 : Subchain merge procedure
        S_max = findBestSubchain(G, C, Z, t, h - C.size)

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

        # STEP 4 : Chain update procedure
        Z =  updateChains(G, C, S_max, h - C.size, Z)



    return best_graph,best_lsm,C