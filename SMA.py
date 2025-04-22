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


def CLSM(G,C,T,v):
    core = set(C)
    chain = set(T)
    if v is not None:
        chain.add(v)
    G1 = core.union(chain)
    core_subgraph = G.subgraph(core).copy()
    chain_subgraph = G.subgraph(chain).copy()
    G1 = G.subgraph(G1)
    internal_edge = G1.number_of_edges() - core_subgraph.number_of_edges() - chain_subgraph.number_of_edges()
    in_degree = chain_subgraph.number_of_edges()
    out_degree = sum(G.degree(u) for u in chain_subgraph.nodes()) - 2 * in_degree - internal_edge
    if out_degree == 0:
        return float('inf')
    CLSM_score = (in_degree + internal_edge) / (out_degree)
    return CLSM_score



def getChains(G, Core, h1):
    P = set(cu.get_neighbour(G, Core.graph.nodes()))
    Z = []
    C_nodes = set(Core.graph.nodes())
    for u in P:
        chain = [u]
        clsm = CLSM(G, C_nodes, chain,None)
        while len(chain) < h1:
            neighbors = set(cu.get_neighbour(G, set(chain))) - C_nodes
            if len(neighbors) == 0:
                break
            v = max(neighbors, key=lambda w: (CLSM(G, Core.graph.nodes(),chain, w), -nx.shortest_path_length(G,chain[0], w), -int(w)))

            clsm_append = CLSM(G, Core.graph.nodes(),chain, v)
            if clsm_append < clsm:
                break
            chain.append(v)
            clsm = clsm_append
        Z.append(chain)
    return Z, P





def updateChains(G, C, S, Z, P_prime, h1): # G: 전체 그래프, C: 현재 그래프, SS: subchain, h1: 최대 길이, Z: chain set
    C_nodes = set(C.graph.nodes())
    Z_prime = []


    for chain in Z:
        new_chain = []
        index = len(chain)
        pivot = chain[0]


        # Update rule 2: Remove chain if pivot is in S
        if pivot in S:
            continue

        # Update rule 1: Truncate chain if it exceeds the size constraint h-|C|
        if len(chain) > h1:
            chain = chain[:h1]

        # Update rule 3: violate the disjointness conditions
        overlap_nodes = [i for i, v in enumerate(chain) if v in S]
        if overlap_nodes:
            index = min(index, overlap_nodes[0])
            new_chain = chain[:index]


        # Update rule 5: candidate nodes with improved CLSM gain
        N_Z = set(cu.get_neighbour_include_set(G, set(chain))) #N(Z,G)
        N_S = set(cu.get_neighbour(G, S)) #N(S,G)
        R = N_Z.intersection(N_S) - C_nodes
        if R:
            N_R = set(cu.get_neighbour_include_set(G, R)) - C_nodes
            for u in chain:
                if u in N_R:
                    index_tmp = chain.index(u)
                    index = min(index, index_tmp+1)
                    break

            new_chain = chain[:index]
        chain_set = set(chain)
        R_prime = chain_set.intersection(N_S) - C_nodes
        if R_prime:
            for u in chain:
                if u in R_prime:
                    index_tmp = chain.index(u)
                    index = min(index, index_tmp+1)
                    break
            new_chain = chain[:index]


        if len(new_chain) == 0:
            Z_prime.append(chain)
            continue
        clsm = CLSM(G, C_nodes, new_chain,None)
        while len(new_chain) < h1:  # Extend chain greedily
            candidates = set(cu.get_neighbour(G, set(new_chain))) - C_nodes
            if not candidates:
                break
            best_v = max(candidates, key=lambda w: (cu.CLSM(G, C_nodes, new_chain, w), -nx.shortest_path_length(G, pivot, w), -int(w)))
            clsm_append = CLSM(G, C_nodes, new_chain, best_v)
            if clsm_append < clsm:
                break
            clsm = clsm_append
            new_chain.append(best_v)

        Z_prime.append(new_chain)

    # Update Rule 6: Add new chains from new pivots
    new_pivots = set(cu.get_neighbour(G, S)) - C_nodes - P_prime
    P_prime = P_prime.union(new_pivots)

    for pivot in new_pivots:
        new_chain = [pivot]
        clsm = CLSM(G, C_nodes, new_chain, None)
        while len(new_chain) < h1:
            candidates = set(cu.get_neighbour(G, set(new_chain))) - C_nodes
            if not candidates:
                break
            best_v = max(candidates, key=lambda w: (CLSM(G, C_nodes, new_chain, w), -nx.shortest_path_length(G, pivot, w), -int(w)))
            clsm_append = CLSM(G, C_nodes, new_chain, best_v)
            if clsm_append < clsm:
                break
            clsm = clsm_append
            new_chain.append(best_v)
        Z_prime.append(new_chain)

    return Z_prime, P_prime





def findBestSubchain(G, C, Z, t, h_prime, q_nodes):# find best sequence
    C = C.graph.nodes()
    S_max = []
    lsm_max = -1

    for index, T in enumerate(Z):
        Subchain = []
        Core = set(C)
        for v in T:
            Core = G.subgraph(set(Core).union({v}))
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
    P = set()
    # STEP 2 : Chain identification procedure
    Z,P = getChains(G, C, h - C.size)


    while C.size < h:


        # STEP 3 : Subchain merge procedure
        S_max = findBestSubchain(G, C, Z, t, h - C.size,q)

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
        if C.size == h:
            break

        # STEP 4 : Chain update procedure
        Z,P =  updateChains(G, C, S_max,Z,P, h - C.size)



    return best_graph,best_lsm,C