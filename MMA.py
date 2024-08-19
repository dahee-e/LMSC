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
    core = set(C.nodes())
    chain = set(T)
    chain.add(v)
    G1 = core.union(chain)
    G1 = G.subgraph(G1)
    in_degree = G1.number_of_edges()
    out_degree = sum(G.degree(u) for u in G1.nodes()) - 2 * in_degree
    LM = (in_degree) / (out_degree) if (out_degree) > 0 else 0
    return LM


def getChains(G, Core, SS, h1, M):
    if M is None:
        candidates = cu.get_neighbour(G, Core.graph.nodes())
    else:
        candidates = set(cu.get_neighbour(G, SS)) - set(Core.graph.nodes())

    candidates = candidates - set(Core.graph.nodes())
    sequence_set = set()

    for u in candidates:
        chain = [u]
        lm = 0
        while len(chain) < h1:
            v = max(candidates, key=lambda w: LM(G, Core.graph,chain, w))
            lm_append = LM(G, Core.graph,chain, v)
            if lm_append >= lm:
                chain.append(v)
                lm = lm_append
        if len(chain) > 1:
            M.append(chain)
    return M



def updateChains(G, C, SS, h1, M):
    one_hop = cu.get_neighbour(G, SS)
    for chain in M:
        R = one_hop - set(chain.graph.nodes())
        two_hop = cu.get_neighbour(G, one_hop)
        T_new = []
        candidates = set(M[0])
        i = 0
        while len(T_new) < h1:
            u = max(R, key=lambda x: LM(G, chain.graph, T_new, x))
            T_new = T_new + [u]
            candidates = candidates - {u}
            if M[i] in two_hop:
                candidates = candidates.union(one_hop.intersection(cu.get_neighbour(G, M[i].graph.nodes())))
            elif u != M[i]:
                candidates = candidates.union(cu.get_neighbour(G, {u}) - C.nodes())
            else:
                candidates = candidates.union(M[i+1])
                i += 1
    return M



def findBestSequence(G, C, M, t, h_prime):# find best sequence , best sequence에 포함된 chain은 삭제하고 chain set return
    C = C.graph
    S_max = []
    lsm_max = 0

    for index, T in enumerate(M):
        Subseq = []
        Core = C.copy()
        for v in T.sequence:
            Core = G.subgraph(set(Core.nodes()).union({v}))
            lsm_current = cu.LSM(G, Core, t)
            Subseq.append(v)
            if lsm_max < lsm_current:
                lsm_max = lsm_current
                S_max = Subseq.copy()
            if len(Subseq) == h_prime:
                break

    return S_max


def run(G, q, l, h, t):
    C = []

    #STEP 0 : preprocessing
    if len(q) == 1:
        initial_graph = G.subgraph(q)
    else:
        initial_graph = G.subgraph(cu.steiner_tree(G, q).nodes())
    print("seed V=", initial_graph.number_of_nodes(), "\tE=", initial_graph.number_of_edges())

    C.append(SubgraphData(initial_graph, cu.LSM(G, initial_graph, t), len(initial_graph),list(initial_graph.nodes())))

    i = 0
    while C[i].size != h:

        #STEP 1 : Greedy expansion procedure
        lsm_current = C[i].lsm_value
        while C[i].size < h:

            neighbours = cu.get_neighbour(G, C[i].graph.nodes())

            v = max(neighbours, key=lambda x: cu.LSM(G, C[i].graph, t, x))
            lsm_append = cu.LSM(G, C[i].graph.copy(), t, v)
            if lsm_current < lsm_append:
                sequence = C[i].sequence.copy()
                sequence.append(v)

                community = G.subgraph(sequence).copy()

                C.append(SubgraphData(community, lsm_append, len(community),sequence))
                lsm_current = lsm_append
                i += 1
            else:
                break

        if C[i].size == h-1: # 넣을 수 있는 노드가 하나밖에 없을때, 종료
            break



        while C[i].size < h:
            # STEP 2 : Chain identification procedure
            M = getChains(G, C[i], h - C[i].size, M)

            # STEP 3 : Chain merge procedure
            S_max, _ = findBestSequence(G, C[i], M, t, h - C[i].size)
            sequence = C[i].sequence.copy()
            sequence.extend(S_max)
            community = G.subgraph(sequence).copy()
            lsm_current = cu.LSM(G, community, t)
            C.append(SubgraphData(community, lsm_current, len(community), sequence))
            i += 1

            # STEP 4 : Chain update procedure
            M =  updateChains(G, C, S_max, h - C[i].size, M)





    filtered_sorted = sorted((res for res in C if (l <= res.size) and (res.size <= h)), key=lambda x: x.lsm_value,reverse=True)
    return C, filtered_sorted[0].graph, filtered_sorted[0].lsm_value if filtered_sorted else None