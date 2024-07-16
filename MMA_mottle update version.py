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



def f1(G,C,T,v,t): # mottle 만 고려
    #T is the current mottle
    #v is the node to be added
    core = set(C.nodes())
    mottle = set(T)
    mottle.add(v)
    mottle = G.subgraph(mottle)
    print(mottle.nodes())
    in_degree = mottle.number_of_edges()
    out_degree = sum(G.degree(u) for u in mottle.nodes()) - 2 * in_degree
    G1 = core.union(mottle)
    G1 = G.subgraph(G1)
    internal_edges = sum(G1.degree(v) for v in (G1.nodes()-core)) - 2 * in_degree

    LM = in_degree / (out_degree-internal_edges) if (out_degree-internal_edges) > 0 else 0
    return LM



def f2(G,C,T,v,t): # mottle & core 전부 고려 = 기존과 동일
    # T is the current mottle
    # v is the node to be added
    core = set(C.nodes())
    mottle = set(T)
    mottle.add(v)
    mottle = G.subgraph(mottle)
    print(mottle.nodes())
    in_degree = mottle.number_of_edges()
    in_degree_core = C.number_of_edges()
    out_degree = sum(G.degree(u) for u in mottle.nodes()) - 2 * in_degree
    out_degree_core = sum(G.degree(u) for u in core) - 2 * in_degree_core
    G1 = core.union(mottle)
    G1 = G.subgraph(G1)
    internal_edges = sum(G1.degree(v) for v in (G1.nodes() - core)) - 2 * in_degree
    LM = (in_degree+in_degree_core+internal_edges) / (out_degree+out_degree_core-(2*internal_edges)) if (out_degree+out_degree_core-(2*internal_edges)) > 0 else 0
    return LM



def f3(G,C,T,v,t): # mottle 고려 & mottle과 core간의 edge 고려
    # T is the current mottle
    # v is the node to be added
    core = set(C.nodes())
    mottle = set(T)
    mottle.add(v)
    mottle = G.subgraph(mottle)
    print(mottle.nodes())
    in_degree = mottle.number_of_edges()
    out_degree = sum(G.degree(u) for u in mottle.nodes()) - 2 * in_degree
    G1 = core.union(mottle)
    G1 = G.subgraph(G1)
    internal_edges = sum(G1.degree(v) for v in (G1.nodes() - core)) - 2 * in_degree

    LM = (in_degree +internal_edges)/ (out_degree - internal_edges) if (out_degree - internal_edges) > 0 else 0
    return LM


# def get_neighbour_candidates(G,C,T):
#     neighbours = set()
#     for x in T:
#         neighbours.update(G.neighbors(x))
#     neighbours -= set(C.nodes())
#     neighbours -= set(T.nodes())
#     return list(neighbours)
#

def getMottles(G, Core, h1, t, f):
    M = []
    sequence_set = set()
    core_cmty = Core
    for u in cu.get_neighbour(G, core_cmty.graph):
        lsm = 0
        node = [u]
        while len(node) < h1:
            candidates = set(cu.get_neighbour(G,G.subgraph(node))) - set(core_cmty.graph.nodes())
            v = max(candidates, key=lambda w: f(G, core_cmty.graph,node, w, t))
            lsm_append = f(G, core_cmty.graph,node, v, t)

            if lsm_append <= lsm:
                sequence = tuple(sorted(node))
                if sequence not in sequence_set:
                    community = G.subgraph(node).copy()
                    M.append(SubgraphData(community, lsm, len(community), node.copy()))
                    sequence_set.add(sequence)
                break
            lsm = lsm_append
            node.append(v)
        if len(node) == h1 :
            sequence = tuple(sorted(node))
            if sequence not in sequence_set:
                community = G.subgraph(node).copy()
                M.append(SubgraphData(community, lsm, len(community), node.copy()))
                sequence_set.add(sequence)
    return M


# def getMottles(G, Core, h1, t, f):
#     M = []
#     core_cmty = Core
#     for u in cu.get_neighbour(G, core_cmty.graph):
#         lsm = 0
#         node = [u]
#         while len(node) < h1:
#             candidates = set(cu.get_neighbour(G,G.subgraph(node))) - set(core_cmty.graph.nodes())
#             v = max(candidates, key=lambda w: f(G, core_cmty.graph,node, w, t))
#             lsm_append = f(G, core_cmty.graph,node, v, t)
#
#             if lsm_append <= lsm:
#                 community = G.subgraph(node).copy()
#                 M.append(SubgraphData(community, lsm, len(community), node.copy()))
#                 break
#             lsm = lsm_append
#             node.append(v)
#         if len(node) == h1 :
#             community = G.subgraph(node).copy()
#             M.append(SubgraphData(community, lsm, len(community), node.copy()))
#     return M

# def updateMottles(G, C,C_prev, h1, t, f):
#     M = []
#     neighbours = cu.get_neighbour(G, C) - cu.get_neighbour(G, C_prev)
#     for u in neighbours:
#         T = [u]
#         while len(T) < h1:
#             candidates = cu.get_neighbour(G,T) - set(C.nodes())
#             v = max(candidates, key=lambda w: f(G, C, T, w, t))
#             T.append(v)
#         M.append(T)
#     return M

# def findBestSequence(G, C, M, t): best sequence가 포함된 mottle은 삭제하고 mottle set return
#     C = C.graph
#     S_max = []
#     lsm_max = 0
#     for T in M:
#         Subseq = []
#         lsm_mtle = 0
#         for v in T:
#             if lsm_mtle < cu.LSM(G, C, v, t):
#                 Subseq.append(v)
#             else:
#                 break
#         if lsm_max < lsm_mtle:
#             lsm_max = lsm_mtle
#             S_max = Subseq
#
#     return S_max, M

def findBestSequence(G, C, M, t):
    C = C.graph
    S_max = []
    lsm_max = 0
    for T in M:
        Subseq = []
        lsm_mtle = 0
        for v in T:
            if lsm_mtle < cu.LSM(G, C, v, t):
                Subseq.append(v)
            else:
                break
        if lsm_max < lsm_mtle:
            lsm_max = lsm_mtle
            S_max = Subseq

    return S_max, M

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

            neighbours = cu.get_neighbour(G, C[i].graph)

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


        #STEP 2 : Mottle identification procedure
        M = getMottles(G,C[i], h-C[i].size,t,f2)

        #STEP 3 : Mottle merge procedure


        while C[i].size < h:
            S_max, _ = findBestSequence(G, C[i], M, t)
            sequence = C[i].sequence.copy()
            sequence.extend(S_max)
            community = G.subgraph(sequence).copy()
            lsm = cu.LSM(G, community, t)
            C.append(SubgraphData(community, lsm, len(community), sequence))
            M_prime = updateMottles(G, C[i+1], C[i], h - C[i+1].size, t, f1)
            M = M + M_prime
            i += 1


    filtered_sorted = sorted((res for res in C if (l <= res.size) and (res.size <= h)), key=lambda x: x.lsm_value,reverse=True)
    return C, filtered_sorted[0].graph if filtered_sorted else None