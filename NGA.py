import networkx as nx

import matplotlib.pyplot as plt
import common_utility as cu
class SubgraphData:
    def __init__(self, graph, lsm_value, size,sequence):
        self.graph = graph
        self.lsm_value = lsm_value
        self.size = size
        self.sequence = sequence



# def LSM(G, C, t):
#     in_degree = C.number_of_edges()
#     out_degree = 0
#     for u in C.nodes():
#         out_degree += G.degree(u)
#     out_degree -= in_degree * 2
#     LM = in_degree / out_degree
#
#     return LM * (1 / pow(len(C), t))





def run(G, q, l, h, t):
    C = []
    if len(q) == 1:
        initial_graph = G.subgraph(q)
    else:
        initial_graph = G.subgraph(cu.steiner_tree(G, q).nodes())
    print("seed V=", initial_graph.number_of_nodes(), "\tE=", initial_graph.number_of_edges())

    C.append(SubgraphData(initial_graph, cu.LSM(G, initial_graph, t), len(initial_graph),list(initial_graph.nodes())))

    i = 0
    while C[i].size < h:

        neighbours = cu.get_neighbour(G, C[i].graph)

        v = max(neighbours, key=lambda x: cu.LSM(G, C[i].graph, t, x))

        sequence = C[i].sequence.copy()
        sequence.append(v)

        community = G.subgraph(sequence).copy()
        lsm = cu.LSM(G, community, t)

        C.append(SubgraphData(community, lsm, len(community),sequence))
        i += 1

    # Filter and sort the C by LSM value, returning the best graph that satisfies the size constraints [l, h]
    filtered_sorted = sorted((res for res in C if (l <= res.size) and (res.size<= h)), key=lambda x: x.lsm_value, reverse=True)
    return C,filtered_sorted[0].graph, filtered_sorted[0].lsm_value if filtered_sorted else None