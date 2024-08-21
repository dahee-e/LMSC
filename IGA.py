import networkx as nx

import matplotlib.pyplot as plt
import common_utility as cu
class SubgraphData:
    def __init__(self, graph, lsm_value, size,sequence):
        self.graph = graph
        self.lsm_value = lsm_value
        self.size = size
        self.sequence = sequence # save to dictionary

def run(G, q, l, h, t, weak):
    C = []
    best_lsm = 0
    best_graph = None
    if len(q) == 1:
        initial_graph = G.subgraph(q)
    else:
        initial_graph = G.subgraph(cu.steiner_tree(G, q).nodes())
    print("seed V=", initial_graph.number_of_nodes(), "\tE=", initial_graph.number_of_edges())
    initial_sequence = {0: list(initial_graph.nodes())}
    lsm_current, in_degree, out_degree = cu.LSM(G, initial_graph, t)
    C = SubgraphData(initial_graph.copy(), lsm_current, len(initial_graph), initial_sequence)
    i = 0
    in_degree = 0
    out_degree = 0
    while C.size < h:

        neighbours = cu.get_neighbour(G, C.graph.nodes())

        v = max(neighbours, key=lambda x: (cu.LSM(G, C.graph, t, x)[0], -min(nx.shortest_path_length(G, q_node, x) for q_node in q), -int(x)))
        current_C = list(C.graph.nodes())
        current_C.append(v)

        C.sequence[i+1] = v

        community = G.subgraph(current_C).copy()
        lsm,in_degree,out_degree = cu.LSM(G, community, t)

        C = SubgraphData(community, lsm, len(community),C.sequence)
        i += 1
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




    return best_graph,best_lsm,C