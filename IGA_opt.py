import networkx as nx
import common_utility as cu
import time

class SubgraphData:
    def __init__(self, graph, lsm_value, size, sequence):
        self.graph = graph
        self.lsm_value = lsm_value
        self.size = size
        self.sequence = sequence


def run(G, q, l, h, t,steiner_type):
    if len(q) == 1:
        initial_graph = G.subgraph(q)
    else:
        initial_graph = G.subgraph(cu.steiner_tree(G, q,steiner_type).nodes())

    print("seed V=", initial_graph.number_of_nodes(), "\tE=", initial_graph.number_of_edges())

    start_time = time.time()
    lsm_calc = cu.LSMCalculator(G, initial_graph.nodes(), t)
    print("LSM init time:", time.time() - start_time)

    initial_sequence = {0: list(initial_graph.nodes())}
    C = SubgraphData(
        initial_graph.copy(),
        lsm_calc.get_lsm(),
        len(initial_graph),
        initial_sequence
    )

    best_lsm = -1
    best_graph = None
    i = 0
    lsm_current = C.lsm_value

    while C.size < h:
        neighbours = cu.get_neighbour(G, C.graph.nodes())
        if not neighbours:
            break

        #compute best neighbour node by the lsm_calc

        v = max(neighbours, key=lambda x: (lsm_calc.get_lsm_append_node(x), -min(nx.shortest_path_length(G, q_node, x) for q_node in q), -int(x)))


        current_C = list(C.graph.nodes())
        current_C.append(v)
        C.sequence[i + 1] = v
        lsm_calc.add_node(v)

        community = G.subgraph(current_C).copy()
        C = SubgraphData(
            community,
            lsm_calc.get_lsm(),
            len(community),
            C.sequence
        )
        i += 1

        if C.size >= l:
            if best_lsm < C.lsm_value:
                best_lsm = C.lsm_value
                best_graph = C.graph.copy()

    return best_graph, best_lsm, C