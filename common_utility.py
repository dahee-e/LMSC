import networkx as nx
from networkx.algorithms import approximation as ax
import matplotlib.pyplot as plt


# def steiner_tree(G, q):
#     seed = ax.steinertree.steiner_tree(G, q, method='mehlhorn')
#     # #visualize the graph
#     # nx.draw(seed, with_labels=True)
#     # plt.show()
#     return seed

def steiner_tree(G, q):
    seed = ax.steinertree.steiner_tree(G, q, method='mehlhorn')
    # #visualize the graph
    # nx.draw(seed, with_labels=True)
    # plt.show()
    return seed


def get_neighbour(G, C):
    neighbours = set()
    for x in G.subgraph(C.nodes()):
        neighbours.update(G.neighbors(x))
    neighbours -= set(C.nodes())
    return list(neighbours)


def LSM(G, C, t, v=None):
    C_sub = set(C.nodes())
    if v is not None:
        C_sub.add(v)
    C_sub = G.subgraph(C_sub)
    in_degree = C_sub.number_of_edges()
    out_degree = sum(G.degree(u) for u in C_sub.nodes()) - 2 * in_degree
    LM = in_degree / out_degree if out_degree != 0 else 0
    return LM * (1 / pow(len(C_sub),int(t)))