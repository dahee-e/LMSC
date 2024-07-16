import networkx as nx
from networkx.algorithms import approximation as ax
import matplotlib.pyplot as plt


def steiner_tree(G, q):
    seed = ax.steinertree.steiner_tree(G, q, method='mehlhorn')
    print("seed V=", seed.number_of_nodes(), "\tE=", seed.number_of_edges())

    '''visualize the graph
    # nx.draw(seed, with_labels=True)
    # plt.show()'''

    return seed


# def LSM(G, C, t):
#     in_degree = C.number_of_edges()
#     out_degree = 0
#     for u in C.nodes():
#         out_degree += G.degree(u)
#     out_degree -= in_degree * 2
#     LM = in_degree / out_degree
#
#     return LM * (1 / pow(len(C), t))
#
#
#

def LSM(G, C, t, v=None):
    C = nx.subgraph(G, set(C))
    if v is not None:
        C = nx.subgraph(G, set(C).union({v}))
    in_degree = C.number_of_edges()
    out_degree = sum(G.degree(u) for u in C.nodes()) - 2 * in_degree
    LM = in_degree / out_degree if out_degree != 0 else 0
    return LM * (1 / pow(len(C),int(t)))






def run(G, q, l, h, t):
    C = list()
    C_LSM = list()
    C_size = list()
    i = 0
    if len(q) == 1:
        C.append(G.subgraph(q))
    else:
        C.append(steiner_tree(G, q))

    C_LSM.append(LSM(G, C[i], t))
    C_size.append(len(C[i]))

    while len(C[i]) < h:
        neighbors = set()
        for x in C[i].nodes():
            neighbors.update(G.neighbors(x))
        neighbors -= set(C[i].nodes())

        v = max(neighbors, key=lambda x: LSM(G, C[i], x, t))
        C.append(C[i].copy())
        i += 1
        C[i].add_node(v)
        C_LSM.append(LSM(G, C[i], t))
        C_size.append(len(C[i]))



    # best_C is sorted by C_LSM and best_C satisfy [l,h] size, and return best one
    sort_C = sorted(zip(C, C_LSM, C_size), key=lambda x: x[1], reverse=True)
    size_C= [x for x in sort_C if l <= x[2] <= h]
    best_C = size_C[0][0] if size_C else None
    return best_C