import networkx as nx
from networkx.algorithms import approximation as ax
import matplotlib.pyplot as plt


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


def steiner_tree(G, q):
    G1 = G.copy()
    # Add weights to the graph (using max degree)
    for u, v in G1.edges():
        weight = max(G1.degree(u), G1.degree(v))
        G1[u][v]['weight'] = weight

    seed = ax.steinertree.steiner_tree(G1, q, method='mehlhorn')
    return seed


def get_neighbour(G, C):
    neighbours = set()
    for x in G.subgraph(C):
        neighbours.update(G.neighbors(x))
    C = set(C)
    neighbours -= C
    return list(neighbours)

def get_neighbour_include_set(G, C):
    neighbours = set()
    for x in G.subgraph(C):
        neighbours.update(G.neighbors(x))
    return list(neighbours)


def LSM(G, C, t, v=None):
    C_sub = set(C.nodes())
    if v is not None:
        C_sub.add(v)
    C_sub = G.subgraph(C_sub)
    in_degree = C_sub.number_of_edges()
    out_degree = sum(G.degree(u) for u in C_sub.nodes()) - 2 * in_degree
    LM = in_degree / out_degree if out_degree != 0 else 0
    return LM * (1 / pow(len(C_sub),t)), in_degree, out_degree