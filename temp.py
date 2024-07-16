#finding a community for a give source v in a graph G
#using the Luo algorithm
#input: G, v
#output: a community C
#intialize C = {v}
#create a new subgraph S with v
#do
#create a new neighbor set N with adjacent vertices of v.

    #do
    #create a new list Q to store new adding vertices
    #sort the vertex in N based on their degree increment
    #for each vertex u in N
    #   compute luo modularity \delta M = \frac{community(S add u) inside node}{community(S add u) connected out edge} - \frac{community(S) in edge}{community(S) out edge}
    #   if \delta M > 0
    #       add u to S
    #       add u list Q
    #   endif
    #endfor
    #while Q is empty

    #do
    #create a new list deleteQ
    #for each vertex u in V_s
    #   compute luo modularity delta M = \frac{community(S remove u) inside node}{community(S remove u) connected out edge} - \frac{community(S) in edge}{community(S) connected out edge}
    #   if \delta M > 0 and removing u from S does not disconnect S
    #       remove v from S
    #       add v to deleteQ
    #   endif
    #endfor
    #while deleteQ is empty

    #create a new neighbor set newN of S
    #while newN is the same as N
    #if modularity of S > 1 and S contains v
    #   return S
    #else
    #   print out "no community found"
#/*end*/
def luo(G, q):
    C = set()
    C.add(q)
    S = nx.Graph()
    S.add_node(q)
    N = set()
    N.update(G.neighbors(q))
    while len(N) > 0:
        newN = set()
        for i in S.nodes():
            newN.update(G.neighbors(i))
        newN -= C
        Q = []
        for u in N:
            deltaM = get_modularity_append(G, S, u) - get_modularity(G, S)
            if deltaM > 0:
                S.add_node(u)
                for v in S.nodes():
                    if G.has_edge(u, v):
                        S.add_edge(u, v)
                C.add(u)
                Q.append(u)

        deleteQ = []
        for u in Q:
            deltaM = get_modularity_remove(G, S, u) - get_modularity(G, S)
            S_temp = S.copy()
            S_temp.remove_node(u)
            if deltaM > 0 and nx.is_connected(S_temp):
                S.remove_node(u)
                C.remove(u)
                deleteQ.append(u)

        N = set()
        for i in S.nodes():
            N.update(G.neighbors(i))
        N -= C
        if N == newN:
            break
    if get_modularity(G, S) > 0 and S.has_node(v):
        return sorted(C)
    else:
        print("no community found")
        return None

import networkx as nx
import time


#modularity = \frac{community indegree}{community outdegree}
def get_modularity(G, S):
    if len(S) <= 1:
        return 0
    in_degree = S.number_of_edges()
    out_degree = 0
    for u in S.nodes():
        out_degree += G.degree(u)
    out_degree -= in_degree*2
    return in_degree / out_degree
def get_modularity_append(G, S1, u):
    #append_indegree is the number of edges in graph G that connect u to S
    S_temp = S1.copy()
    for v in S1.nodes():
        if G.has_edge(u, v):
            S_temp.add_edge(u, v)
    S1 = S_temp
    in_degree = S1.number_of_edges()
    out_degree = 0
    for u in S1.nodes():
        out_degree += G.degree(u)
    out_degree -= in_degree*2
    return in_degree / out_degree
def get_modularity_remove(G, S, u):
    S1 = S.copy()
    S1.remove_node(u)
    if len(S1) <= 1:
        return 0
    in_degree = S1.number_of_edges()
    out_degree = 0
    for u in S1.nodes():
        out_degree += G.degree(u)
    out_degree -= in_degree*2
    return in_degree / out_degree





#read the graph from the file
fileName = "./dataset/network.dat"
G = nx.read_edgelist(fileName)
print("V=", G.number_of_nodes(), "\tE=", G.number_of_edges())
G.remove_edges_from(nx.selfloop_edges(G))
start_time = time.time()
#query = input("Enter a node: ")
query = "8"
print("The community of node", query, "is", luo(G, query))
print("Running time:", time.time() - start_time, "seconds")

