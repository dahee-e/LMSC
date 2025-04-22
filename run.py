import argparse
import networkx as nx
import os
import time
import sys
import IGA
import SMA


sys.setrecursionlimit(10000)


def get_base(file_path) :
    path = os.path.dirname(file_path)
    return path+"/"



#############################################################################
parser = argparse.ArgumentParser(description='Run algorithm on network')
parser.add_argument('--q',type=str, nargs='*', default=['5','13'],
                    help='a set of query nodes q')

parser.add_argument('--l', type=int, default=1,
                    help='lower bound l')

parser.add_argument('--h', type=int, default=16,
                    help='upper bound h')

parser.add_argument('--t', type=float, default=1,
                    help='threshold tau for local sketch modularity`')

parser.add_argument('--network', default="./dataset/karate/network.dat",
                    help='a folder name containing network.dat')

parser.add_argument('--algorithm', default="SMA",
                    help='specify algorithm name')
parser.add_argument('--weak', type=bool, default=False)
# parser.add_argument('--scalability', type=bool, default=False,
#                     help='for scalability test')






args = parser.parse_args()
print("network ", args.network)
print("algorithm ", args.algorithm)


params  = dict()
params['q'] = args.q
params['l'] = args.l
params['h'] = args.h
params['t'] = args.t

output = get_base(args.network)
output = output + args.algorithm
for key in params.keys() :
    print("params ", key, params[key])
    output = output + "_"+str(key)+"_"+str(params[key])
output = output+".dat"
print("output", output)

#############################################################################
#############################################################################

#############################################################################
# Global Parameter
G = None
C = None
#############################################################################


# read network





G = nx.read_edgelist(args.network)
G.remove_edges_from(nx.selfloop_edges(G))
print("entire graph V=", G.number_of_nodes(), "\tE=", G.number_of_edges())
#############################################################################

if nx.is_connected(G) == False:
    print("graph is not connected")
    sys.exit(0)

start_time = time.time()

if args.algorithm == 'IGA':
    best_graph,best_lsm,C = IGA.run(G, args.q, args.l, args.h, args.t)
elif args.algorithm == 'SMA':
    best_graph,best_lsm,C = SMA.run(G, args.q, args.l, args.h, args.t)


run_time = time.time() - start_time


with open(output, 'w') as f:
    f.write("seconds" + "\t" + str(run_time) + '\n')
    f.write("nodes" + "\t")
    if best_graph is not None:
        for node in best_graph.nodes():
            f.write(str(node) + " ")
        f.write("\n")
        f.write("lsm\t" + str(best_lsm) + '\n')
        for i, node in C.sequence.items():
            f.write(str(i) + "\t" + str(node) + '\n')
    else :
        f.write("Can't find community")
f.close()
