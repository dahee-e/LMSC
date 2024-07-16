import argparse
import networkx as nx
import os
import time
import sys
import NGA
import MMA

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

parser.add_argument('--algorithm', default="NGA",
                    help='specify algorithm name')

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

start_time = time.time()




if args.algorithm == 'NGA':
    process ,C, C_lsm = NGA.run(G, args.q, args.l, args.h, args.t)
elif args.algorithm == 'MMA':
    process ,C, C_lsm = MMA.run(G, args.q, args.l, args.h, args.t)

run_time = time.time() - start_time

result = list()
if C != None :
    result = list(C)
print("----------------------------------------------------------")
# for comp in result:
#     comp = [int(x) for x in comp]
#     comp = sorted(comp, reverse=False)
#     for u in list(comp):
#         print(u, ' ', end="")
#     print("")
print("----------------------------------------------------------")

with open(output, 'w') as f:
    f.write("seconds" + "\t" + str(run_time) + '\n')
    f.write("nodes" + "\t")
    for node in result:
        f.write(str(node) + " ")
    f.write("\n")

    f.write("lsm\t"+str(C_lsm) + '\n')
    f.write("----------------------------------------------------------" + '\n')
    f.write("----------------------------------------------------------" + '\n')
    f.write("PROCESS : \n")

#process data format
# class SubgraphData:
#     def __init__(self, graph, lsm_value, size):
#         self.graph = graph
#         self.lsm_value = lsm_value
#         self.size = size
    i = 0
    for p in process:
        f.write("["+str(i)+"]\nlsm_value : " + str(p.lsm_value) + "\nsize : " + str(p.size) + "\nnode:")
        for node in p.sequence:
             f.write(str(node) + " ")
        f.write("\n")
        f.write("----------------------------------------------------------" + '\n')
        i += 1
f.close()
