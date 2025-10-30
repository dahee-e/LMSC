import networkx as nx
from networkx.algorithms import approximation as ax
import matplotlib.pyplot as plt


class LSMCalculator:
    def __init__(self, G, initial_community, t):
        self.G = G
        self.t = t
        self.community = set(initial_community)

        # 각 노드별 통계 저장 (커뮤니티 + 이웃만)
        self.node_stats = {}

        # 초기화
        self._initialize_node_stats()


        self.l_C = sum(self.node_stats[v]['in'] for v in self.community) // 2
        self.o_C = sum(self.node_stats[v]['out'] for v in self.community)

    def _initialize_node_stats(self):
        """커뮤니티와 그 이웃들만 초기화"""
        # 1단계: 커뮤니티 노드들 초기화
        for node in self.community:
            neighbour = set(self.G.neighbors(node))
            degree = len(neighbour)
            in_count = len(neighbour & self.community)
            self.node_stats[node] = {
                'in': in_count,
                'out': degree - in_count
            }
            n = set(neighbour) - self.community
            for v in n:
                v_neighbour = set(self.G.neighbors(v))
                v_degree = len(v_neighbour)
                v_in_count = len(v_neighbour & self.community)
                self.node_stats[v] = {
                    'in': v_in_count,
                    'out': v_degree - v_in_count
                }
    # def _initialize_node_stats(self):
    #     """커뮤니티와 그 이웃들만 초기화"""
    #     # 1단계: 커뮤니티 노드들 초기화
    #     for node in self.community:
    #         degree = self.G.degree(node)
    #         self.node_stats[node] = {
    #             'degree': degree,
    #             'in': 0,
    #             'out': degree
    #         }
    #
    #     # 2단계: 커뮤니티 노드들의 이웃 초기화 및 in/out 계산
    #     for node in self.community:
    #         for neighbor in self.G.neighbors(node):
    #             if neighbor in self.community:
    #                 # 커뮤니티 내부 연결
    #                 self.node_stats[node]['in'] += 1
    #                 self.node_stats[node]['out'] -= 1
    #             else:
    #                 # 커뮤니티 외부 이웃 - 없으면 초기화
    #                 if neighbor not in self.node_stats:
    #                     neighbor_degree = self.G.degree(neighbor)
    #                     self.node_stats[neighbor] = {
    #                         'degree': neighbor_degree,
    #                         'in': 0,
    #                         'out': neighbor_degree
    #                     }
    #                 # 외부 이웃 업데이트
    #                 self.node_stats[neighbor]['in'] += 1
    #                 self.node_stats[neighbor]['out'] -= 1

    def get_lsm(self):
        if self.o_C == 0:
            return 0
        return (self.l_C / self.o_C) * (1 / pow(len(self.community), self.t))

    def get_lsm_append_node(self, v):
        """노드 v를 추가했을 때의 LSM 계산"""
        if v in self.community:
            return self.get_lsm()


        w_v_C = self.node_stats[v]['in']
        new_l_C = self.l_C + w_v_C
        new_o_C = self.o_C + self.node_stats[v]['out'] - w_v_C
        new_size = len(self.community) + 1

        if new_o_C == 0:
            return 0
        else:
            return (new_l_C / new_o_C) * (1 / pow(new_size, self.t))

    def add_node(self, v):
        """노드 v를 커뮤니티에 추가하고 새로운 이웃들 초기화"""
        if v in self.community:
            return

        # v가 아직 초기화되지 않았다면 초기화
        if v not in self.node_stats:
            v_degree = self.G.degree(v)
            v_in = sum(1 for neighbor in self.G.neighbors(v)
                      if neighbor in self.community)
            self.node_stats[v] = {
                'in': v_in,
                'out': v_degree - v_in
            }

        w_v_C = self.node_stats[v]['in']

        # 전체 통계 업데이트
        self.l_C += w_v_C
        self.o_C += self.node_stats[v]['out'] - w_v_C

        # v를 커뮤니티에 추가
        self.community.add(v)



        for neighbor in self.G.neighbors(v):
            if neighbor != v and neighbor not in self.community:
                # neighbor가 없으면 초기화
                if neighbor not in self.node_stats:
                    neighbor_degree = self.G.degree(neighbor)
                    neighbor_in = sum(1 for n in self.G.neighbors(neighbor)
                                     if n in self.community)
                    self.node_stats[neighbor] = {
                        'in': neighbor_in,
                        'out': neighbor_degree - neighbor_in
                    }
                else:
                    # 이미 있으면 업데이트만
                    self.node_stats[neighbor]['in'] += 1
                    self.node_stats[neighbor]['out'] -= 1

    def get_stats(self):
        return self.l_C, self.o_C, len(self.community)

    def get_node_stats(self, v):
        return self.node_stats.get(v, None)

    def append_node_stats(self, v):
        if v not in self.node_stats:
            v_degree = self.G.degree(v)
            v_in = sum(1 for neighbor in self.G.neighbors(v)
                      if neighbor in self.community)
            self.node_stats[v] = {
                'in': v_in,
                'out': v_degree - v_in
            }

    def add_node_set(self, v_set):
        for v in v_set:
            self.add_node(v)


class CLSMCalculator:
    # CLSM = (l_Z + l_CZ) / (o_Z - l_CZ)
    def __init__(self, G, C, Z):
        self.G = G
        self.C = set(C) # current community
        self.Z = set(Z)

        self.node_stats = {} # 각 노드별 통계 저장 (체인 + 체인의 이웃)  chian과의 연결, 외부와의 연결, C와의 연결

        self._initialize_node_stats()
        self.l_Z = sum(self.node_stats[v]['inZ'] for v in self.Z) // 2
        self.l_CZ = sum(self.node_stats[v]['inC'] for v in self.Z)
        self.o_Z = sum(self.node_stats[v]['out'] for v in self.Z)
    def _initialize_node_stats(self):
        """체인과 그 이웃들만 초기화"""
        # 1단계: 체인 노드들 초기화
        for node in self.Z:
            neighbour = set(self.G.neighbors(node))
            degree = len(neighbour)
            inZ_count = len(neighbour & self.Z)
            inC_count = len(neighbour & self.C)
            self.node_stats[node] = {
                'inZ': inZ_count,
                'inC': inC_count,
                'out': degree - inZ_count
            }
            n = set(neighbour) - self.Z
            for v in n:
                if v in self.C:
                    continue
                v_neighbour = set(self.G.neighbors(v))
                v_degree = len(v_neighbour)
                v_inZ_count = len(v_neighbour & self.Z)
                v_inC_count = len(v_neighbour & self.C)
                self.node_stats[v] = {
                    'inZ': v_inZ_count,
                    'inC': v_inC_count,
                    'out': v_degree - v_inZ_count
                }

    def get_clsm(self):
        if self.o_Z - self.l_CZ == 0:
            return float('inf')
        return (self.l_Z + self.l_CZ) / (self.o_Z - self.l_CZ)
    def get_clsm_append_node(self, v):
        """노드 v를 추가했을 때의 CLSM 계산"""

        w_v_Z = self.node_stats[v]['inZ']
        w_v_C = self.node_stats[v]['inC']
        w_v_O = self.node_stats[v]['out']
        new_l_Z = self.l_Z + w_v_Z
        new_l_CZ = self.l_CZ + w_v_C
        new_o_Z = self.o_Z + w_v_O - w_v_Z

        if new_o_Z - new_l_CZ == 0:
            return float('inf')
        else:
            return (new_l_Z + new_l_CZ) / (new_o_Z - new_l_CZ)

    def add_node(self, v):
        """노드 v를 체인에 추가하고 새로운 이웃들 초기화"""
        if v in self.Z:
            return

        # v가 아직 초기화되지 않았다면 초기화
        if v not in self.node_stats:
            v_neighbour = set(self.G.neighbors(v))
            v_degree = len(v_neighbour)
            v_inZ_count = len(v_neighbour & self.Z)
            v_inC_count = len(v_neighbour & self.C)
            self.node_stats[v] = {
                'inZ': v_inZ_count,
                'inC': v_inC_count,
                'out': v_degree - v_inZ_count
            }

        w_v_Z = self.node_stats[v]['inZ']
        w_v_C = self.node_stats[v]['inC']
        w_v_O = self.node_stats[v]['out']

        # 전체 통계 업데이트
        self.l_Z += w_v_Z
        self.l_CZ += w_v_C
        self.o_Z += w_v_O - w_v_Z

        # v를 체인에 추가
        self.Z.add(v)

        for neighbor in self.G.neighbors(v):
            if neighbor != v and neighbor not in self.Z and neighbor not in self.C:
                # neighbor가 없으면 초기화
                if neighbor not in self.node_stats:
                    neighbor_nbr = set(self.G.neighbors(neighbor))
                    neighbor_degree = len(neighbor_nbr)
                    neighbor_inZ = len(neighbor_nbr & self.Z)
                    neighbor_inC = len(neighbor_nbr & self.C)
                    self.node_stats[neighbor] = {
                        'inZ': neighbor_inZ,
                        'inC': neighbor_inC,
                        'out': neighbor_degree - neighbor_inZ
                    }
                else:
                    self.node_stats[neighbor]['inZ'] += 1
                    self.node_stats[neighbor]['out'] -= 1



    def get_stats(self):
        return self.l_Z, self.l_CZ, self.o_Z

    def get_node_stats(self, v):
        return self.node_stats.get(v, None)

    def append_node_stats(self, v):
        if v not in self.node_stats:
            v_neighbour = set(self.G.neighbors(v))
            v_degree = len(v_neighbour)
            v_inZ_count = len(v_neighbour & self.Z)
            v_inC_count = len(v_neighbour & self.C)
            self.node_stats[v] = {
                'inZ': v_inZ_count,
                'inC': v_inC_count,
                'out': v_degree - v_inZ_count
            }






def CLSM(G, C, T, v):
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


def steiner_tree(G, q,type):
    G1 = G.copy()
    if type == 'unit':
        for u, v in G1.edges():
            G1[u][v]['weight'] = 1
    elif type == 'min':
        for u, v in G1.edges():
            weight = min(G1.degree(u), G1.degree(v))
            G1[u][v]['weight'] = weight
    elif type == 'max':
        for u, v in G1.edges():
            weight = max(G1.degree(u), G1.degree(v))
            G1[u][v]['weight'] = weight
    elif type == 'sum':
        for u, v in G1.edges():
            weight = G1.degree(u) + G1.degree(v)
            G1[u][v]['weight'] = weight
    elif type =='inverse':
        for u, v in G1.edges():
            weight = max(G1.degree(u), G1.degree(v))
            G1[u][v]['weight'] = 1.0 / weight
    seed = ax.steinertree.steiner_tree(G1, q, method='mehlhorn')
    return seed


def get_neighbour(G, C):
    neighbours = set()
    for x in C:
        neighbours.update(G.neighbors(x))
    neighbours -= set(C)
    return list(neighbours)


def get_neighbour_include_set(G, C):
    neighbours = set()
    for x in C:
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
    return LM * (1 / pow(len(C_sub), t)), in_degree, out_degree