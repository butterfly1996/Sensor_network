import itertools
import multiprocessing
import numpy as np
from igraph import *
from circle_util import Geometry

from wbg import WBG
geom = Geometry()
from utils import *
import global_vars


def fill_mat(params):
    # global global_vars.misfit
    ii, jj = params
    if global_vars.dbg.vid_to_vs[ii].overlap(global_vars.dbg.vid_to_vs[jj]):
        global_vars.misfit[ii, jj] = global_vars.dbg.rotated_angle(global_vars.dbg.vid_to_as[jj], global_vars.dbg.vid_to_vs[jj])

class DBG(WBG):
    def __init__(self, lenght, height):
        WBG.__init__(self, lenght, height, mode='strong')
        self.H = height
        self.L = lenght

    def build_virtual_nodes(self):
        for i in range(len(self.sensors_list)):
            for j in range(i + 1, len(self.sensors_list)):
                u = self.sensors_list[i]
                v = self.sensors_list[j]
                if len(geom.circle_intersection((u.xi, u.yi, u.r), (v.xi, v.yi, v.r))) > 0 and i != j:
                    u.set_virtual_nodes(v)
                    v.set_virtual_nodes(u)
        # for s and t
        for i in range(len(self.sensors_list)):
            u = self.sensors_list[i]
            if u.xi <= u.r:
                u.add_virtual_node(angle(+np.arccos(-u.xi / u.r) - u.alpha))
                u.add_virtual_node(angle(-np.arccos(-u.xi / u.r) + u.alpha))
                u.add_virtual_node(angle(+np.arccos(-u.xi / u.r) + u.alpha))
                u.add_virtual_node(angle(-np.arccos(-u.xi / u.r) - u.alpha))
        for i in range(len(self.sensors_list)):
            u = self.sensors_list[i]
            if u.xi >= self.L - u.r:
                u.add_virtual_node(angle(+np.arccos((self.L - u.xi) / u.r) - u.alpha))
                u.add_virtual_node(angle(-np.arccos((self.L - u.xi) / u.r) + u.alpha))
                u.add_virtual_node(angle(+np.arccos((self.L - u.xi) / u.r) + u.alpha))
                u.add_virtual_node(angle(-np.arccos((self.L - u.xi) / u.r) - u.alpha))
        self.num_virtuals = np.sum([len(s.virtual_nodes) for s in self.sensors_list]) + 2

    def rotated_angle(self, actual_node, virtual_node):
        actual = actual_node.betai
        virtual = virtual_node.betai
        return min(abs(angle(actual - virtual)), 2 * np.pi - abs(angle(actual - virtual)))

    def build_dbg(self):
        num_virtuals = self.num_virtuals

        self.adj_matrix = np.full((num_virtuals, num_virtuals), np.inf)

        print ('nvir %d' % num_virtuals)
        # print ('detail ', [len(s.virtual_nodes) for s in self.sensors_list])

        self.s.vid = self.s.id = 0
        self.t.vid = num_virtuals - 1
        self.t.id = len(self.sensors_list) + 1

        self.vid_to_as = {self.s.vid: self.s, self.t.vid: self.t}
        self.vid_to_vs = {self.s.vid: self.s, self.t.vid: self.t}
        temp = 1
        for i, u in enumerate(self.sensors_list):
            u.id = i + 1
            for ii, uu in enumerate(u.virtual_nodes):
                uu.vid = temp
                self.vid_to_as[uu.vid] = u
                self.vid_to_vs[uu.vid] = uu
                temp += 1

        # print ('s->')
        # s->
        aj = 1
        for j, v in enumerate(self.sensors_list):
            if v.xi <= v.r:
                for jj, vv in enumerate(v.virtual_nodes):
                    if np.equal(super().ds(self.s, vv), 0):
                        global_vars.misfit[0][aj + jj] = self.rotated_angle(v, vv)
            aj += len(v.virtual_nodes)

        # print ('t->')
        # -> t
        ai = 1
        for i, u in enumerate(self.sensors_list):
            if u.xi >= self.L - u.r:
                for ii, uu in enumerate(u.virtual_nodes):
                    if np.equal(super().ds(uu, self.t), 0):
                        global_vars.misfit[ai + ii][num_virtuals - 1] = 0 #self.rotated_angle(u, uu)
            ai += len(u.virtual_nodes)

        pool = multiprocessing.Pool(24)
        # print ('aaaa')
        num_inter = 0
        for i, u in enumerate(self.sensors_list):
            for j, v in enumerate(self.sensors_list):
                if len(geom.circle_intersection((u.xi, u.yi, u.r), (v.xi, v.yi, v.r))) > 0 and u.xi < v.xi and i != j:
                    # print ('i j ', i, j)
                    ii_list = [vn.vid for vn in u.virtual_nodes]
                    jj_list = [vn.vid for vn in v.virtual_nodes]

                    paramlist = list(itertools.product(ii_list, jj_list))
                    # print (len(paramlist))
                    num_inter += 1
                    pool.map(fill_mat, paramlist)
        pool.close()
        self.adj_matrix = global_vars.misfit
        # print ('aaaa')

        '''
        for ii in ii_list:
            for jj in jj_list:
                self.fill_mat([ii, jj])
        '''
        # print ("inter count %d" % num_inter)

    def show_matrix(self):
        print(self.adj_matrix)

    def dfs2(self, s, t):
        # print ('remove')
        temp_adj = self.adj_matrix.copy()
        temp_adj[temp_adj == np.inf] = 0
        # print ('remove')

        g = Graph.Weighted_Adjacency(temp_adj.tolist(), ADJ_DIRECTED)
        # print ('zero length')
        xs, ys = np.where(self.adj_matrix == 0)
        for x, y in zip(xs, ys):
            g.add_edge(x, y)
        # print ('zero length')

        s = set(g.subcomponent(s, mode="out"))
        t = set(g.subcomponent(t, mode="in"))
        return True if len(s.intersection(t)) > 0 else False

    def dijkstra2(self, s, t):
        # print ('remove')
        temp_adj = self.adj_matrix.copy()
        temp_adj[temp_adj == np.inf] = 0
        # print ('remove')

        g = Graph.Weighted_Adjacency(temp_adj.tolist(), ADJ_DIRECTED)
        # print ('zero length')
        xs, ys = np.where(self.adj_matrix == 0)
        for x, y in zip(xs, ys):
            g.add_edge(x, y)
        # print ('zero length')

        temp_w = g.es['weight']
        for x, y in zip(xs, ys):
            temp_w[g.get_eid(x, y)] = 0
        g.es['weight'] = temp_w

        return g.shortest_paths_dijkstra(s, t, 'weight')[0][0]

    def max_flow(self):
        new_adj = np.full((2 * (self.num_virtuals-2) + 2, 2 * (self.num_virtuals-2) + 2), np.inf)
        new_adj[:1, :(self.num_virtuals-2) + 2] = self.adj_matrix[:1, :(self.num_virtuals-2) + 2].copy()
        new_adj[(self.num_virtuals-2) + 2:, 1:(self.num_virtuals-2) + 2] = self.adj_matrix[1:(self.num_virtuals-2) + 1, 1:].copy()
        eye = np.eye((self.num_virtuals - 2))
        eye[eye == 0] = np.inf
        new_adj[1:(self.num_virtuals-2) + 1, (self.num_virtuals-2) + 2:] = eye


        # correct Graph
        # print ('remove')
        temp_adj = new_adj.copy()
        temp_adj[temp_adj == np.inf] = 0
        # print ('remove')

        # print ('zero length')
        temp_adj[new_adj == 0] = 0.001
        # print ('zero length')
        g = Graph.Weighted_Adjacency(temp_adj.tolist(), ADJ_DIRECTED)

        # round 1
        s = 0
        t = self.adj_matrix.shape[0] - 1
        # res1 = g.maxflow(s, t)
        # print ('res 1, ', res1.value)

        # round 1
        res1 = g.maxflow(s, t)
        print ('res 1, ', res1.value)

        # round 2
        g.es['weight'] = res1.flow  # set weight to flow value
        g.delete_edges(np.where(np.array(g.es['weight']) == 0)[0].tolist())  # remove unnecessary edges
        paths = []
        assert len(g.incident(s)) == res1.value, "failed assert"
        for i, e in enumerate(g.incident(s)):
            path = []
            u = g.get_edgelist()[e][1]
            #print (e, u)
            path.append(u)
            while g.get_edgelist()[g.incident(u)[0]][1] != t:
                u = g.get_edgelist()[g.incident(u)[0]][1]
                #print (u)
                path.append(u)
            paths.append(path)
        #print (paths)
        #for p in paths:
            #print ('act paths', [self.vid_to_as[v].id for v in p if v < t])
            #print ('vir paths', [v for v in p if v < t])
        def get_conflict_path_ids(paths, id):
            path = paths[id]
            ids = [True if len(set([self.vid_to_as[v].id for v in p if v<t]).intersection(set([self.vid_to_as[v].id for v in path if v<t]))) > 0 else False
                 for p in paths]
            ids = np.where(np.array(ids)==True)[0]
            # for p in paths:
            #     print ('all paths', [self.vid_to_as[v].id for v in p if v < t])
            # print  ('current path', [self.vid_to_as[v].id for v in path if v < t])
            # print  ('paths ', paths)
            # print  ('ids ', ids.tolist())
            # print  ('id ', id)

            # ids.tolist().remove(id)

            return ids
        new_paths = []
        while len(paths) > 0:
            argmin = np.argmin([len(get_conflict_path_ids(paths, i)) for i in range(len(paths))])
            new_paths.append(argmin)
            cids = get_conflict_path_ids(paths, argmin)
            for cid in sorted(cids,reverse=True):
                del paths[cid]
        return len(new_paths)