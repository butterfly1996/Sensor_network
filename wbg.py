import numpy as np
from scipy.optimize import linear_sum_assignment
import distance

from sensors_field import Sensors_field
from sensor import Sensor

class WBG(Sensors_field):
    def __init__(self, lenght, height, mode='weak'):
        Sensors_field.__init__(self, lenght, height)

        # sensor ung voi le trai va le phai va ko thuoc sensors_list
        self.s = Sensor(xi=0, name='s')
        self.t = Sensor(xi=lenght, name='t')
        if mode not in ['weak', 'strong']:
            raise ValueError
        self.mode = mode

    def build_WBG(self):
        ## ma tran ke
        ## index 0 la s
        ## index n+1 la t
        ## tu 1 den n la chi so cua cac sensor
        self.adj_matrix = np.full((len(self.sensors_list) + 2, len(self.sensors_list) + 2), np.inf)
        for i, si in enumerate([self.s] + self.sensors_list + [self.t]):
            for j, sj in enumerate([self.s] + self.sensors_list + [self.t]):
                if i < j:
                    if i == 0 and j == len(self.sensors_list) + 1:
                        continue
                    self.adj_matrix[i][j] = self.w(si, sj)
                elif i == j:
                    self.adj_matrix[i][j] = 0
                else:
                    self.adj_matrix[i][j] = self.adj_matrix[j][i]
        self.o_adj_matrix = self.adj_matrix.copy()

    def dw(self, vi, vj):  # weak distance
        if vi == self.s:
            return vj.xL
        if vi == self.t:
            return self.L - vj.xR
        if vj == self.s:
            return vi.xL
        if vj == self.t:
            return self.L - vi.xR
        elif (vi.xL <= vj.xL and vj.xL <= vi.xR) or (vj.xL <= vi.xL and vi.xL <= vj.xR):
            return 0
        elif vj.xL - vi.xR > 0:
            return vj.xL - vi.xR
        else:
            return vi.xL - vj.xR

    def ds(self, vi, vj):  # strong distance
        if vi == self.s or vi == self.t or vj == self.s or vj == self.t:
            return self.dw(vi, vj)
        elif vi.overlap(vj):
            return 0
        else:
            return distance.minimum__sectors_distance(vi, vj)[2]

    def w(self, vi, vj, normalize=True):  # weight
        try:
            lr = vi.lr
        except:
            lr = vj.lr
        func = np.ceil

        if not normalize:
            lr = 1.0
            func = lambda x: x

        if self.mode == 'weak':
            return func(self.dw(vi, vj) / lr)
        else:
            return func(self.ds(vi, vj) / lr)

    def show_matrix(self):
        print(self.adj_matrix)

    def dijkstra(self):
        dist = [np.inf] * self.adj_matrix.shape[0]
        prev = [None] * self.adj_matrix.shape[0]
        Q = []
        if len(self.adj_matrix) == 0:
            return None
        for i in range(self.adj_matrix.shape[0]):
            # dist[i] = np.inf
            # prev[i] = None
            Q.append(i)
        dist[0] = 0
        while len(Q) is not 0:
            u = Q[np.argmin([dist[q] for q in Q])]
            Q.remove(u)
            for i in range(self.adj_matrix.shape[0]):
                if i == u:
                    continue
                if (u == 0 and i == self.adj_matrix.shape[0] - 1) or (u == self.adj_matrix.shape[0] - 1 and i == 0):
                    continue
                # print("debug:: ", u, ":::", i)
                alt = dist[u] + self.adj_matrix[u][i]
                if alt < dist[i]:
                    dist[i] = alt
                    prev[i] = u

        j = self.adj_matrix.shape[0] - 1
        p = [j]
        while prev[j] != None:
            p.append(prev[j])
            j = prev[j]
        p = list(reversed(p))
        return p

    def length(self, p, pso=False):
        res = 0
        if not pso:
            if len(p) == 2 and p[0] == 0 and p[1] == self.adj_matrix.shape[0] - 1:  # direct barrier
                lr = self.sensors_list[0].lr
                return np.ceil(self.L / lr)
        for i in range(len(p) - 1):
            # print ('add', self.adj_matrix[p[i]][p[i + 1]])
            res += self.adj_matrix[p[i]][p[i + 1]]
        return res

    def min_num_mobile_greedy(self, k):
        Pk = []
        q = 0
        while q < k:
            p = self.dijkstra()
            if len(p) < 2:  ## empty path
                # print('%')
                break
            if self.length(p) <= np.ceil(self.L / self.sensors_list[0].lr):
                Pk.append(p)
                q += 1
                for i in range(1, len(p) - 1):
                    self.adj_matrix[p[i], :] = np.inf
                    self.adj_matrix[:, p[i]] = np.inf
            else:
                break
        Nm = 0
        for p in Pk:
            Nm += self.length(p)
        if q < k:
            for i in range(k - q):
                Pk.append([0, len(self.sensors_list) + 1])
            Nm = 0
            for p in Pk:
                Nm += self.length(p)
        return Pk, Nm

    def max_num_barrier_greedy(self, tau):  ## tau: number of deployed mobile sensors
        q = 0
        Pq = []
        while True:
            p = self.dijkstra()
            if len(p) < 2:  ## empty path
                return q + np.floor(
                    (tau - np.sum([self.length(pq) for pq in Pq])) / np.ceil(self.L / self.sensors_list[0].lr)), Pq
            else:
                if self.length(p) <= np.ceil(self.L / self.sensors_list[0].lr):
                    if np.sum([self.length(pq) for pq in Pq]) + self.length(p) <= tau:
                        q += 1
                        Pq.append(p)
                        for i in range(1, len(p) - 1):
                            self.adj_matrix[p[i], :] = np.inf
                            self.adj_matrix[:, p[i]] = np.inf
                    else:
                        if np.sum([self.length(pq) for pq in Pq]) + self.length(p) == tau:
                            return q + 1, Pq
                        else:
                            return q, Pq
                else:
                    return q + np.floor(
                        (tau - np.sum([self.length(pq) for pq in Pq])) / np.ceil(self.L / self.sensors_list[0].lr)), Pq

    def calculate_target_location(self, s1, s2):
        if s1 not in [0, len(self.sensors_list) + 1] and s2 not in [0, len(self.sensors_list) + 1]:
            pa, pb, d = distance.minimum__sectors_distance(self.sensors_list[s1 - 1], self.sensors_list[s2 - 1])
        elif s1 == 0:
            xL, yL = self.sensors_list[s2 - 1].xL, self.sensors_list[s2 - 1].yL
            pa, pb, d = np.array([0, yL]), np.array([xL, yL]), xL
        elif s2 == 0:
            xL, yL = self.sensors_list[s1 - 1].xL, self.sensors_list[s1 - 1].yL
            pa, pb, d = np.array([0, yL]), np.array([xL, yL]), xL
        elif s1 == len(self.sensors_list) + 1:
            xR, yR = self.sensors_list[s2 - 1].xR, self.sensors_list[s2 - 1].yR
            pa, pb, d = np.array([self.L, yR]), np.array([xR, yR]), self.L - xR
        elif s2 == len(self.sensors_list) + 1:
            xR, yR = self.sensors_list[s1 - 1].xR, self.sensors_list[s1 - 1].yR
            pa, pb, d = np.array([self.L, yR]), np.array([xR, yR]), self.L - xR
        if self.o_adj_matrix[s1, s2] > 0:
            dv = d / self.o_adj_matrix[s1, s2]
            phi = distance.arctan(pb[1] - pa[1], pb[0] - pa[0])
            alpha = self.sensors_list[0].alpha
            lr = self.sensors_list[0].lr
            r = self.sensors_list[0].r

            res = []
            if lr == r and 2 * alpha < np.pi:
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(np.array([pa[0] + i * dv * np.cos(phi), pa[1] + i * dv * np.sin(phi), phi]))
            elif lr == 2 * r * np.sin(alpha) and 2 * alpha < np.pi:
                h = np.sqrt(r ** 2 - (lr / 2) ** 2)
                l = np.sqrt(h ** 2 + (dv / 2) ** 2)
                lamda = distance.arctan(2 * h, dv)
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(np.array([pa[0] + i * dv * np.cos(phi) + l * np.cos(phi + lamda),
                                         pa[1] + i * dv * np.sin(phi) + l * np.sin(phi + lamda),
                                         angle(phi + 3 / 2 * np.pi)]))
            elif lr == 2 * r and 2 * alpha > np.pi:
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(
                        np.array([pa[0] + i * dv * np.cos(phi), pa[1] + i * dv * np.sin(phi), angle(phi + np.pi / 2)]))
            return res
        else:
            return []

    def mcbf(self, k, dynamic_sens):  # k la so barierr, dynamic_sens list cac sensor dong duoc trien khai
        tau = len(dynamic_sens)
        Pk, Nm = self.min_num_mobile_greedy(k)
        print ("Nm %d" % Nm)
        locs = []
        for path in Pk:
            print(path)
            for i in range(len(path) - 1):
                locs.extend(self.calculate_target_location(path[i], path[i + 1]))
        if tau < len(locs):
            return (None, None), np.inf
        d = np.zeros((tau, Nm))
        for i, dynamic_sen in enumerate(dynamic_sens):
            for j, loc in enumerate(locs):
                d[i, j] = np.sqrt((dynamic_sen.xi - loc[0]) ** 2 + (dynamic_sen.yi - loc[1]) ** 2)
        row_ind, col_ind = linear_sum_assignment(d)
        min_cost = d[row_ind, col_ind].sum()
        return locs, (row_ind, col_ind), min_cost
        # locs la vi tri cac target, (row_ind, col_ind) la ghep cap giua dynamic sensor den target, min_cost chi phi minimum