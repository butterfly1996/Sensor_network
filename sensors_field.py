from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random
from matplotlib import collections  as mc
import matplotlib.patches as patches
random.seed(149)

# import distance
def angle(value):
    # chuan hoa gia tri goc tu -pi den pi
    return (value+np.pi) % (2*np.pi) - np.pi

class Sensor():
    def __init__(self, xi=None, yi=None, betai=None, r=None, alpha=None):
        if xi == None or yi == None or betai == None or r == None or alpha == None:
            return
        self.xi = xi
        self.yi = yi
        betai = angle(betai)
        alpha = angle(alpha)
        self.betai = betai
        self.alpha = alpha
        self.r = r
        self.alpha = alpha
        self.lr = 2*r if self.alpha>np.pi/2 else max(r, 2*r*np.sin(self.alpha))

        if -np.pi/2 <= angle(self.betai - self.alpha) <= np.pi/2 and -np.pi/2 <= angle(self.betai + self.alpha) <= np.pi/2:
            self.xL = self.xi
        elif angle(self.betai - self.alpha)*angle(self.betai + self.alpha) < 0 and np.abs(angle(self.betai - self.alpha)-angle(self.betai + self.alpha)) > np.pi:
            self.xL = self.xi - self.r
        else:
            self.xL = min(self.xi + r*np.cos(angle(self.betai-self.alpha)), self.xi + r*np.cos(angle(self.betai+self.alpha)))

        if (-np.pi <= angle(self.betai-self.alpha) <= -np.pi/2 or np.pi/2 <= angle(self.betai-self.alpha) <= np.pi)\
                and (-np.pi <= angle(self.betai+self.alpha) <= -np.pi/2 or np.pi/2 <= angle(self.betai+self.alpha) <= np.pi):
            self.xR = self.xi
        elif angle(self.betai - self.alpha)*angle(self.betai + self.alpha) < 0 and np.abs(angle(self.betai - self.alpha)-angle(self.betai + self.alpha)) < np.pi:
            self.xR = self.xi + self.r
        else:
            self.xR = max(self.xi + r * np.cos(angle(self.betai - self.alpha)),
                          self.xi + r * np.cos(angle(self.betai + self.alpha)))

        # print(self.xL, self.xR)

        # self.lr = 2 * self.r if alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(alpha))
    def overlap(self, s2):
        return True if distance.minimum__sectors_distance(self, s2)[2] == 0 else False
class Sensors_field():
    def __init__(self, lenght, height):
        self.L=lenght
        self.H=height
        # self.n=num_station
        # self.to = num_mobile
        # # alpha la nua goc nhin cua sensor
        # self.alpha = alpha
        #
        # self.lr = 2 * self.r if self.alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(self.alpha))
        # self.xR = max(self.xi, self.xi+self.r*np.cos(self.betai-self.alpha), self.xi+self.r*np.cos(self.betai+self.alpha),
        # self.xi+self.r if -self.alpha <= -self.betai and -self.betai <= self.alpha else 0)
        # self.xL = min(self.xi, self.xi+self.r*np.cos(self.betai-self.alpha), self.xi+self.r*np.cos(self.betai+self.alpha),
        # self.xi+self.r if -self.alpha <= -self.betai and -self.betai <= self.alpha else self.xi+9*self.r)
        self.pointslist = []
        self.sensors_list = []
    def create_sensors_randomly(self, num_sensor = 100, r=3, alpha=60):
        for i in range(0, num_sensor):
            sensor = Sensor(xi=random.uniform(0, self.L), yi = random.uniform(0, self.H), betai= random.uniform(0, 360), r=r, alpha=alpha)
            self.sensors_list.append(sensor)
    def field_show(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplot()
        for i in range(0, len(self.sensors_list)):
            sens = self.sensors_list[i]
        # for sens in self.sensors_list:
        #     fov = Wedge((sens.xi, sens.yi), sens.r, (sens.betai-sens.alpha)/np.pi*180, (sens.betai+sens.alpha)/np.pi*180, color=np.array([random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]), alpha=0.5, label=str(i))

            # ax.add_artist(fov)
            patch = patches.Wedge((sens.xi, sens.yi), sens.r, (sens.betai - sens.alpha) / np.pi * 180,
                          (sens.betai + sens.alpha) / np.pi * 180,
                          color=np.array([random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)]), alpha=0.5,
                          label=str(i + 1))
            ax.add_patch(patch)
        lines = []
        # show = True
        # for cupple_point in self.pointslist:
        #      if cupple_point == None:
        #          show = False
        #          break
        # if show == True:
        for cupple_point in self.pointslist:
            try:
                lc = mc.LineCollection(np.array([cupple_point]), linewidths=2)
                ax.add_collection(lc)
            except ValueError:
                pass
        plt.xlim(xmax = self.L, xmin=0)
        plt.ylim(ymax = self.H, ymin = 0)
        plt.legend()  # <--- here
        plt.show()
        pass
    def add_sensor(self, sensor):
        self.sensors_list.append(sensor)
    def add_dis_2_points(self, points):
        self.pointslist.append(points)
        # print(self.pointslist)
        
class WBG(Sensors_field):
    def __init__(self, lenght, height, mode='weak'):
        Sensors_field.__init__(self, lenght, height)

        # sensor ung voi le trai va le phai va ko thuoc sensors_list
        self.s = Sensor()
        self.t = Sensor()
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
    def dw(self, vi, vj):# weak distance
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
        elif vj.xL-vi.xR > 0:
            return vj.xL-vi.xR
        else:
            return vi.xL-vj.xR
    def ds(self, vi, vj):# strong distance
        if vi == self.s or vi == self.t or vj == self.s or vj == self.t:
            return self.dw(vi, vj)
        elif vi.overlap(vj):
            return 0
        else:
            return distance.minimum__sectors_distance(vi, vj)[2]
    def w(self, vi, vj):# weight
        try:
            lr = vi.lr
        except:
            lr = vj.lr
        if self.mode == 'weak':

            return np.ceil(self.dw(vi, vj)/lr)
        else:
            # print("***")
            # print(self.ds(vi, vj)/lr)
            return np.ceil(self.ds(vi, vj)/lr)
    def show_matrix(self):
        print(self.adj_matrix)
    def dijkstra(self):
        dist = [np.inf]*(len(self.sensors_list)+2)
        prev = [None]*(len(self.sensors_list)+2)
        Q = []
        if len(self.adj_matrix)==0:
            return None
        for i, si in enumerate([self.s] + self.sensors_list + [self.t]):
            # dist[i] = np.inf
            # prev[i] = None
            Q.append(i)
        dist[0] = 0
        while len(Q) is not 0:
            u = Q[np.argmin([dist[q] for q in Q])]
            Q.remove(u)
            for i, _ in enumerate([self.s] + self.sensors_list + [self.t]):
                if i == u:
                    continue
                if (u == 0 and i == len(self.sensors_list)+1) or (u == len(self.sensors_list)+1 and i == 0):
                    continue
                # print("debug:: ", u, ":::", i)
                alt = dist[u]+self.adj_matrix[u][i]
                if alt < dist[i]:
                    dist[i] = alt
                    prev[i] = u

        j = len(self.sensors_list)+1
        p = [j]
        while prev[j] != None:
            p.append(prev[j])
            j = prev[j]
        p = list(reversed(p))
        return p
    def length(self, p):
        res = 0
        if len(p) == 2 and p[0] == 0 and p[1] == len(self.sensors_list)+1: # direct barrier
            lr = self.sensors_list[0].lr
            return np.ceil(self.L/lr)
        for i in range(len(p) - 1):
            res += self.o_adj_matrix[p[i]][p[i + 1]]
        return res
    def min_num_mobile_greedy(self, k):
        Pk = []
        q=0
        while q<k:
            p = self.dijkstra()
            if len(p)<2: ## empty path
                print('%')
                break
            if self.length(p)<=np.ceil(self.L/self.sensors_list[0].lr):
                Pk.append(p)
                q+=1
                for i in range(1, len(p)-1):
                    self.adj_matrix[p[i], :] = np.inf
                    self.adj_matrix[:, p[i]] = np.inf
            else:
                break
        Nm = 0
        for p in Pk:
            Nm += self.length(p)
        if q < k:
            for i in range(k-q):
                Pk.append([0, len(self.sensors_list)+1])
            Nm = 0
            for p in Pk:
                Nm += self.length(p)
        return Pk, Nm
    def max_num_barrier_greedy(self, tau): ## tau: number of deployed mobile sensors
        q = 0
        Pq = []
        while True:
            p = self.dijkstra()
            if len(p) < 2: ## empty path
                return q+np.floor((tau-np.sum([self.length(pq) for pq in Pq]))/np.ceil(self.L/self.sensors_list[0].lr)), Pq
            else:
                if self.length(p) <= np.ceil(self.L/self.sensors_list[0].lr):
                    if np.sum([self.length(pq) for pq in Pq]) + self.length(p) <= tau:
                        q += 1
                        Pq.append(p)
                        for i in range(1, len(p) - 1):
                            self.adj_matrix[p[i], :] = np.inf
                            self.adj_matrix[:, p[i]] = np.inf
                    else:
                        if np.sum([self.length(pq) for pq in Pq]) + self.length(p) == tau:
                            return q+1, Pq
                        else:
                            return q, Pq
                else:
                    return q + np.floor((tau - np.sum([self.length(pq) for pq in Pq])) / np.ceil(self.L / self.sensors_list[0].lr)), Pq
    def calculate_target_location(self, s1, s2):
        pa, pb, d = distance.minimum__sectors_distance(self.sensors_list[s1], self.sensors_list[s2])
        if d > 0:
            dv = d/self.o_adj_matrix[s1, s2]
            phi = distance.arctan(pb[1]-pa[1], pb[0]-pa[0])
            alpha = self.sensors_list[0].alpha
            lr = self.sensors_list[0].lr
            r = self.sensors_list[0].r

            res = []
            if lr == r and 2*alpha<np.pi:
                for i in range(self.o_adj_matrix):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi), pa[1]+i*dv*np.sin(phi), phi]))
            elif lr == 2*r*np.sin(alpha) and 2*alpha<np.pi:
                h = np.sqrt(r**2-(lr/2)**2)
                l = np.sqrt(h**2+(dv/2)**2)
                lamda = distance.arctan(2*h, dv)
                for i in range(self.o_adj_matrix):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi)+l*np.cos(phi+lamda), pa[1]+i*dv*np.sin(phi)+l*np.sin(phi+lamda), angle(phi+3/2*np.pi)]))
            elif lr == 2*r and 2*alpha>np.pi:
                for i in range(self.o_adj_matrix):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi), pa[1]+i*dv*np.sin(phi), angle(phi+np.pi/2)]))
            return res
        else:
            return []

    def mcbf(self, k, dynamic_sens): # k la so barierr, dynamic_sens list cac sensor dong duoc trien khai
        tau = len(dynamic_sens)
        Pk, Nm = self.min_num_mobile_greedy(k)
        locs = []
        for path in Pk:
            for i in range(len(path)-1):
                locs.extend(self.calculate_target_location(path[i], path[i+1]))
        d = np.zeros((tau, Nm))
        for i, dynamic_sen in enumerate(dynamic_sens):
            for j, loc in enumerate(locs):
                d[i,j] = np.sqrt((dynamic_sen.xi-loc[0])**2+(dynamic_sen.yi-loc[1])**2)
        row_ind, col_ind = linear_sum_assignment(d)
        return d[row_ind, col_ind].sum() # chi phi minimum





if __name__ == '__main__':
    import  distance
    wbg = WBG(lenght=10, height=10, mode='strong')
    # sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    #s1 = Sensor(3, 3, np.pi / 5, 2, np.pi / 4)
    #s2 = Sensor(8, 3, 5*np.pi / 6, 2, np.pi / 4)
    wbg.create_sensors_randomly(20)
    wbg.build_WBG()
    wbg.show_matrix()

    # Pk, Nm = wbg.min_num_mobile_greedy(3)
    # print("######################################################")
    # print(Pk)
    # print(Nm)
    # ## debug
    # for path in Pk:
    #     if len(path) > 2:
    #         for i in range(1, len(path)-2):
    #             res = distance.minimum__sectors_distance(wbg.sensors_list[path[i]-1], wbg.sensors_list[path[i+1]-1])
    #             wbg.add_dis_2_points([res[0], res[1]])

    print("######################################################")
    Nb, Pk = wbg.max_num_barrier_greedy(8)
    print(Nb)
    print(Pk)
    ## debug
    for path in Pk:
        if len(path) > 2:
            for i in range(1, len(path)-2):
                res = distance.minimum__sectors_distance(wbg.sensors_list[path[i]-1], wbg.sensors_list[path[i+1]-1])
                wbg.add_dis_2_points([res[0], res[1]])
    wbg.field_show()