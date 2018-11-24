from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random
from matplotlib import collections  as mc
import matplotlib.patches as patches
from pso import *
import codecs
import logging
from pyswarm import pso

logging.basicConfig(filename='tunning.log',level=logging.INFO)

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
        if alpha < 0:
            alpha = -alpha
        self.betai = betai
        self.alpha = alpha
        self.r = r
        self.lr = 2*r if self.alpha>=np.pi/2 else max(r, 2*r*np.sin(self.alpha))

        # leftmost point
        if -np.pi/2 <= angle(self.betai - self.alpha) <= np.pi/2 and -np.pi/2 <= angle(self.betai + self.alpha) <= np.pi/2:
            self.xL = self.xi
            self.yL = self.yi
        elif angle(self.betai - self.alpha)*angle(self.betai + self.alpha) < 0 and np.abs(angle(self.betai - self.alpha)-angle(self.betai + self.alpha)) > np.pi:
            self.xL = self.xi - self.r
            self.yL = self.yi
        else:
            self.xL = min(self.xi + r*np.cos(angle(self.betai-self.alpha)), self.xi + r*np.cos(angle(self.betai+self.alpha)))
            argmin = np.argmin([self.xi + r*np.cos(angle(self.betai-self.alpha)), self.xi + r*np.cos(angle(self.betai+self.alpha))])
            if argmin == 0:
                self.yL = self.yi + r*np.sin(angle(self.betai-self.alpha))
            else:
                self.yL = self.yi + r * np.sin(angle(self.betai+self.alpha))

        # rightmost point
        if (-np.pi <= angle(self.betai-self.alpha) <= -np.pi/2 or np.pi/2 <= angle(self.betai-self.alpha) <= np.pi)\
                and (-np.pi <= angle(self.betai+self.alpha) <= -np.pi/2 or np.pi/2 <= angle(self.betai+self.alpha) <= np.pi):
            self.xR = self.xi
            self.yR = self.yi + r * np.sin(angle(self.betai - self.alpha))
        elif angle(self.betai - self.alpha)*angle(self.betai + self.alpha) < 0 and np.abs(angle(self.betai - self.alpha)-angle(self.betai + self.alpha)) < np.pi:
            self.xR = self.xi + self.r
            self.yR = self.yi + r * np.sin(angle(self.betai - self.alpha))
        else:
            self.xR = max(self.xi + r * np.cos(angle(self.betai - self.alpha)),
                          self.xi + r * np.cos(angle(self.betai + self.alpha)))
            argmax = np.argmax([self.xi + r * np.cos(angle(self.betai - self.alpha)),
                                self.xi + r * np.cos(angle(self.betai + self.alpha))])
            if argmax == 0:
                self.yR = self.yi + r * np.sin(angle(self.betai - self.alpha))
            else:
                self.yR = self.yi + r * np.sin(angle(self.betai + self.alpha))



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
        self.targets = []
        self.dynamics = []
        self.sensors_list = []
        self.mobile_sensors_list = []
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

        x = [p[0] for p in self.dynamics]
        y = [p[1] for p in self.dynamics]
        plt.scatter(x, y, color='red')

        x = [p[0] for p in self.targets]
        y = [p[1] for p in self.targets]
        plt.scatter(x, y, color='green')

        plt.legend()  # <--- here
        plt.show()
        pass
    def add_sensor(self, sensor):
        self.sensors_list.append(sensor)
    def add_dis_2_points(self, points):
        self.pointslist.append(points)
        # print(self.pointslist)
    def add_target(self, point):
        self.targets.append(point)
    def add_dynamic(self, point):
        self.dynamics.append(point)
    def creat_sensors_to_text(self, filename):
        file = codecs.open(filename, mode="w", encoding="utf-8")
        for sensor in self.sensors_list:
            file.write("S"+"\t"+str(sensor.xi)+"\t"+str(sensor.yi)+"\t"+str(sensor.betai)+"\t"+str(sensor.r)+"\t"+str(sensor.alpha)+"\n")
        for sensor in self.mobile_sensors_list:
            file.write("M"+"\t"+str(sensor.xi)+"\t"+str(sensor.yi)+"\t"+str(sensor.betai)+"\t"+str(sensor.r)+"\t"+str(sensor.alpha)+"\n")
    def load_sensors_from_txt(self, filename):
        file = codecs.open(filename, mode="r", encoding="utf-8")
        for line in file.readlines():
            try:
                field_array = line.split("\t")
                if field_array[0] == "S":
                    s = Sensor(xi=float(field_array[1]), yi=float(field_array[2]), betai=float(field_array[3]), alpha=float(field_array[5]), r = float(field_array[4]))
                    self.sensors_list.append(s)
            except:
                print("error")
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
    def length(self, p, pso=False):
        res = 0
        if not pso:
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
        if s1 not in [0, len(self.sensors_list)+1] and s2 not in [0, len(self.sensors_list)+1]:
            pa, pb, d = distance.minimum__sectors_distance(self.sensors_list[s1-1], self.sensors_list[s2-1])
        elif s1 == 0:
            xL, yL = self.sensors_list[s2 - 1].xL, self.sensors_list[s2 - 1].yL
            pa, pb, d = np.array([0, yL]), np.array([xL, yL]), xL
        elif s2 == 0:
            xL, yL = self.sensors_list[s1 - 1].xL, self.sensors_list[s1 - 1].yL
            pa, pb, d = np.array([0, yL]), np.array([xL, yL]), xL
        elif s1 == len(self.sensors_list)+1:
            xR, yR = self.sensors_list[s2 - 1].xR, self.sensors_list[s2 - 1].yR
            pa, pb, d = np.array([self.L, yR]), np.array([xR, yR]), self.L-xR
        elif s2 == len(self.sensors_list)+1:
            xR, yR = self.sensors_list[s1 - 1].xR, self.sensors_list[s1 - 1].yR
            pa, pb, d = np.array([self.L, yR]), np.array([xR, yR]), self.L-xR
        if self.o_adj_matrix[s1, s2] > 0:
            dv = d/self.o_adj_matrix[s1, s2]
            phi = distance.arctan(pb[1]-pa[1], pb[0]-pa[0])
            alpha = self.sensors_list[0].alpha
            lr = self.sensors_list[0].lr
            r = self.sensors_list[0].r

            res = []
            if lr == r and 2*alpha<np.pi:
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi), pa[1]+i*dv*np.sin(phi), phi]))
            elif lr == 2*r*np.sin(alpha) and 2*alpha<np.pi:
                h = np.sqrt(r**2-(lr/2)**2)
                l = np.sqrt(h**2+(dv/2)**2)
                lamda = distance.arctan(2*h, dv)
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi)+l*np.cos(phi+lamda), pa[1]+i*dv*np.sin(phi)+l*np.sin(phi+lamda), angle(phi+3/2*np.pi)]))
            elif lr == 2*r and 2*alpha>np.pi:
                for i in range(int(self.o_adj_matrix[s1, s2])):
                    res.append(np.array([pa[0]+i*dv*np.cos(phi), pa[1]+i*dv*np.sin(phi), angle(phi+np.pi/2)]))
            return res
        else:
            return []

    def mcbf(self, k, dynamic_sens): # k la so barierr, dynamic_sens list cac sensor dong duoc trien khai
        tau = len(dynamic_sens)
        Pk, Nm = self.min_num_mobile_greedy(k)
        print ("Nm %d"%Nm)
        locs = []
        for path in Pk:
            print(path)
            for i in range(len(path)-1):
                locs.extend(self.calculate_target_location(path[i], path[i+1]))
        if tau < len(locs):
            return (None, None), np.inf
        d = np.zeros((tau, Nm))
        for i, dynamic_sen in enumerate(dynamic_sens):
            for j, loc in enumerate(locs):
                d[i,j] = np.sqrt((dynamic_sen.xi-loc[0])**2+(dynamic_sen.yi-loc[1])**2)
        row_ind, col_ind = linear_sum_assignment(d)
        min_cost = d[row_ind, col_ind].sum()
        return locs, (row_ind, col_ind), min_cost
        # locs la vi tri cac target, (row_ind, col_ind) la ghep cap giua dynamic sensor den target, min_cost chi phi minimum

    def add_population(self, num_particles, num_barriers, omega, c1, c2):
        self.population = Population(self, num_particles, num_barriers, omega, c1, c2)




if __name__ == '__main__':
    import  distance
    wbg = WBG(lenght=100, height=20, mode='strong')
    # sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    #s1 = Sensor(3, 3, np.pi / 5, 2, np.pi / 4)
    #s2 = Sensor(8, 3, 5*np.pi / 6, 2, np.pi / 4)
    wbg.create_sensors_randomly(50, r=2, alpha=np.pi)
    wbg.build_WBG()
    wbg.show_matrix()


    # ## debug
    # for path in Pk:
    #     if len(path) > 2:
    #         for i in range(1, len(path)-2):
    #             res = distance.minimum__sectors_distance(wbg.sensors_list[path[i]-1], wbg.sensors_list[path[i+1]-1])
    #             wbg.add_dis_2_points([res[0], res[1]])

    # print("######################################################")
    # Nb, Pk = wbg.max_num_barrier_greedy(8)
    # print(Nb)
    # print(Pk)
    # ## debug
    # for path in Pk:
    #     if len(path) > 2:
    #         for i in range(1, len(path)-2):
    #             res = distance.minimum__sectors_distance(wbg.sensors_list[path[i]-1], wbg.sensors_list[path[i+1]-1])
    #             wbg.add_dis_2_points([res[0], res[1]])
    # wbg.field_show()


    print("######################################################")
    '''
    dynamic_sens = [Sensor(xi=random.uniform(0, wbg.L), yi = random.uniform(0, wbg.H), betai= random.uniform(0, 360), r=3, alpha=60) for _ in range(10)]
    for dynamic_sen in dynamic_sens:
        wbg.add_dynamic(np.array([dynamic_sen.xi, dynamic_sen.yi]))
        ## hien thi cham do tren do thi
        ## moi cham ung voi vi tri sensor dong ban dau

    res = wbg.mcbf(3, dynamic_sens)
    print (res[2])
    for loc in res[0]:
        wbg.add_target(loc[:2])
        ## hien thi cham xanh tren do thi
        ## moi cham ung voi vi tri muc tieu can dat sensor dong
    '''
    Pk, Nm = wbg.min_num_mobile_greedy(3)
    print("######################################################")
    print(Pk)
    print(Nm)
    logging.info('greedy %d'%Nm)
    '''
    for i in range(len(wbg.sensors_list)):
        for j in range(i+1, len(wbg.sensors_list)):
            res = distance.minimum__sectors_distance(wbg.sensors_list[i], wbg.sensors_list[j])
            wbg.add_dis_2_points([res[0], res[1]])'''
    print("######################################################")

    wbg.field_show()
    
    for OMEGA in [0.5,1.0]:
        for C1 in [2.0,0.5,1.0,1.5,0.1]:
            for C2 in [0.1,0.5,1.0,1.5,2.0]:
                print ('CASE: %f, %f, %f'%(OMEGA, C1, C2))
                min_fit = np.inf
                fits = []
                for i in range(30):
                    wbg.add_population(50, 3, OMEGA, C1, C2)  
    #wbg.population.initialize()
    #wbg.population.show()
    #wbg.population.cross_over()
    #wbg.population.clone(5)
    #wbg.population.show()
    #wbg.population.particles[0].fitness(verbose=True)
                    fit = wbg.population.evolve(500)
                    fits.append(fit)
                    if fit < min_fit:
                        min_fit = fit
                    
                    ###START
                    # N = 200
                    # K = 3
                    # def banana(chrome):
                    #     mask = np.ones_like(chrome).astype(float) # mask values, 1 if vertex can be selected, else -np.inf
                    #     mask[0] = -np.inf # s 
                    #     paths = []
                    #     result = 0
                    #     for k in range(K):
                    #         path = [0]
                    #         next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, chrome)])
                    #         if mask[next_v]==-np.inf: # invalid
                    #             path.append(N+1) # direct barrier
                    #             paths.append(path)
                    #             continue
                    #         while next_v <= N:
                    #             mask[next_v] = -np.inf
                    #             path.append(next_v)
                    #             next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, chrome)]) # if negative? => BUG
                    #         mask[next_v] = -np.inf
                    #         path.append(N+1)
                    #         paths.append(path)
                    #     #print (paths)
                    #     result = np.sum([wbg.length(path) for path in paths])
                    #     return result

                    # lb = [1]*(1+N+K)
                    # ub = [100]*(1+N+K)

                    # xopt, fopt = pso(banana, lb, ub)
                    # print ('pyswarm')
                    # print (xopt)
                    # print (fopt)
                    ###END

                logging.info('CASE: %f, %f, %f'%(OMEGA, C1, C2))
                mean_fit = np.mean(fits)
                logging.info('min_fit %d, mean_fit %d'%(min_fit, mean_fit))
                logging.info('#########')



    wbg.field_show()