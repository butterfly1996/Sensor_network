from scipy.optimize import linear_sum_assignment
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random
from matplotlib import collections  as mc
import matplotlib.patches as patches
from pso_ga import *
import codecs
import logging
from circle_util import Geometry
from math import atan2
import itertools
import multiprocessing
import distance
#from pyswarm import pso

logging.basicConfig(filename='tunning.log',level=logging.INFO)
geom = Geometry()

np.random.seed(23)

# import distance
def angle(value):
    # chuan hoa gia tri goc tu -pi den pi
    return (value+np.pi) % (2*np.pi) - np.pi

def is_between_angles(value, bisector, angle1, angle2):
    angle_min = min(angle1, angle2)
    angle_max = max(angle1, angle2)
    if angle_max-angle_min < np.pi:
        if np.isclose(bisector, (angle_max+angle_min)/2):
            return angle_min <= value <= angle_max
        else:
            return not angle_min <= value <= angle_max
    elif value < 0:
        if bisector < 0:
            if np.isclose(bisector, (angle_max-2*np.pi+angle_min)/2):
                return angle_max-2*np.pi <= value <= angle_min
            else:
                return not angle_max-2*np.pi <= value <= angle_min
        else:
            if np.isclose(bisector-2*np.pi, (angle_max-2*np.pi+angle_min)/2):
                return angle_max-2*np.pi <= value <= angle_min
            else:
                return not angle_max-2*np.pi <= value <= angle_min
    else:
        if bisector < 0:
            if np.isclose(bisector, (angle_max-2*np.pi+angle_min)/2):
                return angle_max-2*np.pi <= value-2*np.pi <= angle_min
            else:
                return not angle_max-2*np.pi <= value-2*np.pi <= angle_min
        else:
            if np.isclose(bisector-2*np.pi, (angle_max-2*np.pi+angle_min)/2):
                return angle_max-2*np.pi <= value-2*np.pi <= angle_min
            else:
                return not angle_max-2*np.pi <= value-2*np.pi <= angle_min

class Sensor():
    def __init__(self, xi=None, yi=None, betai=None, r=None, alpha=None, name=None):
        self.name = name
        if self.name in ['s', 't']:
            self.xi = xi
            self.xL = self.xR = self.xi
            self.virtual_nodes = [self]
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
        
        self.virtual_nodes = []

    def overlap(self, s2):
        #return True if np.isclose(distance.minimum__sectors_distance(self, s2)[2], 0) else False #TODO: re implement
        return distance.is_overlap(self,s2)
    def add_virtual_node(self, beta):
        self.virtual_nodes.append(Sensor(self.xi, self.yi, beta, self.r, self.alpha))
    def _get_circle_intersect(self, s2):
        intersections = geom.circle_intersection((self.xi, self.yi, self.r), (s2.xi, s2.yi, s2.r))
        return intersections
    def _get_circle_tangency(self, p):
        tangency = geom.circle_tangency((self.xi, self.yi, self.r), p)
        return tangency
    def _get_boundary_intersect(self, x):
        pass
    def set_virtual_nodes(self, s2):
        d = np.linalg.norm(np.array([self.xi, self.yi])-np.array([s2.xi, s2.yi]))
        r = self.r
        if np.sqrt(2)*r <= d <= 2*r:
            intersections = self._get_circle_intersect(s2)
            
            for intersection in intersections:
                phi = atan2(intersection[1]-self.yi, intersection[0]-self.xi)
                self.add_virtual_node(angle(phi-self.alpha))
                self.add_virtual_node(angle(phi+self.alpha))
        elif r <= d <= np.sqrt(2)*r:
            intersections = self._get_circle_intersect(s2)
            tangencies = s2._get_circle_tangency((self.xi,self.yi))

            for intersection in intersections:
                phi = atan2(intersection[1]-self.yi, intersection[0]-self.xi)
                self.add_virtual_node(angle(phi-self.alpha))
                self.add_virtual_node(angle(phi+self.alpha))

            phi0 = atan2(tangencies[0][1]-self.yi, tangencies[0][0]-self.xi)
            phi1 = atan2(tangencies[1][1]-self.yi, tangencies[1][0]-self.xi)
            s0s2 = atan2(s2.yi-self.yi, s2.xi-self.xi)
            if is_between_angles(phi0-self.alpha, s0s2, phi0, phi1):
                self.add_virtual_node(angle(phi0+self.alpha))
            else:
                self.add_virtual_node(angle(phi0-self.alpha))
            if is_between_angles(phi1-self.alpha, s0s2, phi0, phi1):
                self.add_virtual_node(angle(phi1+self.alpha))
            else:
                self.add_virtual_node(angle(phi1-self.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[0][1]-s2.yi, tangencies[0][0]-s2.xi)-s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[0][1]-s2.yi, tangencies[0][0]-s2.xi)+s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[1][1]-s2.yi, tangencies[1][0]-s2.xi)-s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[1][1]-s2.yi, tangencies[1][0]-s2.xi)+s2.alpha))
        elif d <= r:
            intersections = self._get_circle_intersect(s2)

            for intersection in intersections:
                phi = atan2(intersection[1]-self.yi, intersection[0]-self.xi)
                self.add_virtual_node(angle(phi-self.alpha))
                self.add_virtual_node(angle(phi+self.alpha))
            
            phi = atan2(s2.yi-self.yi, s2.xi-self.xi)
            self.add_virtual_node(angle(phi-self.alpha))
            self.add_virtual_node(angle(phi+self.alpha))

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
    def field_show(self, virtual_node=False):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.subplot()
        plt.gca().set_aspect('equal', adjustable='box')

        for i in range(0, len(self.sensors_list)):
            sens = self.sensors_list[i]
            color = np.array([random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)])
            patch = patches.Wedge((sens.xi, sens.yi), sens.r, (sens.betai - sens.alpha) / np.pi * 180,
                          (sens.betai + sens.alpha) / np.pi * 180,
                          color=color, alpha=0.6,
                          label=str(i + 1))
            ax.add_patch(patch)
            if virtual_node: # for debug:
                for j in range(len(sens.virtual_nodes)):
                    patch = patches.Wedge((sens.xi, sens.yi), sens.r, (sens.virtual_nodes[j] - sens.alpha) / np.pi * 180,
                            (sens.virtual_nodes[j] + sens.alpha) / np.pi * 180,
                            color=color, alpha=0.3,
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


class DBG(WBG):
    def __init__(self, lenght, height):
        WBG.__init__(self, lenght, height, mode='strong')
        self.H = height
        self.L = lenght
    def build_virtual_nodes(self):
        for i in range(len(self.sensors_list)):
            for j in range(i+1, len(self.sensors_list)):
                u = self.sensors_list[i]
                v = self.sensors_list[j]
                if len(geom.circle_intersection((u.xi,u.yi,u.r),(v.xi,v.yi,v.r)))>0 and i!=j:
                    u.set_virtual_nodes(v)
                    v.set_virtual_nodes(u)
        # for s and t
        for i in range(len(self.sensors_list)):
            u = self.sensors_list[i]
            if u.xi <= u.r:
                u.add_virtual_node(angle(+np.arccos(-u.xi/u.r)-u.alpha))
                u.add_virtual_node(angle(-np.arccos(-u.xi/u.r)+u.alpha))
                u.add_virtual_node(angle(+np.arccos(-u.xi/u.r)+u.alpha))
                u.add_virtual_node(angle(-np.arccos(-u.xi/u.r)-u.alpha))
        for i in range(len(self.sensors_list)):
            u = self.sensors_list[i]
            if u.xi >= self.L-u.r:
                u.add_virtual_node(angle(+np.arccos((self.L-u.xi)/u.r)-u.alpha))
                u.add_virtual_node(angle(-np.arccos((self.L-u.xi)/u.r)+u.alpha))
                u.add_virtual_node(angle(+np.arccos((self.L-u.xi)/u.r)+u.alpha))
                u.add_virtual_node(angle(-np.arccos((self.L-u.xi)/u.r)-u.alpha))
    def rotated_angle(self, actual_node, virtual_node):
        actual = actual_node.betai
        virtual = virtual_node.betai
        return min(abs(angle(actual-virtual)),2*np.pi-abs(angle(actual-virtual)))
    def fill_mat(self, params):
        u, v, ai, aj, ii, jj = params
        if u.virtual_nodes[ii].overlap(v.virtual_nodes[jj]):
            self.adj_matrix[ai+ii,aj+jj] = self.rotated_angle(v,v.virtual_nodes[jj])
    def build_dbg(self):
        num_virtuals = np.sum([len(s.virtual_nodes) for s in self.sensors_list]) + 2
        self.adj_matrix = np.full((num_virtuals, num_virtuals), np.inf)
        print ('nvir %d'%num_virtuals)
        
        #print ('s->')
        # s->
        aj = 1
        for j, v in enumerate(self.sensors_list):
            if v.xi <= v.r:
                for jj, vv in enumerate(v.virtual_nodes):
                    if np.equal(super().ds(self.s, vv), 0):
                        self.adj_matrix[0][aj+jj] = self.rotated_angle(v, vv)
            aj += len(v.virtual_nodes)
        
        #print ('t->')
        # -> t
        ai = 1
        for i, u in enumerate(self.sensors_list):
            if u.xi >= self.L-u.r:
                for ii, uu in enumerate(u.virtual_nodes):
                    if np.equal(super().ds(uu, self.t), 0):
                        self.adj_matrix[ai+ii][num_virtuals-1] = self.rotated_angle(u, uu)
            ai += len(u.virtual_nodes)
        
        pool = multiprocessing.Pool(12)
        # actual sensors
        #print ('aaaa')
        ai = aj = 1
        for i, u in enumerate(self.sensors_list):
            aj = 1
            for j, v in enumerate(self.sensors_list):
                print ('i j ', i, j)
                if i != j and u.xi < v.xi and len(geom.circle_intersection((u.xi,u.yi,u.r),(v.xi,v.yi,v.r)))>0:
                    '''
                    for ii, uu in enumerate(u.virtual_nodes):
                        for jj, vv in enumerate(v.virtual_nodes):
                            if uu.overlap(vv):
                                self.adj_matrix[ai+ii][aj+jj] = rotated_angle(v, vv)
                    '''
                    ii_list = range(len(u.virtual_nodes))
                    jj_list = range(len(v.virtual_nodes))
                    paramlist = list(itertools.product([u],[v],[ai],[aj],ii_list,jj_list))
                    pool.map(self.fill_mat,paramlist)
                aj += len(v.virtual_nodes)
            ai += len(u.virtual_nodes)
    def show_matrix(self):
        print(self.adj_matrix)
    def dfs(self, s, visited, t):
        if s not in visited:
            visited.append(s)
            if s == t:
                return True
            for u in range(len(self.adj_matrix[s])):
                if self.adj_matrix[s][u] != np.inf:
                    self.dfs(u, visited, t)
            del visited[-1]
        return False


if __name__ == '__main__':
    for n in np.arange(40,90,10):
        rate = 0.0
        for simulation in range(20):
            print ('sim %d'%simulation)
            dbg = DBG(lenght=500, height=100)
            dbg.create_sensors_randomly(num_sensor=n, r=50, alpha=np.pi/6)
            #dbg.add_sensor(Sensor(4,5,0,2,np.pi/6))
            #dbg.add_sensor(Sensor(5,3,-np.pi/2,2,np.pi/6))
            dbg.field_show()
            #dbg.build_WBG()
            print ('build vir')
            dbg.build_virtual_nodes()
            print ('build dbg')
            dbg.build_dbg()
            #dbg.show_matrix()            
            rate += dbg.dfs(0, [], np.sum([len(s.virtual_nodes) for s in dbg.sensors_list]) + 1)
            print (rate)
        rate/=20.0
        print ("rate after %d: %f"%(n, rate))