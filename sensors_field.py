from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random
from matplotlib import collections  as mc

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

        print(self.xL, self.xR)

        # self.lr = 2 * self.r if alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(alpha))
    def overlap(self, s2):
        return True if distance.minimum__sectors_distance(self, s2) == 0 else False
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
        for sens in self.sensors_list:
            fov = Wedge((sens.xi, sens.yi), sens.r, (sens.betai-sens.alpha)/np.pi*180, (sens.betai+sens.alpha)/np.pi*180, color="r", alpha=0.5)
            ax.add_artist(fov)
        lines = []
        # show = True
        # for cupple_point in self.pointslist:
        #      if cupple_point == None:
        #          show = False
        #          break
        # if show == True:
        try:
            lc = mc.LineCollection(np.array(self.pointslist), linewidths=2)
            ax.add_collection(lc)
        except ValueError:
            pass
        plt.xlim(xmax = self.L, xmin=0)
        plt.ylim(ymax = self.H, ymin = 0)
        plt.show()
        pass
    def add_sensor(self, sensor):
        self.sensors_list.append(sensor)
    def add_dis_2_points(self, points):
        self.pointslist.append(points)
        print(self.pointslist)
        
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
            return distance.minimum__sectors_distance(vi, vj)
    def w(self, vi, vj):# weight
        try:
            lr = vi.lr
        except:
            lr = vj.lr
        if self.mode == 'weak':
            return np.ceil(self.dw(vi, vj)/lr)
        else:
            return np.ceil(self.ds(vi, vj)/lr)
    def show_matrix(self):
        print(self.adj_matrix)

if __name__ == '__main__':
    import  distance
    wbg = WBG(lenght=10, height=10)
    # sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    #s1 = Sensor(3, 3, np.pi / 5, 2, np.pi / 4)
    #s2 = Sensor(8, 3, 5*np.pi / 6, 2, np.pi / 4)
    wbg.create_sensors_randomly(10)
    wbg.build_WBG()
    wbg.show_matrix()
    wbg.field_show()