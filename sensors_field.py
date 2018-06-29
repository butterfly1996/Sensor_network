from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random
from matplotlib import collections  as mc
# import distance
class Sensor():
    def __init__(self, xi, yi, betai, r, alpha):
        self.xi = xi
        self.yi = yi
        while(betai>np.pi or betai<-np.pi):
            if(betai>np.pi):
                betai-=np.pi
            if (betai < -np.pi):
                betai += np.pi
        while (alpha > np.pi or alpha < -np.pi):
            if (alpha > np.pi):
                alpha -= np.pi
            if (alpha < -np.pi):
                alpha += np.pi
        self.betai = betai
        self.alpha = alpha
        self.r = r
        self.alpha = alpha

        # self.lr = 2 * self.r if alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(alpha))
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
        # for cupple_point in self.pointslist:
        #     lines.append(tuple(cupple_point[0], cupple_point[1]))
        lc = mc.LineCollection(np.array(self.pointslist), linewidths=2)
        ax.add_collection(lc)
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
    def build_WBG(self):
        pass
    def dw(self, vi, vj):# weak distance
        if (vi.xL <= vj.xL and vj.xL <= vi.xR) or (vj.xL <= vi.xL and vi.xL <= vj.xR):
            return 0
        #else
    def ds(self, vi, vj):# strong distance
        pass
    def w(self, vi, vj):# weight
        pass

if __name__ == '__main__':
    import  distance
    sensor_field = Sensors_field(lenght=10, height=10)
    # sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    sensor_field.add_sensor(Sensor(3, 3, np.pi/3, 2, np.pi/4))
    sensor_field.add_sensor(Sensor(8, 3, np.pi/2, 2, np.pi/4))
    res = distance.minimum__sectors_distance(Sensor(3, 3, np.pi / 3, 2, np.pi / 4), Sensor(8, 3, np.pi / 2, 2, np.pi / 4))
    print ("Min: %f"%res[2])
    sensor_field.add_dis_2_points([res[0], res[1]])
    sensor_field.field_show()


