from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random

class Sensor():
    def __init__(self, xi, yi, betai, r, alpha):
        self.xi = xi
        self.yi = yi
        self.betai = betai
        self.alpha = alpha
        self.r = r
        self.lr = 2 * self.r if self.alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(self.alpha))
        self.xR = max(self.xi, self.xi+self.r*np.cos(self.betai-self.alpha), self.xi+self.r*np.cos(self.betai+self.alpha),
            self.xi+self.r if -self.alpha <= -self.betai and -self.betai <= self.alpha else 0)
        self.xL = min(self.xi, self.xi+self.r*np.cos(self.betai-self.alpha), self.xi+self.r*np.cos(self.betai+self.alpha),
            self.xi+self.r if -self.alpha <= -self.betai and -self.betai <= self.alpha else self.xi+9*self.r)
class Sensors_field():
    def __init__(self, lenght, height, num_station, num_mobile):
        self.L=lenght
        self.H=height
        self.n=num_station
        self.to = num_mobile
        # alpha la nua goc nhin cua sensor
        # self.alpha = alpha

        # self.r = r

        #sensor la danh sach cac sensor gom cac truong betai la huong voi tia Ox, li la location (xi, yi)
        self.sensors_list = []
        #tap cac duong di dinh khong giao nhau co tong do dai nho nhat
        self.P_sao_q=[]
        #k-auxiliary set cua P_sao_q
        self.P_k_q = []
        #tap k barriers toi uu cua bai toan Min-num-mobile(k)
        self.P_k = []
    def create_sensors_randomly(self, num_sensor = 100, r=3, alpha=60):
        for i in range(0, num_sensor):
            sensor = Sensor(xi=random.uniform(0, self.L), yi = random.uniform(0, self.H), betai= random.uniform(0, 360), r=r, alpha=alpha)
            self.sensors_list.append(sensor)
    def field_show(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for sens in self.sensors_list:
            fov = Wedge((sens.xi, sens.yi), sens.r, sens.betai-sens.alpha, sens.betai+sens.alpha, color="r", alpha=0.5)
            ax.add_artist(fov)
        plt.xlim(xmax = self.L, xmin=0)
        plt.ylim(ymax = self.H, ymin = 0)
        plt.show()
        pass
# class WBG(Sensors_field):
#     def build_WBG():
#         pass
#     def dw(vi, vj):# weak distance
#         if (vi.xL <= vj.xL and vj.xL <= vi.xR) or (vj.xL <= vi.xL and vi.xL <= vj.xR):
#             return 0
#         #else
#     def ds(vi, vj):# strong distance
#         pass
#     def w(vi, vj):# weight
#         pass

if __name__ == '__main__':
    sensor_field = Sensors_field(lenght=30, height=10, num_station=12, num_mobile=20)
    sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    sensor_field.field_show()


