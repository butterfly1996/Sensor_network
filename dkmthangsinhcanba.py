from matplotlib import pyplot as plt
import numpy as np
from matplotlib.patches import Wedge
import random

class Sensor():
    def __init__(self, xi, yi, betai, r, alpha):
        self.xi = xi
        self.yi = yi
        self.betai = betai
        self.r = r
        self.alpha = alpha
        self.lr = 2 * self.r if alpha >= np.pi / 2 else np.max(self.r, 2 * r * np.sin(alpha))
class Sensors_field():
    def __init__(self, lenght, height):
        self.L=lenght
        self.H=height
        # self.n=num_station
        # self.to = num_mobile
        # # alpha la nua goc nhin cua sensor
        # self.alpha = alpha
        #
        # self.r = r

        #sensor la danh sach cac sensor gom cac truong betai la huong voi tia Ox, li la location (xi, yi)
        self.sensors_list = []
        # mobile sensor
        self.mobile_sensor = []
        #tap cac duong di dinh khong giao nhau co tong do dai nho nhat
        self.P_sao_q=[]
        #k-auxiliary set cua P_sao_q
        self.P_k_q = []
        #tap k barriers toi uu cua bai toan Min-num-mobile(k)
        self.P_k = []
    def create_sensors_randomly(self, num_sensor = 100, r=3):
        for i in range(0, num_sensor):
            sensor = Sensor(xi=random.uniform(0, self.L), yi = random.uniform(0, self.H), betai= random.uniform(0, 360), r=r)
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
    def add_sensor(self, sensor):
        self.sensors_list.append(sensor)
if __name__ == '__main__':
    # sensor_field = Sensors_field(lenght=30, height=10, num_station=12, num_mobile=20, alpha=60)
    # sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3)
    sensor_field = Sensors_field(10, 10)
    alpha = 70
    sensor_field.add_sensor(sensor=Sensor(3, 3, 30, 2, alpha))
    sensor_field.add_sensor(sensor=Sensor(8, 3, 30, 2, alpha))
    sensor_field.field_show()


