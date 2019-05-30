import codecs
import matplotlib.patches as patches
import numpy as np
import random
from matplotlib import collections  as mc
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge

from sensor import Sensor


class Sensors_field():
    def __init__(self, lenght, height):
        self.L = lenght
        self.H = height
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

    def create_sensors_randomly(self, num_sensor=100, r=3, alpha=np.pi/3):
        for i in range(0, num_sensor):
            sensor = Sensor(xi=random.uniform(0, self.L), yi=random.uniform(0, self.H), betai=random.uniform(-np.pi, np.pi),
                            r=r, alpha=alpha)
            self.sensors_list.append(sensor)

    def field_show(self, virtual_node=False, filename=None):
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
            if virtual_node:  # for debug:
                for j in range(len(sens.virtual_nodes)):
                    patch = patches.Wedge((sens.xi, sens.yi), sens.r,
                                          (sens.virtual_nodes[j] - sens.alpha) / np.pi * 180,
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
        plt.xlim(xmax=self.L, xmin=0)
        plt.ylim(ymax=self.H, ymin=0)

        x = [p[0] for p in self.dynamics]
        y = [p[1] for p in self.dynamics]
        plt.scatter(x, y, color='red')

        x = [p[0] for p in self.targets]
        y = [p[1] for p in self.targets]
        plt.scatter(x, y, color='green')

        # plt.legend()  # <--- here

        if filename:
            plt.savefig(filename)
        else:
            plt.show()

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
            file.write("S" + "\t" + str(sensor.xi) + "\t" + str(sensor.yi) + "\t" + str(sensor.betai) + "\t" + str(
                sensor.r) + "\t" + str(sensor.alpha) + "\n")
        for sensor in self.mobile_sensors_list:
            file.write("M" + "\t" + str(sensor.xi) + "\t" + str(sensor.yi) + "\t" + str(sensor.betai) + "\t" + str(
                sensor.r) + "\t" + str(sensor.alpha) + "\n")

    def load_sensors_from_txt(self, filename):
        file = codecs.open(filename, mode="r", encoding="utf-8")
        for line in file.readlines():
            try:
                field_array = line.split("\t")
                if field_array[0] == "S":
                    s = Sensor(xi=float(field_array[1]), yi=float(field_array[2]), betai=float(field_array[3]),
                               alpha=float(field_array[5]), r=float(field_array[4]))
                    self.sensors_list.append(s)
            except:
                print("error")

