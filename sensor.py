import numpy as np
from math import atan2
import distance
from circle_util import Geometry

from utils import *
geom = Geometry()


class Sensor():
    def __init__(self, xi=None, yi=None, betai=None, r=None, alpha=None, name=None):
        self.name = name
        self.id = None
        self.vid = None
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
        self.lr = 2 * r if self.alpha >= np.pi / 2 else max(r, 2 * r * np.sin(self.alpha))

        # leftmost point
        if -np.pi / 2 <= angle(self.betai - self.alpha) <= np.pi / 2 and -np.pi / 2 <= angle(
                self.betai + self.alpha) <= np.pi / 2:
            self.xL = self.xi
            self.yL = self.yi
        elif angle(self.betai - self.alpha) * angle(self.betai + self.alpha) < 0 and np.abs(
                angle(self.betai - self.alpha) - angle(self.betai + self.alpha)) > np.pi:
            self.xL = self.xi - self.r
            self.yL = self.yi
        else:
            self.xL = min(self.xi + r * np.cos(angle(self.betai - self.alpha)),
                          self.xi + r * np.cos(angle(self.betai + self.alpha)))
            argmin = np.argmin([self.xi + r * np.cos(angle(self.betai - self.alpha)),
                                self.xi + r * np.cos(angle(self.betai + self.alpha))])
            if argmin == 0:
                self.yL = self.yi + r * np.sin(angle(self.betai - self.alpha))
            else:
                self.yL = self.yi + r * np.sin(angle(self.betai + self.alpha))

        # rightmost point
        if (-np.pi <= angle(self.betai - self.alpha) <= -np.pi / 2 or np.pi / 2 <= angle(
                self.betai - self.alpha) <= np.pi) \
                and (-np.pi <= angle(self.betai + self.alpha) <= -np.pi / 2 or np.pi / 2 <= angle(
            self.betai + self.alpha) <= np.pi):
            self.xR = self.xi
            self.yR = self.yi + r * np.sin(angle(self.betai - self.alpha))
        elif angle(self.betai - self.alpha) * angle(self.betai + self.alpha) < 0 and np.abs(
                angle(self.betai - self.alpha) - angle(self.betai + self.alpha)) < np.pi:
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
        # return True if np.isclose(distance.minimum__sectors_distance(self, s2)[2], 0) else False #TODO: re implement
        return distance.is_overlap(self, s2)

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
        d = np.linalg.norm(np.array([self.xi, self.yi]) - np.array([s2.xi, s2.yi]))
        r = self.r
        if np.sqrt(2) * r <= d <= 2 * r:
            intersections = self._get_circle_intersect(s2)

            for intersection in intersections:
                phi = atan2(intersection[1] - self.yi, intersection[0] - self.xi)
                self.add_virtual_node(angle(phi - self.alpha))
                self.add_virtual_node(angle(phi + self.alpha))
        elif r <= d <= np.sqrt(2) * r:
            intersections = self._get_circle_intersect(s2)
            tangencies = s2._get_circle_tangency((self.xi, self.yi))

            for intersection in intersections:
                phi = atan2(intersection[1] - self.yi, intersection[0] - self.xi)
                self.add_virtual_node(angle(phi - self.alpha))
                self.add_virtual_node(angle(phi + self.alpha))

            phi0 = atan2(tangencies[0][1] - self.yi, tangencies[0][0] - self.xi)
            phi1 = atan2(tangencies[1][1] - self.yi, tangencies[1][0] - self.xi)
            s0s2 = atan2(s2.yi - self.yi, s2.xi - self.xi)
            if is_between_angles(phi0 - self.alpha, s0s2, phi0, phi1):
                self.add_virtual_node(angle(phi0 + self.alpha))
            else:
                self.add_virtual_node(angle(phi0 - self.alpha))
            if is_between_angles(phi1 - self.alpha, s0s2, phi0, phi1):
                self.add_virtual_node(angle(phi1 + self.alpha))
            else:
                self.add_virtual_node(angle(phi1 - self.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[0][1] - s2.yi, tangencies[0][0] - s2.xi) - s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[0][1] - s2.yi, tangencies[0][0] - s2.xi) + s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[1][1] - s2.yi, tangencies[1][0] - s2.xi) - s2.alpha))
            s2.add_virtual_node(angle(atan2(tangencies[1][1] - s2.yi, tangencies[1][0] - s2.xi) + s2.alpha))
        elif d <= r:
            intersections = self._get_circle_intersect(s2)

            for intersection in intersections:
                phi = atan2(intersection[1] - self.yi, intersection[0] - self.xi)
                self.add_virtual_node(angle(phi - self.alpha))
                self.add_virtual_node(angle(phi + self.alpha))

            phi = atan2(s2.yi - self.yi, s2.xi - self.xi)
            self.add_virtual_node(angle(phi - self.alpha))
            self.add_virtual_node(angle(phi + self.alpha))
        elif d > 2 * r:
            phi = atan2(s2.yi - self.yi, s2.xi - self.xi)
            self.add_virtual_node(phi)

    def update_view(self, new_beta):
        self.__init__(self.xi, self.yi, new_beta, self.r, self.alpha)