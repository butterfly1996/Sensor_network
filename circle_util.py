#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from sympy.geometry import Circle, Point

from math import cos, sin, pi, sqrt, atan2, acos
d2r = pi/180

class Geometry(object):
    def circle_intersection(self, circle1, circle2):
        '''
        @summary: calculates intersection points of two circles
        @param circle1: tuple(x,y,radius)
        @param circle2: tuple(x,y,radius)
        @result: tuple of intersection points (which are (x,y) tuple)
        '''
        # return self.circle_intersection_sympy(circle1,circle2)
        x1,y1,r1 = circle1
        x2,y2,r2 = circle2
        # http://stackoverflow.com/a/3349134/798588
        dx,dy = x2-x1,y2-y1
        d = sqrt(dx*dx+dy*dy)
        if d > r1+r2:
            return [] # no solutions, the circles are separate
        if d < abs(r1-r2):
            return [] # no solutions because one circle is contained within the other
        if d == 0 and r1 == r2:
            return [] # circles are coincident and there are an infinite number of solutions

        a = (r1*r1-r2*r2+d*d)/(2*d)
        h = sqrt(r1*r1-a*a)
        xm = x1 + a*dx/d
        ym = y1 + a*dy/d
        xs1 = xm + h*dy/d
        xs2 = xm - h*dy/d
        ys1 = ym - h*dx/d
        ys2 = ym + h*dx/d
        
        if np.isclose(xs1,xs2) and np.isclose(ys1,ys2):
            return [(xs1,ys1)]

        return [(xs1,ys1),(xs2,ys2)]

    def circle_intersection_sympy(self, circle1, circle2):
        
        x1,y1,r1 = circle1
        x2,y2,r2 = circle2
        c1=Circle(Point(x1,y1),r1)
        c2=Circle(Point(x2,y2),r2)
        intersection = c1.intersection(c2)
        return [tuple(map(float, l)) for l in intersection]
    
    def circle_tangency(self, circle, point):
        (Px, Py) = point
        (Cx, Cy, a) = circle

        b = sqrt((Px - Cx)**2 + (Py - Cy)**2)  # hypot() also works here
        if a > b:
            return []
        th = acos(a / b)  # angle theta
        d = atan2(Py - Cy, Px - Cx)  # direction angle of point P from C
        d1 = d + th  # direction angle of point T1 from C
        d2 = d - th  # direction angle of point T2 from C

        T1x = Cx + a * cos(d1)
        T1y = Cy + a * sin(d1)
        T2x = Cx + a * cos(d2)
        T2y = Cy + a * sin(d2)

        if np.isclose(T1x,T2x) and np.isclose(T1y,T2y):
            return [(T1x, T1y)]

        return [(T1x, T1y), (T2x, T2y)]

    def circle_tangency_sympy(self, circle, point):
        point = Point(point)
        (Cx, Cy, a) = circle
        circle = Circle(Point(Cx, Cy), a)
        
        tangency = circle.tangent_lines(point)

        return [tuple(map(float, l.points[1])) for l in tangency]

def test_circle_intersection():
    geom = Geometry()
    print(
        geom.circle_tangency((6,5,3),(8,2)),
        geom.circle_tangency_sympy((6,5,3),(8,2)))
    # np.testing.assert_almost_equal(
    #     geom.circle_intersection((2,0,1),(0,0,1)),
    #     ((1,0),(1,0)))
    # np.testing.assert_almost_equal(
    #     geom.circle_intersection((1,1,1),(3,1,1)),
    #     ((2,1),(2,1)))
    # np.testing.assert_almost_equal(
    #     geom.circle_intersection((0,0,1),(cos(d2r*45)*2,0,1)),
    #     ((cos(d2r*45),-sin(d2r*45)),
    #      (cos(d2r*45),+sin(d2r*45))))
#test_circle_intersection()