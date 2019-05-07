import numpy as np


def arctan(y, x):
    # if x >= 0:
    #     return np.arctan(y/x)
    if x == 0:
        if y > 0:
            return np.pi / 2
        else:
            return np.pi * 3 / 2
    if x > 0:
        return np.arctan(y/x)
    else:
        return np.pi+np.arctan(y/x)


def angle(value):
    # chuan hoa gia tri goc tu -pi den pi
    return (value) % (2*np.pi)
    # return (value+np.pi) % (2*np.pi) - np.pi


def is_between_angles(value, angle1, angle2):
    angle_min = min(angle1, angle2)
    angle_max = max(angle1, angle2)
    # angle_min = min(angle(angle1), angle(angle2))
    # angle_max = max(angle(angle1), angle(angle2))
    # return angle_min <= angle(value) <= angle_max
    if angle_max-angle_min < np.pi:
        return angle_min <= value <= angle_max
    elif value < 0:
        return angle_max-2*np.pi <= value <= angle_min
    else:
        return angle_max-2*np.pi <= value-2*np.pi <= angle_min

def is_between_anglesHai(value, beta, alpha):
    return abs(angle(beta - value)) <= alpha

SMALL_NUM = 0.00000001
def closestDistanceBetweenLinesHai(a0,a1,b0,b1):
    u = a1 - a0
    v = b1 - b0
    w = a0 - b0
    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a * c - b * b
    sc, sN, sD = D, D, D
    tc, tN, tD = D, D, D

    # compute the line parameters of the two closest points
    if D < SMALL_NUM:# the lines are almost parallel
        sN = 0.0# force using point P0 on segment S1
        sD = 1.0# to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:# get the closest points on the infinite lines
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if (sN < 0.0): # sc < 0 = > the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD: # sc > 1  = > the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0: # tc < 0 = > the t=0 edge is visible
        tN = 0.0
        #recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD: # tc > 1  = > the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if ((-d + b) < 0.0):
            sN = 0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d +  b)
            sD = a
    # finally do the division to get sc and tc
    sc = 0 if abs(sN) < SMALL_NUM else sN / sD
    tc = 0 if abs(tN) < SMALL_NUM else tN / tD

    # get the difference of the two closest points Vector
    dP = w + (sc * u) - (tc * v) # =  S1(sc) - S2(tc)

    return np.linalg.norm(dP)# return the closest distance


def closestDistanceBetweenLines(a0,a1,b0,b1,clampAll=True,clampA0=True,clampA1=True,clampB0=True,clampB1=True):

    ''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
        Return the closest points on each segment and their distance
    '''
    a0 = np.array([a0[0], a0[1], 0])
    b0 = np.array([b0[0], b0[1], 0])
    a1 = np.array([a1[0], a1[1], 0])
    b1 = np.array([b1[0], b1[1], 0])

    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True


    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)

    _A = A / magA
    _B = B / magB

    cross = np.cross(_A, _B);
    denom = np.linalg.norm(cross)**2


    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A,(b0-a0))

        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))

            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return [a0[:2], b0[:2], np.linalg.norm(a0-b0)]
                    return [a0[:2], b1[:2], np.linalg.norm(a0-b1)]


            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return [a1[:2], b0[:2], np.linalg.norm(a1-b0)]
                    return [a1[:2], b1[:2], np.linalg.norm(a1-b1)]


        # Segments overlap, return distance between parallel segments
        return [None, None, np.linalg.norm(((d0*_A)+a0)-b0)]



    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0);
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])

    t0 = detA/denom;
    t1 = detB/denom;

    pA = a0 + (_A * t0) # Projected closest point on segment A
    pB = b0 + (_B * t1) # Projected closest point on segment B


    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1

        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1

        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B,(pA-b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)

        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A,(pB-a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)


    return [pA[:2], pB[:2], np.linalg.norm(pA-pB)]


def closesDistanceBetweenArcs(center1, center2, beta1, beta2, r, alpha):
    res = []
    # arc1
    endpoint1 = center1 + np.array([r * np.cos(beta1 - alpha), r * np.sin(beta1 - alpha)])
    endpoint2 = center1 + np.array([r * np.cos(beta1 + alpha), r * np.sin(beta1 + alpha)])
    # arc2
    endpoint3 = center2 + np.array([r * np.cos(beta2 - alpha), r * np.sin(beta2 - alpha)])
    endpoint4 = center2 + np.array([r * np.cos(beta2 + alpha), r * np.sin(beta2 + alpha)])
    #1: endpoints of two arcs
    res.append([endpoint1, endpoint3, np.linalg.norm(endpoint1-endpoint3)])
    res.append([endpoint1, endpoint4, np.linalg.norm(endpoint1-endpoint4)])
    res.append([endpoint2, endpoint3, np.linalg.norm(endpoint2-endpoint3)])
    res.append([endpoint2, endpoint4, np.linalg.norm(endpoint2-endpoint4)])
    #2: endpoint to interior
    phi = arctan((endpoint3[1]-center1[1]), (endpoint3[0]-center1[0]))
    if (endpoint3[0]-center1[0]) == 0:
        phi = np.sign(endpoint3[1]-center1[1])*np.pi/2
    if is_between_angles(phi, angle(beta1-alpha), angle(beta1+alpha)):
        d = np.linalg.norm(endpoint3-center1)-r
        if d > 0:
            res.append([endpoint3, np.array([center1[0]+r*np.cos(phi), center1[1]+r*np.sin(phi)]), d])
        else:
            res.append([None, None, 0])

    phi = arctan((endpoint4[1] - center1[1]) , (endpoint4[0] - center1[0]))
    if (endpoint4[0]-center1[0]) == 0:
        phi = np.sign(endpoint4[1]-center1[1])*np.pi/2
    if is_between_angles(phi, angle(beta1-alpha), angle(beta1+alpha)):
        d = np.linalg.norm(endpoint4 - center1) - r
        if d > 0:
            res.append([endpoint4, np.array([center1[0] + r * np.cos(phi), center1[1] + r * np.sin(phi)]), d])
        else:
            res.append([None, None, 0])

    phi = arctan((endpoint1[1] - center2[1]) , (endpoint1[0] - center2[0]))
    if (endpoint1[0] - center2[0]) == 0:
        phi = np.sign(endpoint1[1] - center2[1])*np.pi/2
    if is_between_angles(phi, angle(beta2 - alpha), angle(beta2 + alpha)):
        d = np.linalg.norm(endpoint1 - center2) - r
        if d > 0:
            res.append([endpoint1, np.array([center2[0] + r * np.cos(phi), center2[1] + r * np.sin(phi)]), d])
        else:
            res.append([None, None, 0])
    phi = arctan((endpoint2[1] - center2[1]) , (endpoint2[0] - center2[0]))
    if (endpoint2[0] - center2[0]) == 0:
        phi = np.sign(endpoint2[1] - center2[1])*np.pi/2
    if is_between_angles(phi, angle(beta2 - alpha), angle(beta2 + alpha)):
        d = np.linalg.norm(endpoint2 - center2) - r
        if d > 0:
            res.append([endpoint2, np.array([center2[0] + r * np.cos(phi), center2[1] + r * np.sin(phi)]), d])
        else:
            res.append([None, None, 0])
    #3: interiors
    phi = arctan((center2[1]-center1[1]),(center2[0]-center1[0]))
    if center2[0]-center1[0] == 0:
        phi = np.sign(center2[1]-center1[1])*np.pi/2
    phi3 = arctan((endpoint3[1] - center1[1]) , (endpoint3[0] - center1[0]))
    if endpoint3[0] - center1[0] == 0:
        phi3 = np.sign(endpoint3[1] - center1[1])*np.pi/2
    phi4 = arctan((endpoint4[1] - center1[1]) , (endpoint4[0] - center1[0]))
    if endpoint4[0] - center1[0] == 0:
        phi4 = np.sign(endpoint4[1] - center1[1])*np.pi/2
    phi1 = arctan((endpoint1[1] - center2[1]) , (endpoint1[0] - center2[0]))
    if endpoint1[0] - center2[0] == 0:
        phi1 = np.sign(endpoint1[1] - center2[1])*np.pi/2
    phi2 = arctan((endpoint2[1] - center2[1]) , (endpoint2[0] - center2[0]))
    if endpoint2[0] - center2[0] == 0:
        phi2 = np.sign(endpoint2[1] - center2[1])*np.pi/2

    def sign(value):
        if value >= 0:
            return 1
        else:
            return -1
    if is_between_angles(phi, angle(beta1-alpha), angle(beta1+alpha)) and is_between_angles(angle(phi-sign(phi)*np.pi), angle(beta2-alpha), angle(beta2+alpha))\
        and (is_between_angles(phi3, angle(beta1-alpha), angle(beta1+alpha)) or is_between_angles(phi4, angle(beta1-alpha), angle(beta1+alpha)))\
        and (is_between_angles(phi1, angle(beta2-alpha), angle(beta2+alpha)) or is_between_angles(phi2, angle(beta2-alpha), angle(beta2+alpha))):
        d = np.linalg.norm(center1-center2)-2*r
        if d > 0:
            res.append([center1+r*np.array([np.cos(phi), np.sin(phi)]), center2+r*np.array([np.cos(phi-np.pi), np.sin(phi-np.pi)]), d])
        else:
            res.append([None, None, 0])
    #4: intersection
    dx, dy = center2[0] - center1[0], center2[1] - center1[1]
    d = np.sqrt(dx * dx + dy * dy)
    if d < 2*r:
        a = d / 2
        h = np.sqrt(r**2 - a**2)
        xm = center1[0] + a * dx / d
        ym = center1[1] + a * dy / d
        xs1 = xm + h * dy / d
        xs2 = xm - h * dy / d
        ys1 = ym - h * dx / d
        ys2 = ym + h * dx / d
        if is_between_angles(arctan((ys1-center1[1]),(xs1-center1[0])), angle(beta1-alpha), angle(beta1+alpha)) \
                and is_between_angles(arctan((ys1-center2[1]),(xs1-center2[0])), angle(beta2-alpha), angle(beta2+alpha)):
            res.append([np.array([xs1, ys1]), np.array([xs1, ys1]), 0])
        if is_between_angles(arctan((ys2-center1[1]),(xs2-center1[0])), angle(beta1-alpha), angle(beta1+alpha)) \
                and is_between_angles(arctan((ys2-center2[1]),(xs2-center2[0])), angle(beta2-alpha), angle(beta2+alpha)):
            res.append([np.array([xs2, ys2]), np.array([xs2, ys2]), 0])
    res = np.array(res)
    return res[np.argmin(res[:, 2])]


def closestDistanceBetweenArcAndLine(a0, a1, center, beta, r, alpha):
    ### a0 a1 la 2 dau mut cua doan thang
    ### center beta r alpha la thong so cua day cung
    res = []
    # ends
    endpoint1 = center + np.array([r * np.cos(beta - alpha), r * np.sin(beta - alpha)])
    endpoint2 = center + np.array([r * np.cos(beta + alpha), r * np.sin(beta + alpha)])
    res.append([a0, endpoint1, np.linalg.norm(a0-endpoint1)])
    res.append([a0, endpoint2, np.linalg.norm(a0-endpoint2)])
    res.append([a1, endpoint1, np.linalg.norm(a1-endpoint1)])
    res.append([a1, endpoint2, np.linalg.norm(a1-endpoint2)])
    # interiors
    #x = (a0*(a1-b1)**2+c0*(b0-a0)**2+(a1-b1)*(b0-a0)*(c1-a1)) / ((a1-b1)**2+(a0-b0)**2)
    #y = c1-(c0-x)(b0-a0)/(a1-b1)

    x = (a0[0] * (a0[1] - a1[1]) ** 2 + center[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (center[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = center[1] - (center[0] - x)*(a1[0] - a0[0]) / (a0[1] - a1[1])
    d = np.linalg.norm(np.array([x, y]) - center) - r
    phi = arctan((y - center[1]) , (x - center[0]))
    if x - center[0] == 0:
        phi = np.sign(y - center[1])*np.pi/2
    if (a0[0]<= x <=a1[0] or a1[0] <= x <= a0[0]) and (a0[1]<= y <=a1[1] or a1[1] <= y <= a0[1]) and is_between_angles(phi, angle(beta -alpha), angle(beta + alpha)):
        if d > 0:
            res.append([center+r*np.array([np.cos(phi), np.sin(phi)]), np.array([x,y]), d])
        else:
            res.append([None, None, 0]) # intersection
    # end of arc to interior of line
    x = (a0[0] * (a0[1] - a1[1]) ** 2 + endpoint1[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (
    endpoint1[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = endpoint1[1] - (endpoint1[0] - x)*(a1[0] - a0[0]) / (a0[1] - a1[1])
    if (a0[0]<= x <=a1[0] or a1[0] <= x <= a0[0]) and (a0[1]<= y <=a1[1] or a1[1] <= y <= a0[1]):
        res.append([endpoint1, np.array([x, y]), np.linalg.norm(np.array([x, y]) - endpoint1)])

    x = (a0[0] * (a0[1] - a1[1]) ** 2 + endpoint2[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (
        endpoint2[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = endpoint2[1] - (endpoint2[0] - x)*(a1[0] - a0[0]) / (a0[1] - a1[1])
    if (a0[0] <= x <= a1[0] or a1[0] <= x <= a0[0]) and (a0[1] <= y <= a1[1] or a1[1] <= y <= a0[1]):
        res.append([endpoint2, np.array([x, y]), np.linalg.norm(np.array([x, y]) - endpoint2)])
    # end of line to interior of arc
    phi0 = arctan((a0[1] - center[1]) , (a0[0] - center[0]))
    if a0[0] - center[0] == 0:
        phi0 = np.sign(a0[1] - center[1])*np.pi/2
    phi1 = arctan((a1[1] - center[1]) , (a1[0] - center[0]))
    if a1[0] - center[0] == 0:
        phi1 = np.sign(a1[1] - center[1])*np.pi/2
    d = np.linalg.norm(a0 - center) - r
    if is_between_angles(phi0, angle(beta -alpha), angle(beta + alpha)):
        if d > 0:
            res.append([a0, center+r*np.array([np.cos(phi0), np.sin(phi0)]), d])
        else:
            res.append([None, None, 0])
    d = np.linalg.norm(a1 - center) - r
    if is_between_angles(phi1, angle(beta -alpha), angle(beta + alpha)):
        if d > 0:
            res.append([a1, center+r*np.array([np.cos(phi1), np.sin(phi1)]), d])
        else:
            res.append([None, None, 0])
    # intersection
    # A = 1+(a1[1]-a0[1])**2/(a1[0]-a0[0])**2
    # B = -2*center[0]-(2*a0[0]*(a1[1]-a0[1])**2+2*(a1[1]-a0[1])*(a0[1]-center[1])*(a1[0]-a0[0]))/(a1[0]-a0[0])**2
    # C = center[0]*2+(a0[0]**2*(a1[1]-a0[1])**2+(a0[1]-center[1])**2*(a1[0]-a0[0])**2-2*a0[0]*(a1[1]-a0[1])*(a0[1]-center[1])*(a1[0]-a0[0]))/(a1[0]-a0[0])**2
    # delta = B**2-4*A*C
    # if delta >= 0:
    #     res.append(0)
    res = np.array(res)
    return res[np.argmin(res[:, 2])]

SMALL_ANGLE = 0.0001
def is_point_in_fan(point, center, beta, r, alpha):
    if np.linalg.norm(point - center) > r:
        return False
    phi = np.arctan((point[1]-center[1])/(point[0]-center[0]))
    if point[1]-center[1] < 0:
        phi += np.pi
    if alpha >= np.pi - SMALL_ANGLE: # de phong sai so truong hop gan voi pi (Tuc ca duong tron)
        return True
    arg_min = min((beta - alpha)%(2*np.pi), (beta + alpha)%(2*np.pi))
    arg_max = max((beta - alpha)%(2*np.pi), (beta + alpha)%(2*np.pi))
    return arg_min <= phi%(2*np.pi) <= arg_max

def closest_distance_between_point2arc(point, center, beta, r, alpha, code1=-1, code2=-1):
    if is_point_in_fan(point, center, beta, r, alpha):
        return 0
        # return [None, None, 0]
    end_point1 = np.array([center[0] + r*np.cos(beta - alpha), center[1] + r*np.sin(beta - alpha)])
    end_point2 = np.array([center[0] + r*np.cos(beta + alpha), center[1] + r*np.sin(beta + alpha)])
    res = [np.linalg.norm(point - end_point1), np.linalg.norm(point - end_point2)]
    # intersection of radius and arc
    t = r / np.linalg.norm(point - center)
    small_distance = np.array([SMALL_NUM, SMALL_NUM])
    intersec_point1 = np.array([center[0] + (point[0] - center[0])*t, center[1] + (point[1] - center[1])*t])
    intersec_point2 = np.array([center[0] - (point[0] - center[0])*t, center[1] - (point[1] - center[1])*t])
    u = point - center
    phi = arctan(u[1], u[0])
    if is_between_anglesHai(phi, beta, alpha):
        return (abs(t-1))*np.linalg.norm(point - center)
    if is_between_anglesHai(phi + np.pi, beta, alpha):
        return (abs(t+1))*np.linalg.norm(point - center)
    return min(res)

def closesDistanceBetweenArcsHai(center1, center2, beta1, beta2, r, alpha, code1=-1, code2=-1):
    # arc1
    endpoint1 = center1 + np.array([r * np.cos(beta1 - alpha), r * np.sin(beta1 - alpha)])
    endpoint2 = center1 + np.array([r * np.cos(beta1 + alpha), r * np.sin(beta1 + alpha)])
    # arc2
    endpoint3 = center2 + np.array([r * np.cos(beta2 - alpha), r * np.sin(beta2 - alpha)])
    endpoint4 = center2 + np.array([r * np.cos(beta2 + alpha), r * np.sin(beta2 + alpha)])
    res = [
        closest_distance_between_point2arc(endpoint1, center2, beta2, r, alpha),
        closest_distance_between_point2arc(endpoint2, center2, beta2, r, alpha),
        closest_distance_between_point2arc(endpoint3, center1, beta1, r, alpha),
        closest_distance_between_point2arc(endpoint4, center1, beta1, r, alpha)
    ]
    u = center2 - center1
    phi = arctan(u[1], u[0])
    if code1 == 2 and code2 == 30:
        print(phi, end="\t")
        print(beta1 - alpha, end="\t")
        print(beta1 + alpha, end="\t")
        print(beta2 - alpha, end="\t")
        print(beta2 + alpha)
        print(angle(phi), end="\t")
        print(angle(beta1 - alpha), end="\t")
        print(angle(beta1 + alpha), end="\t")
        print(angle(beta2 - alpha), end="\t")
        print(angle(beta2 + alpha))
        print(is_between_anglesHai(angle(np.pi + phi), beta2, alpha), is_between_anglesHai(phi, beta1, alpha))
    if is_between_anglesHai(angle(np.pi + phi), beta2, alpha) and is_between_anglesHai(phi, beta1, alpha):
        res.append(max(closest_distance_between_point2arc(center2, center1, beta1, r, alpha) - r, 0, code1, code2))
        res.append(max(closest_distance_between_point2arc(center1, center2, beta2, r, alpha) - r, 0, code2, code2))
    return min(res)


def minimum__sectors_distance(sensor1, sensor2, code1=-1, code2=-1):
    # TODO: thang hai lam not phan sensor cua may di, k/c 2 sensor = min(khoang cach canh-canh, canh-cung, cung-cung)
    resultarray = []
    a0 = np.array([sensor1.xi, sensor1.yi])
    a1 = np.array([sensor1.xi + sensor1.r * np.cos(sensor1.betai - sensor1.alpha), sensor1.yi+sensor1.r*np.sin(sensor1.betai - sensor1.alpha)])
    a2 = np.array([sensor1.xi + sensor1.r * np.cos(sensor1.betai + sensor1.alpha), sensor1.yi+sensor1.r*np.sin(sensor1.betai + sensor1.alpha)])
    b0 = np.array([sensor2.xi, sensor2.yi])
    b1 = np.array([sensor2.xi + sensor2.r * np.cos(sensor2.betai - sensor2.alpha), sensor2.yi + sensor2.r * np.sin(sensor2.betai - sensor2.alpha)])
    b2 = np.array([sensor2.xi + sensor2.r * np.cos(sensor2.betai + sensor2.alpha), sensor2.yi + sensor2.r * np.sin(sensor2.betai + sensor2.alpha)])
    resultarray.append(closestDistanceBetweenLinesHai(a0, a1, b0, b1))
    resultarray.append(closestDistanceBetweenLinesHai(a0, a1, b0, b2))
    resultarray.append(closestDistanceBetweenLinesHai(a0, a2, b0, b1))
    resultarray.append(closestDistanceBetweenLinesHai(a0, a2, b0, b2))
    resultarray.append(closesDistanceBetweenArcsHai(a0, b0, sensor1.betai, sensor2.betai, sensor1.r, sensor1.alpha, code1, code2))
    resultarray.append(closestDistanceBetweenArcAndLine(a0, a1, b0, sensor2.betai, sensor2.r, sensor2.alpha)[2])
    resultarray.append(closestDistanceBetweenArcAndLine(a0, a2, b0, sensor2.betai, sensor2.r, sensor2.alpha)[2])
    resultarray.append(closestDistanceBetweenArcAndLine(b0, b1, a0, sensor1.betai, sensor1.r, sensor1.alpha)[2])
    resultarray.append(closestDistanceBetweenArcAndLine(b0, b2, a0, sensor1.betai, sensor1.r, sensor1.alpha)[2])
    # print(resultarray)
    resultarray = np.array(resultarray)
    # print(resultarray)
    return np.min(resultarray)
    # return resultarray[np.argmin(resultarray[:, 2])]