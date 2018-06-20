import numpy as np
from dkmthangsinhcanba import Sensor
def closestDistanceBetweenLines(a0,a1,b0,b1,clampAll=True,clampA0=True,clampA1=True,clampB0=True,clampB1=True):

    ''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
        Return the closest points on each segment and their distance
    '''

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
                        return a0,b0,np.linalg.norm(a0-b0)
                    return a0,b1,np.linalg.norm(a0-b1)


            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1,b0,np.linalg.norm(a1-b0)
                    return a1,b1,np.linalg.norm(a1-b1)


        # Segments overlap, return distance between parallel segments
        return None,None,np.linalg.norm(((d0*_A)+a0)-b0)



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


    return pA,pB,np.linalg.norm(pA-pB)

def closesDistanceBetweenArcs(center1, center2, beta1, beta2, r, alpha):
    res = []
    # arc1
    endpoint1 = center1 + np.array([r*np.cos(beta1+alpha), r*np.sin(beta1+alpha)])
    endpoint2 = center1 + np.array([r * np.cos(beta1 - alpha), r * np.sin(beta1 - alpha)])
    # arc2
    endpoint3 = center2 + np.array([r * np.cos(beta2 + alpha), r * np.sin(beta2 + alpha)])
    endpoint4 = center2 + np.array([r * np.cos(beta2 - alpha), r * np.sin(beta2 - alpha)])
    #1: endpoints of two arcs
    res.append(np.linalg.norm(endpoint1-endpoint3))
    res.append(np.linalg.norm(endpoint1-endpoint4))
    res.append(np.linalg.norm(endpoint2-endpoint3))
    res.append(np.linalg.norm(endpoint2-endpoint4))
    #2: endpoint to interior
    phi = np.arctan((endpoint3[1]-center1[1])/(endpoint3[0]-center1[0]))
    if phi <= beta1+alpha and phi >= beta1-alpha:
        res.append(np.linalg.norm(endpoint3-np.array([center1[0]+r*np.cos(phi), center1[1]+r*np.sin(phi)])))

    phi = np.arctan((endpoint4[1] - center1[1]) / (endpoint4[0] - center1[0]))
    if phi <= beta1 + alpha and phi >= beta1 - alpha:
        res.append(np.linalg.norm(endpoint4 - np.array([center1[0] + r * np.cos(phi), center1[1] + r * np.sin(phi)])))

    phi = np.arctan((endpoint1[1] - center2[1]) / (endpoint1[0] - center2[0]))
    if phi <= beta1 + alpha and phi >= beta1 - alpha:
        res.append(np.linalg.norm(endpoint1 - np.array([center1[0] + r * np.cos(phi), center1[1] + r * np.sin(phi)])))

    phi = np.arctan((endpoint2[1] - center2[1]) / (endpoint2[0] - center2[0]))
    if phi <= beta1 + alpha and phi >= beta1 - alpha:
        res.append(np.linalg.norm(endpoint2 - np.array([center1[0] + r * np.cos(phi), center1[1] + r * np.sin(phi)])))
    #3: interiors
    phi = np.arctan((center2[1]-center1[1])/(center2[0]-center1[0]))
    if beta1-alpha <= phi <= beta1+alpha and beta2-alpha <= phi <= beta2+alpha:
        res.append(np.linalg.norm(center1+r*np.array([np.cos(phi), np.sin(phi)])-center2-r*np.array([np.cos(phi-np.pi), np.sin(phi-np.pi)])))
    #4: intersection
    a = 2*r*(center2[1]-center1[1])
    b = 2*r*(center2[0]-center1[0])
    d = np.linalg.norm(center2-center1)
    c = -d**2
    if b**2*c**2 >= (a**2+b**2)*(c**2-a**2):
        res.append(0)
    return np.min(np.array(res))

def closestDistanceBetweenArcAndLine(a0, a1, center, beta, r, alpha):
    ### a0 a1 la 2 dau mut cua doan thang
    ### center beta r alpha la thong so cua day cung
    res = []
    # ends
    endpoint1 = center + np.array([r * np.cos(beta + alpha), r * np.sin(beta + alpha)])
    endpoint2 = center + np.array([r * np.cos(beta - alpha), r * np.sin(beta - alpha)])
    res.append(np.linalg.norm(a0, endpoint1))
    res.append(np.linalg.norm(a0, endpoint2))
    res.append(np.linalg.norm(a1, endpoint1))
    res.append(np.linalg.norm(a1, endpoint2))
    # interiors
    #x = (a0*(a1-b1)**2+c0*(b0-a0)**2+(a1-b1)*(b0-a0)*(c1-a1)) / ((a1-b1)**2+(a0-b0)**2)
    #y = c1-(c0-x)(b0-a0)/(a1-b1)

    x = (a0[0] * (a0[1] - a1[1]) ** 2 + center[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (center[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = center[1] - (center[0] - x)(a1[0] - a0[0]) / (a0[1] - a1[1])
    d = np.linalg.norm(np.array([x, y]) - center) - r
    if d > 0:
        res.append(d)
    # end of arc to interior of line
    x = (a0[0] * (a0[1] - a1[1]) ** 2 + endpoint1[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (
    endpoint1[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = endpoint1[1] - (endpoint1[0] - x)(a1[0] - a0[0]) / (a0[1] - a1[1])
    res.append(np.linalg.norm(np.array([x, y]) - endpoint1))

    x = (a0[0] * (a0[1] - a1[1]) ** 2 + endpoint2[0] * (a1[0] - a0[0]) ** 2 + (a0[1] - a1[1]) * (a1[0] - a0[0]) * (
        endpoint2[1] - a0[1])) / (
            (a0[1] - a1[1]) ** 2 + (a0[0] - a1[0]) ** 2)
    y = endpoint2[1] - (endpoint2[0] - x)(a1[0] - a0[0]) / (a0[1] - a1[1])
    res.append(np.linalg.norm(np.array([x, y]) - endpoint2))
    # end of line to interior of arc
    d = np.linalg.norm(a0 - center) - r
    if d > 0:
        res.append(d)
    d = np.linalg.norm(a1 - center) - r
    if d > 0:
        res.append(d)
    # intersection
    A = 1+(a1[1]-a0[1])**2/(a1[0]-a0[0])**2
    B = -2*center[0]-(2*a0[0]*(a1[1]-a0[1])**2+2*(a1[1]-a0[1])*(a0[1]-center[1])*(a1[0]-a0[0]))/(a1[0]-a0[0])**2
    C = center[0]*2+(a0[0]**2*(a1[1]-a0[1])**2+(a0[1]-center[1])**2*(a1[0]-a0[0])**2-2*a0[0]*(a1[1]-a0[1])*(a0[1]-center[1])*(a1[0]-a0[0]))/(a1[0]-a0[0])**2
    delta = B**2-4*A*C
    if delta >= 0:
        res.append(0)
    return np.min(np.array(res))

def minimum__sectors_distance(sensor1, sensor2):
    resultarray = []
    a0 = np.array(sensor1.xi, sensor1.yi)
    a1 = np.array(sensor1.xi + sensor1.r * np.cos(sensor1.betai - sensor1.alpha), sensor2.yi+sensor2.r*np.sin(sensor2.betai - sensor2.alpha))
    b0 = np.array(sensor2.xi, sensor2.yi)
    b1 = np.array(sensor2.xi + sensor2.r * np.cos(sensor2.betai - sensor2.alpha), sensor2.yi+sensor2.r*np.sin(sensor2.betai - sensor2.alpha))
    result = closestDistanceBetweenLines(a0, a1, b0, b1)
    resultarray.append(result)
    ############################################################################################################3
    a0 = np.array(sensor1.xi, sensor1.yi)
    a1 = np.array(sensor1.xi + sensor1.r * np.cos(sensor1.betai - sensor1.alpha),
                  sensor2.yi + sensor2.r * np.sin(sensor2.betai - sensor2.alpha))
    b0 = np.array(sensor2.xi, sensor2.yi)
    b1 = np.array(sensor2.xi + sensor2.r * np.cos(sensor2.betai - sensor2.alpha),
                  sensor2.yi + sensor2.r * np.sin(sensor2.betai - sensor2.alpha))
    result = closestDistanceBetweenLines(a0, a1, b0, b1)
    resultarray.append(result)
    ############################################################################################################3
    pass