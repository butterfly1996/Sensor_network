import numpy as np
from src.sensors.sensors_field import Sensor, Sensors_field
from sympy.solvers import solve
from sympy import Symbol
import sympy
def get_min_dis_line2line(sen1, sen2, x, y):
    result = []
    #R>r>0, u=beta+alpha, R>s>0, v=beta+alpha
    u = sen1.betai+sen1.alpha
    v = sen2.betai+sen2.alpha
    A=np.array([[2, -2*np.cos(u-v)],
                [-2*np.cos(u-v), 2]])
    if np.linalg.det(A) != 0:
        B = np.array([-2*x*np.cos(u)-2*y*np.sin(u),
                      2 * x * np.cos(v) + 2 * y * np.sin(v)])
        solution = np.linalg.solve(A, B)
        # print(solution)
        result.append([[sen1.xi+solution[0]*np.cos(u), sen1.yi+solution[0]*np.sin(u)],
                   [sen2.xi+solution[1] * np.cos(v), sen2.yi+solution[1] * np.sin(v)]])
    ##################################################################
    # R>r>0, u=beta-alpha, R>s>0, v=beta+alpha
    u = sen1.betai - sen1.alpha
    v = sen2.betai + sen2.alpha
    A = np.array([[2, -2 * np.cos(u - v)],
                  [-2 * np.cos(u - v), 2]])
    if np.linalg.det(A) != 0:
        B = np.array([-2 * x * np.cos(u) - 2 * y * np.sin(u),
                      2 * x * np.cos(v) + 2 * y * np.sin(v)])
        solution = np.linalg.solve(A, B)
        # print(solution)
        result.append([[sen1.xi + solution[0] * np.cos(u), sen1.yi + solution[0] * np.sin(u)],
                   [sen2.xi + solution[1] * np.cos(v), sen2.yi + solution[1] * np.sin(v)]])
    ##################################################################
    # R>r>0, u=beta+alpha, R>s>0, v=beta-alpha
    u = sen1.betai + sen1.alpha
    v = sen2.betai - sen2.alpha
    A = np.array([[2, -2 * np.cos(u - v)],
                  [-2 * np.cos(u - v), 2]])
    if np.linalg.det(A) != 0:
        B = np.array([-2 * x * np.cos(u) - 2 * y * np.sin(u),
                      2 * x * np.cos(v) + 2 * y * np.sin(v)])
        solution = np.linalg.solve(A, B)
        # print(solution)
        result.append([[sen1.xi + solution[0] * np.cos(u), sen1.yi + solution[0] * np.sin(u)],
                       [sen2.xi + solution[1] * np.cos(v), sen2.yi + solution[1] * np.sin(v)]])
    ##################################################################
    # R>r>0, u=beta-alpha, R>s>0, v=beta-alpha
    u = sen1.betai - sen1.alpha
    v = sen2.betai - sen2.alpha
    A = np.array([[2, -2 * np.cos(u - v)],
                  [-2 * np.cos(u - v), 2]])
    if np.linalg.det(A) != 0:
        B = np.array([-2 * x * np.cos(u) - 2 * y * np.sin(u),
                      2 * x * np.cos(v) + 2 * y * np.sin(v)])
        solution = np.linalg.solve(A, B)
        # print(solution)
        result.append([[sen1.xi + solution[0] * np.cos(u), sen1.yi + solution[0] * np.sin(u)],
                       [sen2.xi + solution[1] * np.cos(v), sen2.yi + solution[1] * np.sin(v)]])
    ##################################################################
    return result
#cong thuc chi dung trong truong hop 2 sensor chung ban kinh
def get_min_dis_arg2arg(sen1, sen2, x, y):
    result = []
    phi = np.arctan(y/x)
    modul = np.sqrt(x**2+y**2)
    if(np.abs(np.sin(phi)-x/modul)>0.0001 and np.abs(np.cos(phi)-y/modul)>0.0001):
        phi+=np.pi
    if (np.abs(np.sin(phi) - x / modul) > 0.0001 and np.abs(np.cos(phi) - y / modul) < 0.0001):
        phi = -phi
    if (np.abs(np.sin(phi) - x / modul) < 0.0001 and np.abs(np.cos(phi) - y / modul) > 0.0001):
        phi = np.pi-phi
    print("phi::: ", phi)
    #r=s=R, beta+alpha>u, v>beta-alpha
    #u=v
    u = Symbol('u')
    solution = solve(-2*x*sympy.sin(u)+2*y*sympy.cos(u))
    print(solution)
    #u=pi+2phi-v
    # u = Symbol('u')
    # solution = solve(-2 * x * sympy.sin(u) + 2 * y * sympy.cos(u) - 2 * sen1.r *sympy.sin(2*u - 2*phi))
    solution = [phi, phi-np.pi]
    if modul<2*sen1.r:
        solution.append(np.arccos(-modul/2/sen1.r)+phi)
        solution.append(-np.arccos(-modul/2/sen1.r)+phi)
    print(solution)
    return result
def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step
def get_min_dis_arg2line(sen1, sen2, x, y, num = 10000):
    # s=R, u=beta+alpha
    u = sen1.betai+sen1.alpha
    c = x*np.cos(u)+y*np.sin(u)
    print("beta: ", sen2.betai, "\talpha: ", sen2.alpha)
    bound_up = sen2.betai+sen2.alpha
    bound_low = sen2.betai - sen2.alpha
    result = []
    # v = Symbol('v')
    v_values = np.linspace(bound_low, bound_up, num)

    # print(v_values)
    solution = [x*np.sin(v)-y*np.cos(v)-(sen2.r*np.cos(u-v)-c)*np.sin(u-v)
                for v in v_values]
    # plt.plot(solution, v_values)
    # plt.show()
    index = 0
    # result = 0
    solution = np.abs(np.array(solution))
    min = solution[0]
    # print(min)
    # print(solution)
    # print(len(solution))

    for i in range(0, len(solution)):
        if (min > solution[i]):
            min = solution[i]
            index=i
            # print(i)
        # if(solution[i]<0.0001):
        #     print(i, "::", solution)
    v = v_values[index]
    r = sen2.r*np.cos(u-v)-c
    result.append([[sen1.xi+r*np.cos(u), sen1.yi+r*np.sin(u)], [sen2.xi+sen2.r*np.cos(v), sen2.yi+sen2.r*np.sin(v)]])
    #############################################################################################################
    # s=R, u=beta-alpha
    u = sen1.betai - sen1.alpha
    c = x * np.cos(u) + y * np.sin(u)
    # bound_up = sen2.betai + sen2.alpha
    # bound_low = sen2.betai - sen2.alpha
    # # v = Symbol('v')
    # v_values = np.linspace(bound_low, bound_up, num)

    # print(v_values)
    solution = [x * np.sin(v) - y * np.cos(v) - (sen2.r * np.cos(u - v) - c) * np.sin(u - v)
                for v in v_values]
    # plt.plot(solution, v_values)
    # plt.show()
    index = 0
    solution = np.abs(np.array(solution))
    min = solution[0]
    # print(min)
    # print(solution)
    # print(len(solution))

    for i in range(0, len(solution)):
        if (min > solution[i]):
            min = solution[i]
            index = i
            # print(i)
        # if (solution[i] < 0.0001):
        #     print(i, "::", solution)
    v = v_values[index]
    r = sen2.r * np.cos(u - v) - c
    result.append([[sen1.xi + r * np.cos(u), sen1.yi + r * np.sin(u)],
                   [sen2.xi + sen2.r * np.cos(v), sen2.yi + sen2.r * np.sin(v)]])
    ###################################################################################################################
    # r=R, v=beta+alpha
    v = sen2.betai + sen2.alpha
    c = x * np.cos(v) + y * np.sin(v)
    bound_up = sen1.betai + sen1.alpha
    bound_low = sen1.betai - sen1.alpha
    u_values = np.linspace(bound_low, bound_up, num)

    # print(v_values)
    solution = [-x * np.sin(u) + y * np.cos(u) + (sen1.r * np.cos(u - v) + c) * np.sin(u - v)
                for u in u_values]
    # plt.plot(solution, v_values)
    # plt.show()
    index = 0
    solution = np.abs(np.array(solution))
    min = solution[0]

    for i in range(0, len(solution)):
        if (min > solution[i]):
            min = solution[i]
            index = i
            # print(i)
        # if (solution[i] < 0.0001):
        #     print(i, "::", solution)
    u = u_values[index]
    s = sen1.r * np.cos(u - v) + c
    result.append([[sen1.xi + sen1.r * np.cos(u), sen1.yi + sen1.r * np.sin(u)],
                   [sen2.xi + s * np.cos(v), sen2.yi + s * np.sin(v)]])
    ###################################################################################################################
    # r=R, v=beta-alpha
    v = sen2.betai - sen2.alpha
    c = x * np.cos(v) + y * np.sin(v)
    bound_up = sen1.betai + sen1.alpha
    bound_low = sen1.betai - sen1.alpha
    u_values = np.linspace(bound_low, bound_up, num)

    # print(v_values)
    solution = [-x * np.sin(u) + y * np.cos(u) + (sen1.r * np.cos(u - v) + c) * np.sin(u - v)
                for u in u_values]
    # plt.plot(solution, v_values)
    # plt.show()
    index = 0
    solution = np.abs(np.array(solution))
    min = solution[0]
    # print(min)
    # print(solution)
    # print(len(solution))

    for i in range(0, len(solution)):
        if (min > solution[i]):
            min = solution[i]
            index = i
            # print(i)
        # if (solution[i] < 0.0001):
        #     print(i, "::", solution)
    u = u_values[index]
    s = sen1.r * np.cos(u - v) + c
    result.append([[sen1.xi + sen1.r * np.cos(u), sen1.yi + sen1.r * np.sin(u)],
                   [sen2.xi + s * np.cos(v), sen2.yi + s * np.sin(v)]])
    ###################################################################################################################
    return result
    # print("bound", bound_low, ":::", bound_up)
    # print("solution::: ", solution[index], ":::", bound_low+epsilon*index, ":::", index)
def get_min_dis_point2line():
    
    pass
def get_min_distance(sen1, sen2):
    x = sen1.xi-sen2.xi
    y = sen1.yi - sen2.yi
    result = get_min_dis_line2line(sen1, sen2, x, y)
    # get_min_dis_arg2arg(sen1, sen2, x, y)
    # print("arg2line: ", get_min_dis_arg2line(sen1, sen2, x, y))
    result = get_min_dis_arg2line(sen1, sen2, x, y, num=5000)
    return result

if __name__ == '__main__':
    print("running ....")
    # sen1 = Sensor(xi=3, yi=3, alpha=70.0/180*np.pi, betai=0.0/180*np.pi, r=2)
    # sen2 = Sensor(xi=8, yi=3, alpha=70.0/180*np.pi, betai=0.0/180*np.pi, r=2)
    sen1 = Sensor(3, 3, np.pi/3, 2, np.pi/4)
    sen2 = Sensor(8, 3, np.pi/2, 2, np.pi/4)

    poinits = get_min_distance(sen1, sen2)
    print(poinits)
    print("plot")
    sensor_fields = Sensors_field(lenght=10, height=10)
    sensor_fields.add_sensor(sen1)
    sensor_fields.add_sensor(sen2)
    sensor_fields.add_dis_2_points(poinits[0])
    sensor_fields.add_dis_2_points(poinits[1])
    sensor_fields.add_dis_2_points(poinits[2])
    sensor_fields.add_dis_2_points(poinits[3])
    sensor_fields.field_show()
