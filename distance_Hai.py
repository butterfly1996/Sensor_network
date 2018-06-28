import numpy as np
from sympy.solvers import solve
from sympy import Symbol
from sensors_field import Sensor, Sensors_field
import time
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
    #r=s=R, beta+alpha>u, v>beta-alpha
    #u=v
    u = Symbol('u')
    solution = solve(-2*x*sympy.sin(u)+2*y*sympy.cos(u))
    print(solution)
    #u=pi+2phi-v
    u = Symbol('u')
    solution = solve(-2 * x * sympy.sin(u) + 2 * y * sympy.cos(u) - 2 * sen1.r *sympy.sin(2*u - 2*phi))
    print(solution)

    return result
def get_min_distance(sen1, sen2):
    x = sen1.xi-sen2.xi
    y = sen1.yi - sen2.yi
    result = get_min_dis_line2line(sen1, sen2, x, y)
    get_min_dis_arg2arg(sen1, sen2, x, y)
    return result

if __name__ == '__main__':
    print("running ....")
    sen1 = Sensor(xi=3, yi=3, alpha=70.0/180*np.pi, betai=0, r=2)
    sen2 = Sensor(xi=8, yi=3, alpha=70.0/180*np.pi, betai=180.0/180*np.pi, r=2)
    print(get_min_distance(sen1, sen2))
    print(np.sin(np.pi/2))
    sensor_fields = Sensors_field(lenght=10, height=10)
    sensor_fields.add_sensor(sen1)
    sensor_fields.add_sensor(sen2)
    sensor_fields.field_show()
