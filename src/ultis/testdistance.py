import numpy as np
from src.sensors.sensors_field import Sensor, Sensors_field
from src.ultis import distance


def test(sen1, sen2):
    sens_field = Sensors_field(10, 10)
    sens_field.add_sensor(sen1)
    sens_field.add_sensor(sen2)
    res = distance.minimum__sectors_distance(sen1, sen2)
    print("Min: %f" % res[2])
    sens_field.add_dis_2_points([res[0], res[1]])
    sens_field.field_show()
if __name__ == '__main__':
    testcases =[]
    sen1 = Sensor(xi=3, yi=3, alpha=70.0 / 180 * np.pi, betai=0, r=2)

    sen2 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=0.0 / 180 * np.pi, r=2)
    sen3 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=15.0 / 180 * np.pi, r=2)
    sen4 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=-15.0 / 180 * np.pi, r=2)
    sen5 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=90.0 / 180 * np.pi, r=2)
    sen6 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=-90.0 / 180 * np.pi, r=2)
    sen7 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=120.0 / 180 * np.pi, r=2)
    sen8 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=225.0 / 180 * np.pi, r=2)

    sen9 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=0.0 / 180 * np.pi, r=2)
    sen10 = Sensor(xi=5, yi=3, alpha=70.0 / 180 * np.pi, betai=15.0 / 180 * np.pi, r=2)
    sen11 = Sensor(xi=4, yi=3, alpha=70.0 / 180 * np.pi, betai=-15.0 / 180 * np.pi, r=2)
    sen12 = Sensor(xi=3, yi=3, alpha=70.0 / 180 * np.pi, betai=90.0 / 180 * np.pi, r=2)
    sen13 = Sensor(xi=3, yi=7, alpha=70.0 / 180 * np.pi, betai=-90.0 / 180 * np.pi, r=2)
    sen14 = Sensor(xi=2, yi=6, alpha=70.0 / 180 * np.pi, betai=120.0 / 180 * np.pi, r=2)
    sen15 = Sensor(xi=2, yi=5, alpha=70.0 / 180 * np.pi, betai=225.0 / 180 * np.pi, r=2)

    case = [sen1, sen2]
    testcases.append(case)
    case = [sen1, sen3]
    testcases.append(case)
    case = [sen1, sen4]
    testcases.append(case)
    case = [sen1, sen5]
    testcases.append(case)
    case = [sen1, sen6]
    testcases.append(case)
    case = [sen1, sen7]
    testcases.append(case)
    case = [sen1, sen8]
    testcases.append(case)
    case = [sen1, sen9]
    testcases.append(case)
    case = [sen1, sen10]
    testcases.append(case)
    case = [sen1, sen11]
    testcases.append(case)
    case = [sen1, sen12]
    testcases.append(case)
    case = [sen1, sen13]
    testcases.append(case)
    case = [sen1, sen14]
    testcases.append(case)
    case = [sen1, sen15]
    testcases.append(case)


    for c in testcases:
        test(c[0], c[1])