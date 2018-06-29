import numpy as np
from sensors_field import Sensor, Sensors_field
import distance

def test(sen1, sen2):
    sens_field = Sensors_field()
    sens_field.add_sensor(sen1)
    sens_field.add_sensor(sen2)
    sens_field.field_show()
if __name__ == '__main__':
    testcases =[]
    sen1 = Sensor(xi=3, yi=3, alpha=70.0 / 180 * np.pi, betai=0, r=2)
    sen2 = Sensor(xi=8, yi=3, alpha=70.0 / 180 * np.pi, betai=0.0 / 180 * np.pi, r=2)
    case = [sen1, sen2]
    testcases.append(case)
    for c in testcases:
        test(c[0], c[1])