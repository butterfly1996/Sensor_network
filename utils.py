import numpy as np


def angle(value):
    # chuan hoa gia tri goc tu -pi den pi
    return (value + np.pi) % (2 * np.pi) - np.pi

def is_between_angles(value, bisector, angle1, angle2):
    angle_min = min(angle1, angle2)
    angle_max = max(angle1, angle2)
    if angle_max - angle_min < np.pi:
        if np.isclose(bisector, (angle_max + angle_min) / 2):
            return angle_min <= value <= angle_max
        else:
            return not angle_min <= value <= angle_max
    elif value < 0:
        if bisector < 0:
            if np.isclose(bisector, (angle_max - 2 * np.pi + angle_min) / 2):
                return angle_max - 2 * np.pi <= value <= angle_min
            else:
                return not angle_max - 2 * np.pi <= value <= angle_min
        else:
            if np.isclose(bisector - 2 * np.pi, (angle_max - 2 * np.pi + angle_min) / 2):
                return angle_max - 2 * np.pi <= value <= angle_min
            else:
                return not angle_max - 2 * np.pi <= value <= angle_min
    else:
        if bisector < 0:
            if np.isclose(bisector, (angle_max - 2 * np.pi + angle_min) / 2):
                return angle_max - 2 * np.pi <= value - 2 * np.pi <= angle_min
            else:
                return not angle_max - 2 * np.pi <= value - 2 * np.pi <= angle_min
        else:
            if np.isclose(bisector - 2 * np.pi, (angle_max - 2 * np.pi + angle_min) / 2):
                return angle_max - 2 * np.pi <= value - 2 * np.pi <= angle_min
            else:
                return not angle_max - 2 * np.pi <= value - 2 * np.pi <= angle_min