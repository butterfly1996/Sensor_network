from src.sensors.sensors_field import Sensors_field

if __name__ == '__main__':
    sensor_field = Sensors_field(lenght=30, height=10, num_station=2, num_mobile=20)
    sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
    sensor_field.field_show()

