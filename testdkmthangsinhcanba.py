#from pyswarm import pso
from psopy import minimize
import numpy as np
from sensors_field import Sensor
from sensors_field import WBG

N = 80
K = 3
def banana(chrome):
    print ("call")
    wbg.update_views(chrome)
    result = wbg.loss()
    if False: #verbose:
        print ('\nfitness: %d'%result)
    return result

import  distance
global distance
wbg = WBG(lenght=500, height=100, mode='strong')
# sensor_field.create_sensors_randomly(num_sensor=sensor_field.n, r=3, alpha=60)
#s1 = Sensor(3, 3, np.pi / 5, 2, np.pi / 4)
#s2 = Sensor(8, 3, 5*np.pi / 6, 2, np.pi / 4)
wbg.create_sensors_randomly(N, r=50, alpha=np.pi/3)
wbg.build_disBG()
wbg.setup_loss()
print (wbg.dis_matrix)

wbg.field_show()
lb = [-np.pi]*(N)
ub = [np.pi]*(N)

x0 = np.random.uniform(-np.pi, np.pi, (30,N))

#xopt, fopt = pso(banana, lb, ub, swarmsize=50, maxiter=10, debug=True)
import time
start = time.time()
res = minimize(banana, x0, options={'max_iter':50, 'stable_iter':10, 'verbose':True})
print ("time: ", (time.time()-start))
xopt, fopt = res.x, res.fun
print (xopt)
print (fopt)
wbg.update_views(xopt)
wbg.field_show()

