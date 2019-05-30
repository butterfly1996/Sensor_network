import logging
import numpy as np
import multiprocessing
import ctypes
from igraph import *
from time import time
import click

from dbg import DBG
import global_vars


logging.basicConfig(filename='log.txt', level=logging.DEBUG)
global_vars.initialize()

MODE_EXIST_BARRIER = 0 # bai toan xac dinh su ton tai bao phu rao chan (1a)
MODE_MIN_COST_BARRIER = 1 # bai toan xay dung mot hang rao voi chi phi toi thieu (1b)
MODE_MAX_NUM_BARRIER = 2 # bai toan toi da so rao chan dat duoc (2)

@click.command()
@click.option('--mode', default=MODE_EXIST_BARRIER)
def run_exp(mode=MODE_EXIST_BARRIER):
    if mode == MODE_EXIST_BARRIER:
        for n in np.arange(10, 90, 10):
            NUM_SIM = 100
            simulation = 0

            r = 50
            rate = 0.0
            while simulation < NUM_SIM:  # or count < 10:
                simulation += 1
                start = time()

                global_vars.dbg = DBG(lenght=500, height=100)
                global_vars.dbg.create_sensors_randomly(num_sensor=n, r=r, alpha=np.pi / 6)
                global_vars.dbg.build_virtual_nodes()

                shared_array_base = multiprocessing.Array(ctypes.c_double, [np.inf] * global_vars.dbg.num_virtuals ** 2)
                shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
                global_vars.misfit = shared_array.reshape(global_vars.dbg.num_virtuals, global_vars.dbg.num_virtuals)
                global_vars.dbg.build_dbg()

                res = global_vars.dbg.dfs2(0, global_vars.dbg.adj_matrix.shape[0] - 1)

                print('time ', time() - start, ', rate ', res)
                rate += res
            rate /= NUM_SIM
            print('rate n=%d : ' % (n), rate)
            logging.info('rate n=%d: %f' % (n, rate))
    elif mode == MODE_MIN_COST_BARRIER:
        for r in np.arange(20, 30, 10):
            NUM_SIM = 100
            simulation = 0
            count = 0

            n = 40
            rate = 0.0
            while simulation < NUM_SIM or count < NUM_SIM:
                count += 1
                start = time()

                global_vars.dbg = DBG(lenght=500, height=100)
                global_vars.dbg.create_sensors_randomly(num_sensor=n, r=r, alpha=np.pi / 6)
                global_vars.dbg.build_virtual_nodes()

                shared_array_base = multiprocessing.Array(ctypes.c_double, [np.inf] * global_vars.dbg.num_virtuals ** 2)
                shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
                global_vars.misfit = shared_array.reshape(global_vars.dbg.num_virtuals, global_vars.dbg.num_virtuals)
                global_vars.dbg.build_dbg()

                res = global_vars.dbg.dijkstra2(0, global_vars.dbg.adj_matrix.shape[0] - 1)
                if res != np.inf:
                    rate += res
                    simulation += 1
                    print(res)
                else:
                    print('bugg')
                    continue

                print('time ', time() - start, ', total angle ', res)
            rate /= simulation
            print('total angle r=%d : ' % (r), rate)
            logging.info('total angle r=%d: %f' % (n, rate))
    elif mode == MODE_MAX_NUM_BARRIER:
        for n in np.arange(10, 50, 10):
            NUM_SIM = 50
            simulation = 0

            r = 30
            rate = 0.0
            while simulation < NUM_SIM:  # or count < 10:
                simulation += 1
                start = time()

                global_vars.dbg = DBG(lenght=100, height=50)
                global_vars.dbg.create_sensors_randomly(num_sensor=n, r=r, alpha=np.pi / 3)
                global_vars.dbg.build_virtual_nodes()

                shared_array_base = multiprocessing.Array(ctypes.c_double, [np.inf] * global_vars.dbg.num_virtuals ** 2)
                shared_array = np.ctypeslib.as_array(shared_array_base.get_obj())
                global_vars.misfit = shared_array.reshape(global_vars.dbg.num_virtuals, global_vars.dbg.num_virtuals)
                global_vars.dbg.build_dbg()

                res = global_vars.dbg.max_flow()

                print('time ', time() - start, ', num bar ', res)
                rate += res
            rate /= NUM_SIM
            print('num bar n=%d : ' % (n), rate)
            logging.info('num bar n=%d: %f' % (n, rate))

if __name__ == '__main__':
    run_exp()