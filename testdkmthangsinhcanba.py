#!/usr/bin/env python3
import itertools
import multiprocessing
import defs

def func(params):
    a = params[0]
    b = params[1]
    #print (a*b*c*d)
    return a*b

if __name__ == '__main__':
    #Generate values for each parameter
    a = range(1,5)
    b = range(1,5)

    #Generate a list of tuples where each tuple is a combination of parameters.
    #The list will contain all possible combinations of parameters.
    paramlist = list(itertools.product(a,b))

    #Generate processes equal to the number of cores
    pool = multiprocessing.Pool(4)

    #Distribute the parameter sets evenly across the cores
    res  = pool.map(func,paramlist)
    print (res[:10])