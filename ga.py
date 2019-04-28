import numpy as np
import copy
from scipy import stats
import time

class Population():
    def __init__(self, wbg, num_particles, K, omega, c1, c2):
        print ('Init population ...')
        self.wbg = wbg
        self.num_bins = 100 # TODO: refactor
        self.omega = omega
        self.c1 = c1
        self.c2 = c2
        self.num_particles = num_particles
        self.particles = [Particle(self, wbg, K, self.omega) for i in range(num_particles)]
        self.parents = []
        self.p_cross = 0.5
        self.p_mutation = 0.3
    def initialize(self):
        for i in range(len(self.particles)):
            self.particles[i].initialize()
    def evolve(self, MAX_ITER=500):
        self.initialize()
        print('START, fitness=%d' % (sorted(self.particles, key=lambda x: x.fitness())[0].fitness()))
        for iter in range(MAX_ITER):
            start = time.time()
            # cross_over
            print ('cross')
            cross_ids = np.random.choice(np.arange(len(self.particles)), int(self.p_cross*len(self.particles)), False)
            for i in range(0,len(cross_ids)//2*2,2):
                child1, child2 = self.particles[cross_ids[i]].cross_over(self.particles[cross_ids[i+1]])
                self.particles.append(child1)
                self.particles.append(child2)
            # mutation
            print ('mut')
            mutation_ids = np.random.choice(np.arange(len(self.particles)), int(self.p_mutation*len(self.particles)), False)
            for i in range(len(mutation_ids)):
                self.particles[mutation_ids[i]].mutation()
            # selection
            self.particles = sorted(self.particles, key=lambda x: x.fitness())[:self.num_particles]
            aTime = time.time() - start
            if (iter + 1) % 1 == 0:
                print('Iter %d, fitness=%d, time' % (iter, self.particles[0].fitness()), aTime)
        return self.particles[0].chrome, self.particles[0].fitness()

    def show(self):
        for i, particle in enumerate(self.particles):
            print('\nParticle %d, fitness %d,  affinity %f, density %f, prob %f:'%(i, particle.fitness(), particle.affinity(), particle.density(), particle.ps(self.num_bins)))
            print(particle.chrome)
            print(particle._paths)

class Particle():
    def __init__(self, population, wbg, K, omega):
        self.population = population
        self.N = len(wbg.sensors_list)
        self.K = K
        self.omega = omega
        self.wbg = wbg
        self.chrome = np.random.randint(-np.pi, np.pi, size=(self.N,)).astype(float)
    def __str__(self):
        return self.chrome
    def initialize(self):
        # TODO: refactor setter
        self.compute_fitness()
    def compute_fitness(self, verbose=False):
        self.wbg.update_views(self.chrome)
        result = self.wbg.loss()
        self._fitness = result
    def fitness(self):
        return self._fitness
    def cross_over(self, other_particle):
        child1 = Particle(self.population, self.wbg, self.K, self.omega)
        child1.chrome = self.chrome.copy()
        child2 = Particle(other_particle.population, other_particle.wbg, other_particle.K, other_particle.omega)
        child2.chrome = other_particle.chrome.copy()

        point1, point2 = np.random.choice(np.arange(len(child1.chrome)+1), 2, replace=False)
        point1, point2 = min(point1, point2), max(point1, point2)
        temp = copy.deepcopy(child1.chrome[point1:point2])
        child1.chrome[point1:point2] = child2.chrome[point1:point2]
        child2.chrome[point1:point2] = temp

        # chrome1 = self.chrome.copy()
        # chrome2 = other_particle.chrome.copy()
        # self.chrome = [np.random.choice([chrome1[i], chrome2[i]]) for i in range(len(self.chrome))]
        # other_particle.chrome = [np.random.choice([chrome1[i], chrome2[i]]) for i in range(len(self.chrome))]

        child1.compute_fitness()
        child2.compute_fitness()

        return child1, child2
    def mutation(self):
        self.chrome += np.random.normal(loc=0.0, scale=np.pi/3, size=(self.N,))
        self.compute_fitness()