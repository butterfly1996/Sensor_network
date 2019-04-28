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
    def initialize(self):
        for i in range(len(self.particles)):
            self.particles[i].initialize()
        # for i in range(len(self.particles)):
        #     self.particles[i].compute_density(self.num_bins)
        self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
        self.best_fit = np.min([particle.fitness() for particle in self.particles])
    def rollete(self, id_list, prob_list, num_choices):
        return np.random.choice(id_list, size=num_choices,p=prob_list, replace=True)
    def clone(self, num_chosen_immune):
        argsort = np.argsort([particle.affinity() for particle in self.particles])[::-1] # desc
        immune_lib = argsort[:int(len(argsort)/5)] # get top 20%
        chosen_immune = np.random.choice(immune_lib, min(num_chosen_immune, len(immune_lib)), replace=False)
        antibodies = list(chosen_immune)+list(argsort)
        probs = [self.particles[i].ps(self.num_bins) for i in antibodies]
        next_gen_ids = self.rollete(antibodies, prob_list=np.array(probs)/np.sum(probs), num_choices=len(argsort))
        self.particles = [self.particles[i] for i in next_gen_ids]
    def cross_over(self):
        argsort = np.argsort([particle.affinity() for particle in self.particles]) # asc
        parents = argsort[:int(len(argsort)/(10/3))] # get top 30%
        self.parents = parents # will be used in update()
        np.random.shuffle(parents)
        for i in range(int(len(parents)/2)):
            self.particles[parents[i]], self.particles[parents[i+1]] = self.particles[parents[i]].cross_over(self.particles[parents[i+1]])
        for i in range(len(self.particles)):
            self.particles[i].compute_density(self.num_bins)
    def update(self):
        for i in range(len(self.particles)):
            self.particles[i].update()
        if np.min([particle.fitness() for particle in self.particles]) < self.best_fit:
            self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
            self.best_fit = np.min([particle.fitness() for particle in self.particles])
    def evolve(self, MAX_ITER=500):
        self.initialize()
        print('START, fitness=%d'%(self.best_fit))
        for iter in range(MAX_ITER):
            #self.show()
            #print('clone...')
            #self.clone(10)
            #print('cross over...')
            #self.cross_over()
            #print('update...')
            start = time.time()
            self.update()
            aTime = time.time()-start
            if (iter+1)%1==0:
                print('Iter %d, fitness=%d, time'%(iter, self.best_fit), aTime)
        return self.best_fit, self.best_g



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
        #self.phi = np.random.randint(1, 100, size=(self.N,)).astype(float)
        # self.rand1 = self.gen_rand()
        # self.rand2 = self.gen_rand()
        self.velocity = np.zeros_like(self.chrome).astype(float)
        #self.velocity = np.random.uniform(-100, 100, size=(self.N+1+self.K,)).astype(float)
        self.c1 = self.population.c1
        self.c2 = self.population.c2
    def __str__(self):
        return self.chrome
    def initialize(self):
        # TODO: refactor setter
        self.compute_fitness()
        # self.compute_affinity()
        # self.compute_density(self.population.num_bins)

        self.best_p = self.chrome.copy()
        self.best_fit = self.fitness()
    def gen_rand(self):
        res = np.random.rand(1)[0]
        return res if res not in [0.25,0.5,0.75] else gen_rand()
    def get_population(self):
        return self.population
    def compute_fitness(self, verbose=False):
        self.wbg.update_views(self.chrome)
        result = self.wbg.loss()
        self._fitness = result
    def fitness(self):
        return self._fitness
    def compute_affinity(self):
        self._affinity = 1.0/(1+self.fitness())
    def affinity(self):
        return self._affinity
    def compute_density(self, num_bins):
        epsilon = 0.001
        max_aff = np.max([particle.affinity() for particle in self.population.particles])+epsilon
        min_aff = np.min([particle.affinity() for particle in self.population.particles])
        bin_size = (max_aff-min_aff)/num_bins
        bin_id = np.floor(self.affinity() / bin_size)
        self._density = np.sum([bin_id*bin_size<=particle.affinity()<(bin_id+1)*bin_size for particle in self.population.particles])
    def density(self):
        return self._density
    def ps(self, num_bins, alpha=0.7):
        return alpha*self.pf()+(1-alpha)*self.pd(num_bins)
    def pf(self):
        return 1.0-self.fitness()/(np.sum(particle.fitness() for particle in self.population.particles))
    def pd(self, num_bins):
        return 1.0-self.density()/(np.sum([particle.density() for particle in self.population.particles]))
    def cross_over(self, other_particle):
        '''
        point1, point2 = np.random.choice(np.arange(len(self.chrome)+1), 2, replace=False)
        point1, point2 = min(point1, point2), max(point1, point2)
        temp = copy.deepcopy(self.chrome[point1:point2])
        self.chrome[point1:point2] = other_particle.chrome[point1:point2]
        other_particle.chrome[point1:point2] = temp
        '''
        chrome1 = self.chrome.copy()
        chrome2 = other_particle.chrome.copy()
        self.chrome = [np.random.choice([chrome1[i], chrome2[i]]) for i in range(len(self.chrome))]
        other_particle.chrome = [np.random.choice([chrome1[i], chrome2[i]]) for i in range(len(self.chrome))]

        # TODO: refactor setter
        self.compute_fitness()
        self.compute_affinity()
        # self.compute_density(self.population.num_bins)

        other_particle.compute_fitness()
        other_particle.compute_affinity()
        # other_particle.compute_density(self.population.num_bins)

        return self, other_particle
    def mutation(self):
        self.rand1 = 4*self.rand1*(1-self.rand1)
        self.rand2 = 4*self.rand2*(1-self.rand2)
    def update(self, mutation=False):
        if mutation:
            self.mutation()
        else:
            self.rand1 = np.random.rand(1)[0]
            self.rand2 = np.random.rand(1)[0]
        self.velocity = self.omega*self.velocity + self.c1*self.rand1*(self.best_p-self.chrome)+self.c2*self.rand2*(self.population.best_g-self.chrome)
        self.chrome += self.velocity # TODO: refactor setter
        # self.chrome = np.clip(self.chrome, 1.0, 100.0)
        self.compute_fitness()
        if self.fitness() < self.best_fit:
            self.best_p = self.chrome.copy()
            self.best_fit = self.fitness()