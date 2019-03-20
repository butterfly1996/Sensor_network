import numpy as np
import copy

class Population():
    def __init__(self, wbg, num_particles, K):
        print ('Init population ...')
        self.wbg = wbg
        self.num_bins = 50 # TODO: refactor
        self.particles = [Particle(self, wbg, K) for i in range(num_particles)]
        self.num_particles = num_particles
        self.crp = 0.5
        self.mup = 0.01
    def initialize(self):
        for i in range(len(self.particles)):
            self.particles[i].initialize()
        for i in range(len(self.particles)):
            self.particles[i].compute_density(self.num_bins)
    def rollete(self, id_list, prob_list, num_choices):
        return np.random.choice(id_list, size=num_choices,p=prob_list)
    def cross_over(self):
        argsort = np.argsort([particle.affinity() for particle in self.particles]) # asc
        parents = argsort[:int(len(argsort)/(10/3))] # get top 30%
        np.random.shuffle(parents)
        self.children = []
        for i in range(int(len(parents)/2)):
            child1, child2 = self.particles[parents[i]].cross_over(self.particles[parents[i+1]])
            child1.compute_density(self.num_bins)
            child2.compute_density(self.num_bins)
            self.particles.append(child1)
            self.particles.append(child2)
            self.children.append(child1)
            self.children.append(child2)
    def mutation(self):
        pass
    def selection(self):
        argsort = np.argsort([particle.fitness() for particle in self.particles])[:self.num_particles] # asc
        self.particles = [self.particles[i] for i in argsort]
    def evolve(self, MAX_ITER=500):
        self.initialize()
        for iter in range(MAX_ITER):
            #self.show()
            print('crossover...')
            self.cross_over();
            print('mutation... ')
            self.mutation()
            print('selection...')
            self.selection()
            print('Iter %d, fitness=%d'%(iter, self.particles[0].fitness()))

    def show(self):
        for i, particle in enumerate(self.particles):
            print('\nParticle %d, fitness %d,  affinity %f, density %f, prob %f:'%(i, particle.fitness(), particle.affinity(), particle.density(), particle.ps(self.num_bins)))
            #print(particle.chrome)
            #print(particle._paths)

class Particle():
    def __init__(self, population, wbg, K):
        self.population = population
        self.N = len(wbg.sensors_list)
        self.K = K
        self.wbg = wbg
        self.chrome = np.random.randint(1, 100, size=(self.N+1+self.K,)).astype(float)
    def initialize(self):
        # TODO: refactor setter
        self.compute_fitness()
        self.compute_affinity()
        # self.compute_density(self.population.num_bins)
        self.best_p = self.chrome.copy()
        self.best_fit = self.fitness()
    def get_population(self):
        return self.population
    def compute_fitness(self, verbose=False):
        mask = np.ones_like(self.chrome).astype(float) # mask values, 1 if vertex can be selected, else -np.inf
        mask[0] = -np.inf # s 
        paths = []
        result = 0
        for k in range(self.K):
            path = [0]
            next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, self.chrome)])
            if mask[next_v]==-np.inf: # invalid
                path.append(self.N+1) # direct barrier
                paths.append(path)
                continue
            while next_v <= self.N:
                mask[next_v] = -np.inf
                path.append(next_v)
                next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, self.chrome)]) # if negative? => BUG
            mask[next_v] = -np.inf
            path.append(self.N+1)
            paths.append(path)
        #print (paths)
        result = np.sum([self.wbg.length(path) for path in paths])
        if verbose:
            print (paths)
            print ('\nfitness: %d'%result)
        self._paths = paths
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
        point1, point2 = np.random.choice(np.arange(len(self.chrome)+1), 2, replace=False)
        point1, point2 = min(point1, point2), max(point1, point2)

        child1 = Particle(self.population, self.wbg, self.K)
        child2 = Particle(self.population, self.wbg, self.K)

        child1.chrome[0:point1] = self.chrome[0:point1]
        child1.chrome[point1:point2] = other_particle.chrome[point1:point2]
        child1.chrome[point2:] = self.chrome[point2:]
        child1.compute_fitness()
        child1.compute_affinity()
		
        child2.chrome[0:point1] = other_particle.chrome[0:point1]
        child2.chrome[point1:point2] = self.chrome[point1:point2]
        child2.chrome[point2:] = other_particle.chrome[point2:]
        child2.compute_fitness()
        child2.compute_affinity()

        return child1, child2
    def mutation(self):
        pass
