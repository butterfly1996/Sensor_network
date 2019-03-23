import numpy as np
import copy

class Population():
    def __init__(self, wbg, num_particles, K, omega, c1, c2):
        print ('Init population ...')
        self.wbg = wbg
        self.c1 = c1
        self.c2 = c2
        self.G0 = 100.0
        self.beta = 20.0
        self.T = 500.0
        self.particles = [Particle(self, wbg, K) for i in range(num_particles)]
        self.num_particles = num_particles
    def initialize(self):
        for i in range(len(self.particles)):
            self.particles[i].initialize()
        
        self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
        self.best_fit = np.min([particle.fitness() for particle in self.particles])
    def G(self, t):
        return self.G0*np.exp(-1.0*self.beta*t/self.T)
    def compute_best(self):
        self._best = np.min([particle.fitness() for particle in self.particles])
    def compute_worst(self):
        self._worst = np.max([particle.fitness() for particle in self.particles])
    def best(self):
        return self._best
    def worst(self):
        return self._worst
    def forceTo(self, particle_id, t):
        pid = particle_id
        return np.sum([np.random.rand()*self.particles[pid].forceBy(self.particles[i], t) for i in range(self.num_particles) if i != pid])
    def calculate_acc(self, particle_id, t):
        #print ("m %f %f"%(self._best, self._worst))
        #print ([particle.m() for particle in self.particles])
        self.particles[particle_id].acc = self.forceTo(particle_id, t)/self.particles[particle_id].mass()
    def update(self, t):
        for i in range(len(self.particles)):
            self.compute_best(); self.compute_worst()
            self.calculate_acc(i, t)
            self.particles[i].update()
        if np.min([particle.fitness() for particle in self.particles]) < self.best_fit:
            self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
            self.best_fit = np.min([particle.fitness() for particle in self.particles])
    def evolve(self, MAX_ITER=500):
        self.initialize()
        for iter in range(MAX_ITER):
            #self.show()
            self.update(iter)
            if (iter+1)%50==0:
                print('Iter %d, fitness=%d'%(iter, self.best_fit))



    def show(self):
        for i, particle in enumerate(self.particles):
            print('\nParticle %d, fitness %d,  affinity %f, density %f, prob %f:'%(i, particle.fitness(), particle.affinity(), particle.density(), particle.ps(self.num_bins)))
            print(particle.chrome)
            print(particle._paths)

class Particle():
    def __init__(self, population, wbg, K):
        self.population = population
        self.N = len(wbg.sensors_list)
        self.K = K
        self.wbg = wbg
        self.chrome = np.random.randint(1, 100, size=(self.N+1+self.K,)).astype(float)
        self.rand1 = self.gen_rand()
        self.rand2 = self.gen_rand()
        #self.velocity = np.random.randn(self.N+1+self.K)
        self.velocity = np.zeros_like(self.N+1+self.K)
        self.acc = np.zeros_like(self.velocity)
        self.c1 = self.population.c1
        self.c2 = self.population.c2
    def initialize(self):
        # TODO: refactor setter
        self.compute_fitness()

        self.best_p = self.chrome.copy()
        self.best_fit = self.fitness()
    def gen_rand(self):
        res = np.random.rand()
        return res
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
    def m(self):
        return (self.fitness()-self.population.worst())/(self.population.best()-self.population.worst())
    def mass(self):
        return self.m()/np.sum([particle.m() for particle in self.population.particles])
    def forceBy(self, other_particle, t):
        return self.population.G(t)*self.mass()*other_particle.mass()/(np.linalg.norm(self.chrome-other_particle.chrome)+0.001)*(other_particle.chrome-self.chrome)
    def mutation(self):
        lower_threshold = 1.0
        upper_threshold = 100.0
        temp = self.chrome[(self.chrome > upper_threshold) | (self.chrome < lower_threshold)]
        self.chrome[(self.chrome > upper_threshold) | (self.chrome < lower_threshold)] = np.random.rand(len(temp))*(upper_threshold-lower_threshold)+lower_threshold
    def update(self):
        self.rand1 = np.random.rand()
        self.rand2 = np.random.rand()
        self.velocity = np.random.rand()*self.velocity + self.c1*self.rand1*(self.best_p-self.chrome)+self.c2*self.rand2*(self.population.best_g-self.chrome)+self.acc
        self.chrome += self.velocity # TODO: refactor setter
        # DEBUG:
        #print("debug: ")
        #print(self.chrome)
        self.mutation()
        self.compute_fitness()

        if self.fitness() < self.best_fit:
            self.best_p = self.chrome.copy()
            self.best_fit = self.fitness()