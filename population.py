import numpy as np
from particle import Particle

class Population:
    def __init__(self, wbg, num_particles, K):
        print ('Init population ...')
        self.wbg = wbg
        self.num_bins = 50 # TODO: refactor
        self.c1 = 0.2
        self.c2 = 0.4
        self.omega = 0.95
        self.particles = [Particle(self, wbg, K) for i in range(num_particles)]
        self.best_fit = np.inf

    def initialize(self):
        for i in range(len(self.particles)):
            self.particles[i].initialize()
        for i in range(len(self.particles)):
            self.particles[i].compute_density(self.num_bins)
        self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
        self.best_fit = np.min([particle.fitness() for particle in self.particles])

    def rollete(self, id_list, prob_list, num_choices):
        return np.random.choice(id_list, size=num_choices,p=prob_list)

    def clone(self, num_chosen_immune):
        argsort = np.argsort([particle.affinity() for particle in self.particles])[::-1] # desc
        immune_lib = argsort[:int(len(argsort)/5)] # get top 20%
        chosen_immune = np.random.choice(immune_lib, min(num_chosen_immune, len(immune_lib)), replace=False)
        antibodies = list(chosen_immune)+list(argsort)
        probs = [self.particles[i].ps(self.num_bins) for i in antibodies]
        next_gen_ids = self.rollete(antibodies, prob_list=np.array(probs)/np.sum(probs), num_choices=self.num_bins)
        self.particles = [self.particles[i] for i in next_gen_ids]
        # print("num particies: ", len(self.particles))

    def cross_over(self):
        """
        lai ghep
        :return:
        """
        argsort = np.argsort([particle.affinity() for particle in self.particles]) # asc
        parents = argsort[:int(len(argsort)*0.3)] # get top 30%
        np.random.shuffle(parents)
        self.children = []
        for i in range(len(parents)//2):
            child1, child2 = self.particles[parents[2*i]].cross_over(self.particles[parents[2*i+1]])
            child1.compute_density(self.num_bins)
            child2.compute_density(self.num_bins)
            self.particles.append(child1)
            self.particles.append(child2)
            self.children.append(len(self.particles)-1)
            self.children.append(len(self.particles)-2)

    def update(self):
        for i in range(len(self.particles)):
            if i in self.children: # mutation after crossover
                self.particles[i].update(mutation=True)
            else:
                self.particles[i].update()
        for i in range(len(self.particles)):
            self.particles[i].compute_density(self.num_bins)
        tmp_fit = np.min([particle.fitness() for particle in self.particles])
        if tmp_fit < self.best_fit:
            self.best_g = self.particles[np.argmin([particle.fitness() for particle in self.particles])].chrome.copy()
            self.best_fit = tmp_fit

    def evolve(self, MAX_ITER=500):
        self.initialize()
        for iter in range(MAX_ITER):
            # self.show()
            # print('clone...')
            self.clone(50)
            # print('cross over...')
            self.cross_over()
            # print('update...')
            self.update()
            # if (iter+1)%10==0:
            #     print('Iter %d, fitness=%d'%(iter, self.best_fit))


    def show(self):
        for i, particle in enumerate(self.particles):
            print('\nParticle %d, fitness %d,  affinity %f, density %f, prob %f:'%(i, particle.fitness(), particle.affinity(), particle.density(), particle.ps(self.num_bins)))
            print(particle.chrome)
            print(particle._paths)