import numpy as np
import copy


class Particle:
	def __init__(self, population, wbg, K):
		self.population = population
		self.N = len(wbg.sensors_list)
		self.K = K
		self.wbg = wbg
		# self.chrome = np.random.randint(1, 100, size=(self.N + 1 + self.K,)).astype(float)
		x = [np.random.uniform(0, 100) for i in range(0, self.N+1)] + [np.random.uniform(0, 30) for i in range(0, self.K)]
		self.chrome = np.array(x)
		self.rand1 = self.gen_rand()
		self.rand2 = self.gen_rand()
		self.mutation_rate = 2
		self.velocity = np.random.uniform(0, 1, self.N+1+self.K)#np.random.randn(self.N + 1 + self.K)  #
		self.c1 = self.population.c1
		self.c2 = self.population.c2
		# print("c1: ", self.c1, "\tc2: ", self.c2)
		self.omega = self.population.omega

	def initialize(self):
		# TODO: refactor setter
		self.compute_fitness()
		self.compute_affinity()
		# self.compute_density(self.population.num_bins)

		self.best_p = self.chrome.copy()
		self.best_fit = self.fitness()

	def gen_rand(self):
		res = np.random.rand(1)[0]
		return res  # if res not in [0.25,0.5,0.75] else self.gen_rand()

	def get_population(self):
		return self.population

	def compute_fitness(self, verbose=False):
		mask = np.ones_like(self.chrome).astype(float)  # mask values, 1 if vertex can be selected, else -np.inf
		mask[0] = -np.inf  # s
		paths = []
		result = 0
		for k in range(self.K):
			path = [0]
			next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, self.chrome)])
			if mask[next_v] == -np.inf:  # invalid
				path.append(self.N + 1)  # direct barrier
				paths.append(path)
				continue
			while next_v <= self.N and mask[next_v] != -np.inf:
				mask[next_v] = -np.inf
				path.append(next_v)
				next_v = np.argmax([ch if m == 1 else m for m, ch in zip(mask, self.chrome)])  # if negative? => BUG
			mask[next_v] = -np.inf
			path.append(self.N + 1)
			# print(path)
			paths.append(path)
		# print ([self.wbg.length(path) for path in paths])
		result = np.sum([self.wbg.length(path) for path in paths])
		# print(paths, result)
		if verbose:
			print(paths)
			# print('\nfitness: %d' % result)
		self._paths = paths
		self._fitness = result

	def fitness(self):
		return self._fitness

	def compute_affinity(self):
		self._affinity = 1.0 / (1 + self.fitness())

	def affinity(self):
		return self._affinity

	def compute_density(self, num_bins):
		epsilon = 0.0001
		max_aff = np.max([particle.affinity() for particle in self.population.particles]) + epsilon
		min_aff = np.min([particle.affinity() for particle in self.population.particles])
		bin_size = (max_aff - min_aff) / num_bins
		bin_id = np.floor(self.affinity() / bin_size)
		self._density = np.sum([bin_id * bin_size <= particle.affinity() < (bin_id + 1) * bin_size for particle in
		                        self.population.particles])

	def density(self):
		return self._density

	def ps(self, num_bins, alpha=0.7):
		return alpha * self.pf() + (1 - alpha) * self.pd(num_bins)

	def pf(self):
		return 1.0 - self.fitness() / (np.sum(particle.fitness() for particle in self.population.particles))

	def pd(self, num_bins):
		return 1.0 - self.density() / (np.sum([particle.density() for particle in self.population.particles]))

	def cross_over(self, other_particle):
		point1, point2 = np.random.choice(np.arange(len(self.chrome) + 1), 2, replace=False)
		point1, point2 = min(point1, point2), max(point1, point2)

		child1 = Particle(self.population, self.wbg, self.K)
		child2 = Particle(self.population, self.wbg, self.K)

		child1.chrome[0:point1] = self.chrome[0:point1]
		child1.chrome[point1:point2] = other_particle.chrome[point1:point2]
		child1.chrome[point2:] = self.chrome[point2:]
		child1.velocity = (self.velocity+other_particle.velocity)*1/2
		child1.initialize()

		child2.chrome[0:point1] = other_particle.chrome[0:point1]
		child2.chrome[point1:point2] = self.chrome[point1:point2]
		child2.chrome[point2:] = other_particle.chrome[point2:]
		child2.velocity = (self.velocity + other_particle.velocity) * 1 / 2
		child2.initialize()

		return child1, child2

	def mutation(self):
		self.rand1 = 4 * self.rand1 * (1 - self.rand1)
		self.rand2 = 4 * self.rand2 * (1 - self.rand2)
		# # print("N: ", self.N)
		# # x = np.arange(self.N+2)
		# # np.random.shuffle(x)

		# mutation_vec = [(1) for i in range(0, self.N+1)]\
		#                         + [-0.5 for i in range(0, self.K)]
		# mutation_vec = (np.array(mutation_vec)) * np.random.uniform(0, 1, self.N + 1 + self.K)
		# # print(mutation_vec * self.mutation_rate)
		# self.velocity *= mutation_vec

	def update(self, mutation=False):
		if mutation:
			self.mutation()
		else:
			self.rand1 = np.random.uniform(0, 1)
			self.rand2 = np.random.uniform(0, 1)
			self.velocity = 0.5 * self.velocity + self.c1 * self.rand1 * (
					self.best_p - self.chrome) + self.c2 * self.rand2 * (self.population.best_g - self.chrome)
		# print("dm sinh oc cho, ", np.linalg.norm(self.velocity), np.linalg.norm((self.best_p-self.chrome))
		#       , np.linalg.norm((self.population.best_g-self.chrome)))
		# self.velocity = self.omega*self.velocity + self.c1*self.rand1*(self.best_p-self.chrome)+self.c2*self.rand2*(self.population.best_g-self.chrome)
		# self.chrome += self.velocity # TODO: refactor setter
		# self.velocity = self.omega * self.velocity + self.c1 * self.rand1 * (
		# 		self.best_p - self.chrome) + self.c2 * self.rand2 * (self.population.best_g - self.chrome)

		self.chrome += self.omega * self.velocity  # * self.chrome # TODO: refactor setter

		# DEBUG:
		# print("debug: ")
		# print(self.chrome)

		self.compute_fitness()
		self.compute_affinity()
		# self.compute_density(self.population.num_bins)
		if self.fitness() < self.best_fit:
			self.best_p = self.chrome.copy()
			self.best_fit = self.fitness()