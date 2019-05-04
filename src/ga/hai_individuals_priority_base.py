import copy
import numpy as np
import random
from src.sensors.sensors_field import WBG
np.random.seed(1234)
random.seed(1234)


class Chromosome:
	def __init__(self, num_sensors, K, x=None):
		self.K = K  # sensor end destination point, not source point
		self.destination = num_sensors + 1
		self.num_sensors = num_sensors
		if x is None:
			# self.x = np.random.uniform(0, K, (K, num_sensors+1))
			self.x = np.array([np.arange(1, num_sensors+2) + np.random.uniform(-1, 1, num_sensors+1) for i in range(0, K)])
			# self.x = np.array([np.arange(1, num_sensors+2) for i in range(0, K)])
			for path in self.x:
				np.random.shuffle(path)
			# self.x[-1] -= 10
			# print(self.x)
		else:
			self.x = x

	def cross_over_hai_position_base(self, chromosome, Pc=0.8):
		# parent1 = self.x
		parent1 = copy.deepcopy(self.x)
		parent2 = copy.deepcopy(chromosome.x)
		for i in range(0, self.num_sensors+1):
			rand1, rand2 = np.random.uniform(0, 1, 2) #* (i+self.num_sensors)/num_sensors
			if rand1 < Pc:
				# print(rand1, rand2)
				if rand2 < Pc:
					z = copy.deepcopy(parent1[:, i])
					parent1[:, i] = copy.deepcopy(parent2[:, i])
					parent2[:, i] = z
				else:
					parent1[:, i] = copy.deepcopy(parent2[:, i])
			elif rand2 < Pc:
				# print(rand1, rand2)
				parent2[:, i] = copy.deepcopy(parent1[:, i])
		return Chromosome(num_sensors=self.destination - 1, K=self.K, x=parent1), \
		       Chromosome(num_sensors=self.destination - 1, K=self.K, x=parent2)

	def cross_over_hai_priority_base(self, chromosome, Pc=0.6, P_skip=0.4):
		# parent1 = self.x
		parent1 = copy.deepcopy(self.x)
		parent2 = chromosome.x
		i1, i2 = 0, 0
		while i1 < self.K and i2 < self.K:
			rand_swap, rand_skip1, rand_skip2 = np.random.uniform(0, 1, 3) #* (i+self.num_sensors)/num_sensors
			# print(rand1, rand2)

			if rand_swap < Pc:
				z = parent1[i1]
				parent1[i1] = parent2[i2]
				parent2[i2] = z
				i1 += 1
				i2 += 1
			else:
				if rand_skip1 < P_skip:
					i1 += 1
				if rand_skip2 < P_skip:
					i2 += 1
		return Chromosome(num_sensors=self.destination - 1, K=self.K, x=parent1), \
		       Chromosome(num_sensors=self.destination - 1, K=self.K, x=parent2)

	def local_mutation(self, include_destination=1):
		# print("mutation:", self.x)
		# rand1 = np.argmax(
		# 	np.random.uniform(0, 1, self.num_sensors+include_destination)*self.x[0])
		# rand2 = np.random.randint(0, self.num_sensors + include_destination)
		np.random.shuffle(self.x)
		for i in self.x:
			rand1, rand2 = np.random.randint(0, self.num_sensors + include_destination, 2)
			while self.x[i][rand1] < self.x[i][-1] and self.x[i][rand2] < self.x[i][-1]:
				rand1, rand2 = np.random.randint(0, self.num_sensors + include_destination, 2)
			# print(rand1, rand2)
			# z = copy.deepcopy(self.x[:, rand1])
			# self.x[:, rand1] = self.x[:, rand2]
			# self.x[:, rand2] = z
			z = copy.deepcopy(self.x[i, rand1])
			self.x[i, rand1] = self.x[i, rand2]
			self.x[i, rand2] = z
		return self.x

	def convert_to_path(self):
		paths = []
		for i in range(0, self.K):
			path = [0]
			args = np.argsort(self.x[i])[::-1]
			len_args = self.destination
			for sensor in range(0, len_args):
				if args[sensor] == self.num_sensors:# num_sensor = destination - 1
					path.append(self.destination)
					paths.append(path)
					break
				elif np.argmax(self.x[:, args[sensor]]) == i:
					path.append(args[sensor]+1)
		# print(paths)
		return paths


if __name__ == '__main__':
	num_sensors = 6
	num_barrier = 3
	wbg = WBG(lenght=50, height=15, mode='strong')
	wbg.create_sensors_randomly(num_sensor=num_sensors, r=2, alpha=np.pi * 4 / 5)
	wbg.build_WBG()
	wbg.o_adj_matrix[0][num_sensors + 1] = np.ceil(wbg.L / wbg.sensors_list[0].lr)
	wbg.o_adj_matrix[num_sensors + 1][0] = np.ceil(wbg.L / wbg.sensors_list[0].lr)
	##############################################################################################
	test_chromosome1 = Chromosome(num_sensors=num_sensors, K=num_barrier, x=None)
	test_chromosome2 = Chromosome(num_sensors=num_sensors, K=num_barrier, x=None)
	# child1, child2 = test_chromosome1.cross_over_hai_position_base(test_chromosome2)
	print("############################################################################")
	print(test_chromosome1.x)
	print("############################################################################")
	test_chromosome1.local_mutation()
	print(test_chromosome1.x)
	print("############################################################################")
	paths = test_chromosome1.convert_to_path()
	print(paths)
	# print("chromosome1: ", test_chromosome1.x)
	# # test_chromosome1.local_mutation()
	# # print(test_chromosome1.x)
	# child1, child2 = test_chromosome1.cross_over_hai(test_chromosome2)
	# # print("###############################################################")
	# print("test_chromosome1", test_chromosome1.x)
	# print("child1: ", child1.x)
	# # print("chromosome2: ", test_chromosome2.x)
	# # print("test_chromosome2", test_chromosome2.x)
	# # print("child2: ", child2.x)