import copy
import numpy as np
import random
from src.sensors.sensors_field import WBG
from heapq import *
np.random.seed(1234)
random.seed(1234)

MAX_VALUE = 9999


class Chromosome:
	def __init__(self, num_sensors, K, x=None, wbg=None):
		self.K = K  # sensor end destination point, not source point
		self.destination = num_sensors + 1
		self.num_sensors = num_sensors
		self.num_features = num_sensors + K
		self.wbg = wbg
		if x is None:
			self.x = np.minimum(np.arange(1, num_sensors + 1 + K), self.destination).astype(int)
			# self.x = np.array([np.arange(1, num_sensors+2) for i in range(0, K)])
			np.random.shuffle(self.x)
			# self.x[-1] -= 10
			# print(self.x)
		else:
			self.x = np.minimum(x, self.destination).astype(int)

	def cross_over_ox(self, chromosome, pc=0.5):
		# child1 = self.x
		passed_gen1 = set()
		passed_gen2 = set()
		des_pos_1 = np.where(self.x == self.destination)[0]
		if len(des_pos_1) != chromosome.K:
			print(self.x)
			print("******", des_pos_1, self.K, len(self.x))
		for i in range(0, self.K):
			self.x[des_pos_1[i]] += i

		des_pos_2 = np.where(chromosome.x == chromosome.destination)
		des_pos_2 = des_pos_2[0]
		if len(des_pos_2) != chromosome.K:
			print(chromosome.x)
			print("***", des_pos_2, chromosome.K, len(chromosome.x))
		for i in range(0, chromosome.K):
			chromosome.x[des_pos_2[i]] += i
		child1 = copy.deepcopy(self.x).astype(int)
		child2 = copy.deepcopy(chromosome.x).astype(int)

		for sensor_index in range(0, self.num_features):
			pass_rate_1, pass_rate_2 = np.random.uniform(0, 1, 2)#* (i+self.num_sensors)/num_sensors
			if pass_rate_1 < pc:
				# pass gen to child1
				passed_gen1.add(self.x[sensor_index])
			else:
				child1[sensor_index] = -1
			if pass_rate_2 < pc:
				# pass gen to child2
				passed_gen2.add(chromosome.x[sensor_index])
			else:
				child2[sensor_index] = -1
		index_children1, index_children2, index1_parent1, index1_parent2 = 0, 0, 0, 0
		while index_children1 < self.num_features:
			if child1[index_children1] < 0:
				while index1_parent2 < chromosome.num_features and chromosome.x[index1_parent2] in passed_gen1:
					index1_parent2 += 1
				if index1_parent2 < chromosome.num_features:
					child1[index_children1] = chromosome.x[index1_parent2]
				index1_parent2 += 1
			index_children1 += 1
		while index_children2 < self.num_features:
			if child2[index_children2] < 0:
				while index1_parent1 < self.num_features and self.x[index1_parent1] in passed_gen2:
					index1_parent1 += 1
				if index1_parent1 < self.num_features:
					child2[index_children2] = self.x[index1_parent1]
				index1_parent1 += 1
			index_children2 += 1
		# child1 = np.minimum(child1, self.destination)
		# child2 = np.minimum(child2, self.destination)

		self.x = np.minimum(self.x, self.destination).astype(int)
		chromosome.x = np.minimum(chromosome.x, chromosome.destination).astype(int)
		return Chromosome(num_sensors=self.num_sensors, K=self.K, x=child1, wbg=self.wbg), \
				Chromosome(num_sensors=self.num_sensors, K=self.K, x=child2, wbg=self.wbg)

	def local_mutation(self, p_swap=0.5):
		mutate_type = np.random.uniform(0, 1)
		rand1, rand2 = np.random.randint(0, int(self.num_sensors + self.K), 2)
		while rand1 == rand2:
			rand1, rand2 = np.random.randint(0, int(self.num_sensors + self.K), 2)
		if mutate_type < p_swap:
			self.x[rand1], self.x[rand2] = self.x[rand2], self.x[rand1]
			self.optimize_path()
		else:  # ko swap se insert
			if rand1 > rand2:
				rand1, rand2 = rand2, rand1
			if rand2 < self.num_features-1:
				self.x = np.concatenate(
					(self.x[0:rand1], [self.x[rand2]], self.x[rand1:rand2], self.x[rand2 + 1:])).astype(int)
			else:
				self.x = np.concatenate(
					(self.x[0:rand1], [self.x[rand2]], self.x[rand1:rand2])).astype(int)
			self.optimize_path()
		return None

	def convert_to_path(self):
		paths = []
		split_point = np.where(self.x == self.destination)[0]
		for i in range(0, self.K):
			start_index = 0 if i == 0 else split_point[i - 1] + 1
			path = [0] + [ver for ver in self.x[start_index:split_point[i]]] + [self.destination]
			paths.append(path)
		return paths

	def dijkstra(self, vertexs, f, t):
		q, seen = [(0, f)], set()
		mins = {x: MAX_VALUE for x in vertexs}
		mins[f] = 0
		prev = {}
		while q:
			(cost, v1) = heappop(q)
			if v1 not in seen:
				seen.add(v1)
				# path = (v1, path)
				if v1 == t:
					path = [t]
					x = prev[t]
					while x != f:
						path.append(x)
						x = prev[x]
					# print("path", path)
					return cost, path[::-1]
				for v2 in vertexs:
				# for c, v2 in self.wbg[v1]:
					if v2 in seen or v2 == v1:
						continue
					# cur_cost = cost + self.wbg.o_adj_matrix[v1)][v2]
					cur_cost = cost + self.wbg.o_adj_matrix[int(v1)][int(v2)]
					if prev is None or cur_cost <= mins[v2]:
						mins[v2] = cur_cost
						prev[v2] = v1
						heappush(q, (cur_cost, v2))
		return float("inf"), None

	def optimize_path(self):
		trash = []
		new_path = []
		split_point = np.where(self.x == self.destination)[0]
		for i in range(0, len(split_point)):
			start_index = 0 if i == 0 else split_point[i-1]+1
			vertexs = [0] + [ver for ver in self.x[start_index:split_point[i]]] + [self.destination]
			cost, path = self.dijkstra(vertexs, 0, self.destination)
			for x in vertexs:
				if x == 0:
					continue
				if x not in path:
					trash.append(x)
			new_path = new_path + path  # skip start point
		if self.x[-1] != self.destination:
			self.x = np.concatenate((new_path, trash, self.x[split_point[-1]+1:])).astype(int)
		else:
			self.x = np.concatenate((new_path, trash)).astype(int)
		return None

	def repair_path(self):

		pass


if __name__ == '__main__':
	num_sensors = 30
	num_barrier = 1
	wbg = WBG(lenght=50, height=15, mode='strong')
	wbg.create_sensors_randomly(num_sensor=num_sensors, r=2, alpha=np.pi * 4 / 5)
	wbg.build_WBG()
	wbg.o_adj_matrix[0][num_sensors + 1] = np.ceil(wbg.L / wbg.sensors_list[0].lr)
	wbg.o_adj_matrix[num_sensors + 1][0] = np.ceil(wbg.L / wbg.sensors_list[0].lr)

	print(wbg.o_adj_matrix)
	##############################################################################################
	test_chromosome1 = Chromosome(num_sensors=num_sensors, K=num_barrier, wbg=wbg,
								x=np.minimum(np.arange(1, num_sensors + 1 + num_barrier), num_sensors + 1))
	test_chromosome2 = Chromosome(num_sensors=num_sensors, K=num_barrier, x=None, wbg=wbg)
	print(test_chromosome1.x)
	print(test_chromosome1.optimize_path(), test_chromosome1.num_features)
	print("F1 1")
	print(test_chromosome1.x)

	print("############################################################################")
	print("F1 2")
	test_chromosome2.optimize_path()
	print(test_chromosome2.x)
	wbg.field_show()
