import copy
import numpy as np


class Chromosome:
	def __init__(self, num_sensors, K, x=None):
		self.num_features = num_sensors + K  # sensor end destination point, not source point
		self.destination = num_sensors + 1
		if x is None:
			self.x = np.array([min(i, self.destination) for i in range(1, num_sensors + 1 + K)])
			np.random.shuffle(self.x)

		else:
			self.x = x

	def cross_over_ox(self, chromosome):
		# print("x: ", len(self.x))
		# print("x: ", self.x)
		cut_point1, cut_point2 = np.random.randint(0, len(self.x), 2)
		if cut_point2 == cut_point1:
			return None, None
		if cut_point1 > cut_point2:
			cut_point1, cut_point2 = cut_point2, cut_point1
		# print("log crossover: ", cut_point1, cut_point2)

		child1 = np.concatenate(
			([self.destination for i in range(0, cut_point1)],
			 self.x[cut_point1:cut_point2],
			 [self.destination for i in range(cut_point2, self.num_features)]), axis=0)
		child2 = np.concatenate(
			([self.destination for i in range(0, cut_point1)],
			 chromosome.x[cut_point1:cut_point2],
			 [self.destination for i in range(cut_point2, self.num_features)]), axis=0)
		i = cut_point1
		j1 = cut_point2
		j2 = cut_point2
		if chromosome.x[i] not in child1:
			child1[j1] = chromosome.x[i]
			j1 = (j1 + 1) % (self.num_features)
		if self.x[i] not in child2:
			# print(j2, i, len(child1), len(child2))
			child2[j2] = self.x[i]
			j2 = (j2 + 1) % (self.num_features)
		i = (i + 1) % (self.num_features)
		while i != int(cut_point1):
			if chromosome.x[i] not in child1:
				child1[j1] = chromosome.x[i]
				j1 = (j1 + 1) % (self.num_features)
			if self.x[i] not in child2:
				child2[j2] = self.x[i]
				j2 = (j2 + 1) % (self.num_features)
			i = (i + 1) % (self.num_features)
		# print(child1)
		# print(child2)
		return Chromosome(num_sensors=self.destination - 1, K=self.num_features - self.destination + 1, x=child1), \
		       Chromosome(num_sensors=self.destination - 1, K=self.num_features - self.destination + 1, x=child2)


def local_mutation(self):
		pos1 = np.random.randint(0, self.num_features)
		if pos1 >= self.num_features-2:
			return None
		pos2 = np.random.randint(pos1+1, len(self.x)-1)
		self.x = np.concatenate((self.x[0:pos1], [self.x[pos2]], self.x[pos1:pos2], self.x[pos2+1:]), axis=0)
		# self.x[pos1], self.x[pos2] = self.x[pos2], self.x[pos1]


def global_mutation(self):
		des_list = np.where(self.x == self.destination)
		tmp_array = np.split(self.x, des_list[0])
		tmp_array = np.concatenate(tmp_array, axis=0)
		des_list = np.where(tmp_array == self.destination)
		if des_list is not None and len(des_list) > 0:
			first_des = des_list[0][0]
			# print("first_des", first_des)
			if first_des==0 or first_des+1 >= self.num_features:
				return None
			pos1, pos2 = np.random.randint(0, first_des), np.random.randint(first_des+1, self.num_features)
			try:
				next_des = np.where(first_des[0] >= pos2)[0][0]
			except:
				next_des = self.num_features
			pos3 = np.random.randint(pos2, next_des)
			x = np.concatenate((tmp_array[0:pos1], tmp_array[pos2:pos3], tmp_array[pos1+1:pos2],
			                    [tmp_array[pos1]], tmp_array[pos3:]),
			                   axis=0)
			# print("global_muation", pos1, pos2, pos3)
			# print(tmp_array)
			# print(x)
			return Chromosome(num_sensors=self.destination - 1, K=self.num_features - self.destination + 1, x=x)
		return None

