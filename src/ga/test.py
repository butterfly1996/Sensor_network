# 02373822922
import heapq
import sys
import numpy as np
from src.sensors.sensors_field import WBG
from src.ga.hai_individuals_priority_base import Chromosome
import copy
np.random.seed(1234)


def initialize(num_sensors, num_barrier, num_individuals):
	return [Chromosome(num_sensors=num_sensors, K=num_barrier) for i in range(0, num_individuals)]


def optimize(path):
	predecessor, Q, visited_nodes = {}, [], set()
	s = path[0]
	t = path[-1]
	d = {s: 0}  # vertex -> minimal distance
	Qd = {}  # vertex -> [d[v], parent_v, v]
	p = {}  # predecessor
	visited_set = {path[0]}

	for v in path:
		if v == s:
			continue
		d[v] = wbg.o_adj_matrix[s][v]
		item = [d[v], s, v]
		heapq.heappush(Q, item)
		Qd[v] = item

	while Q:
		# print(Q)
		cost, parent, u = heapq.heappop(Q)
		if u not in visited_set:
			# print('visit:', u)
			p[u] = parent
			visited_set.add(u)
			if u == t:
				return p, d[u]
			for v in path:
				if d.get(v):
					if d[v] > wbg.o_adj_matrix[u][v] + d[u]:
						d[v] = wbg.o_adj_matrix[u][v] + d[u]
						Qd[v][0] = d[v]  # decrease key
						Qd[v][1] = u  # update predecessor
						heapq._siftdown(Q, 0, Q.index(Qd[v]))
				else:
					d[v] = wbg.o_adj_matrix[u][v] + d[u]
					item = [d[v], u, v]
					heapq.heappush(Q, item)
					Qd[v] = item

	return None


def compute_fitness(chromosome, K, verbose=False):
	# print(chromosome.x)
	# print(K)
	paths = chromosome.convert_to_path()
	result = np.sum([wbg.length(path) for path in paths])
	return result


def crossover(population):
	np.random.shuffle(population)
	loop = len(population)//2
	for i in range(loop):
		# child1, child2 = population[2*i].cross_over_ox(population[2*i+1])
		# child1, child2 = population[2*i].cross_over_hai_position_base(population[2*i+1], Pc=0.4)
		###################################################################################################
		# child1, child2 = population[2*i].cross_over_hai_priority_base(population[2*i+1], Pc=0.4, P_skip=0.4)
		# if child1 is not None:
		# 	population.append(child1)
		# if child2 is not None:
		# 	population.append(child2)
		child1, child2 = population[2 * i].cross_over_ox(population[2 * i + 1], pc=0.1)
		# child1, child2 = population[2 * i].cross_over_hai_priority_base(population[2 * i + 1], Pc=0.3, P_skip=0.6)
		if child1 is not None:
			population.append(child1)
		if child2 is not None:
			population.append(child2)
	# print(population)
	# print(children)
	# population = [population+ children]
	return population


def mutation(population, mutate_rate):
	# pop_size = len(population)

	for chromosome in population:
		rand1 = np.random.uniform(0, 1)
		if rand1 < mutate_rate:
			x = copy.deepcopy(chromosome)
			x.local_mutation()
			population.append(x)
	return population


def selection(population, K, pop_size, verbose=False):
	prob_list = [compute_fitness(chromosome, K) for chromosome in population]
	if verbose:
		# print(prob_list)
		arg_min = np.argmin(prob_list)
		fitness_list = [compute_fitness(chromosome, K, verbose=verbose) for chromosome in population]
		print("best fitness: ", min(fitness_list))
		print(population[arg_min].x)
		print(population[arg_min].convert_to_path())
	rank_list = np.argsort(prob_list)
	new_population = [population[id] for id in rank_list[:pop_size]]
	# prob_list = normalize_prob(prob_list)
	# population = np.random.choice(population, size=pop_size, replace=False, p=prob_list).tolist()

	# population = sorted(population, key=lambda x: compute_fitness(x, K))
	# return population
	return new_population


ITER = 201
EPSILON = 1e-5
POP_SIZE = 100


def normalize_prob(prob_list):
	prob_list = np.array(prob_list)
	prob_list = 1 / (np.array(prob_list))
	# print("sum: ", np.sum(prob_list))
	# pivot = np.mean(prob_list, axis=0)*0.9
	# for i in range(0, len(prob_list)):
	# 	if prob_list[i] < pivot:
	# 		prob_list[i] = EPSILON
	prob_list = prob_list / np.sum(prob_list)
	return prob_list


def main():
	from src.ultis import distance
	global wbg
	num_sensors = 30
	num_barrier = 1
	wbg = WBG(lenght=50, height=10, mode='strong')
	wbg.create_sensors_randomly(num_sensor=num_sensors, r=2, alpha=np.pi * 4 / 5)
	wbg.build_WBG()
	wbg.o_adj_matrix[0][num_sensors+1] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.o_adj_matrix[num_sensors+1][0] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.show_matrix()
	population = initialize(num_sensors=num_sensors, num_barrier=num_barrier, num_individuals=POP_SIZE)
	# for z in population:
	# 	print(z.x)
	# sys.exit(0)
	# z = Chromosome(num_sensors=num_sensors,
	#                K=1,
	#                x=np.array([[1, 2, 3, 4, 25, 6, 26, 8, 9, 10,
	#                             11, 12, 30, 14, 24, 16, 17, 27, 19, 20,
	#                             21, 29, 22, 5, 31, 15, 7, 28, 13, 18, 23]]))
	# print(compute_fitness(z, 1, True))
	# population.append(z)
	# sys.exit(0)
	for iter in range(ITER):
		if iter % 10 == 0:
			print("Iter: ", iter)
		population = crossover(population)
		# print("crossover", len(population))
		population = mutation(population, 0.1)
		# print("mutation", len(population))
		population = selection(population, K=num_barrier, pop_size=POP_SIZE, verbose=iter % 10 == 0)
		# print("seletion", len(population))
	Pk, Nm = wbg.min_num_mobile_greedy(num_barrier)
	print("######################################################")
	print(Pk)
	print(Nm)
	print("######################################################")

	print(wbg.o_adj_matrix[0][24], wbg.o_adj_matrix[24][12], wbg.o_adj_matrix[12][21], wbg.o_adj_matrix[21][17],
	      wbg.o_adj_matrix[17][6], wbg.o_adj_matrix[6][4], wbg.o_adj_matrix[4][14],wbg.o_adj_matrix[14][31])
	wbg.field_show()
	pass


def test():
	p = optimize([0, 1, 6, 17, 31])
	print(p)


if __name__ == '__main__':
	main()
	# test()