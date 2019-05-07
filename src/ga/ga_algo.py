import numpy as np
from src.sensors.sensors_field import WBG
from src.ga.individuals import Chromosome
import copy
np.random.seed(1234)


def initialize(num_sensors, num_barrier, num_individuals):
	return [Chromosome(num_sensors=num_sensors, K=num_barrier) for i in range(0, num_individuals)]


def compute_fitness(chromosome, K, verbose=False):
	# print(chromosome.x)
	# print(K)
	paths = []
	i = 0
	index = 0
	path = [0]
	while i < K:
		if chromosome.x[index] != chromosome.destination:
			path.append(chromosome.x[index])
		else:
			path.append(chromosome.x[index])
			paths.append(path)
			path = [0]
			i += 1
		index += 1
	path_rest = np.concatenate(([0], chromosome.x[index:], [chromosome.destination]), axis=0)
	# print(paths)
	if verbose:
		result = np.sum([wbg.length(path) for path in paths])
	else:
		real_fitness = np.sum([wbg.length(path) for path in paths])
		rest_fitness = wbg.length(path_rest)
		result = real_fitness + rest_fitness/(real_fitness+rest_fitness)*real_fitness*REST_IMPORTANT
	return result


def crossover(population):
	np.random.shuffle(population)
	children = []
	loop = len(population)//2
	for i in range(loop):
		# child1, child2 = population[2*i].cross_over_ox(population[2*i+1])
		child1, child2 = population[2*i].cross_over_hai(population[2*i+1])
		if child1 is not None:
			population.append(child1)
		if child2 is not None:
			population.append(child2)
	# print(population)
	# print(children)
	# population = [population+ children]
	return population


def mutation(population, mutate_rate):
	pop_size = len(population)
	# rand1, rand2 = np.random.uniform(0, 1, 2)
	rand1 = np.random.uniform(0, 1)
	for i in range(pop_size):
		if rand1 < mutate_rate:
			x = copy.deepcopy(population[i])
			x.local_mutation()
			population.append(x)
		# if rand2 < mutate_rate:
		# 	x = copy.deepcopy(population[i])
		# 	result = x.global_mutation()
		# 	if result is not None:
		# 		population.append(result)
	return population


def selection(population, K, pop_size, verbose=False):
	prob_list = [compute_fitness(chromosome, K) for chromosome in population]
	if verbose:
		print(prob_list)
		arg_min = np.argmin(prob_list)
		fitness_list = [compute_fitness(chromosome, K, verbose=verbose) for chromosome in population]
		print("best fitness: ", min(fitness_list))
		print(population[arg_min].x)
	prob_list = normalize_prob(prob_list)
	population = np.random.choice(population, size=pop_size, replace=True, p=prob_list).tolist()
	population = sorted(population, key=lambda x: compute_fitness(x, K))
	return population
ITER = 10
REST_IMPORTANT = 0.5
EPSILON = 1e-5
POP_SIZE = 100


def normalize_prob(prob_list):
	prob_list = np.array(prob_list) ** 2 + EPSILON
	prob_list = 1 / (np.array(prob_list))
	# print("sum: ", np.sum(prob_list))
	pivot = np.mean(prob_list, axis=0)*0.9
	for i in range(0, len(prob_list)):
		if prob_list[i] < pivot:
			prob_list[i] = 0
	prob_list = prob_list / np.sum(prob_list)
	return prob_list


def main():
	from src.ultis import distance
	global wbg
	num_sensors = 30
	num_barrier = 3
	wbg = WBG(lenght=50, height=15, mode='strong')
	wbg.create_sensors_randomly(num_sensor=num_sensors, r=2, alpha=np.pi * 4 / 5)
	wbg.build_WBG()
	wbg.o_adj_matrix[0][num_sensors+1] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.o_adj_matrix[num_sensors+1][0] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.show_matrix()

	population = initialize(num_sensors=num_sensors, num_barrier=num_barrier, num_individuals=POP_SIZE)
	for iter in range(ITER):
		if iter % 10 == 0:
			print("Iter: ", iter)
		population = crossover(population)
		print("crossover", len(population))
		population = mutation(population, 0.3)
		print("mutation", len(population))
		population = selection(population, K=num_barrier, pop_size=POP_SIZE, verbose=iter % 10 == 0)
		print("seletion", len(population))
	Pk, Nm = wbg.min_num_mobile_greedy(num_barrier)
	print("######################################################")
	print(Pk)
	print(Nm)
	print("######################################################")
	wbg.field_show()

	pass


if __name__ == '__main__':
	main()
