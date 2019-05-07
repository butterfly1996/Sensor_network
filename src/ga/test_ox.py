# 02373822922
import copy

import numpy as np

from src.ga.individuals_ox import Chromosome
from src.sensors.sensors_field import WBG

np.random.seed(1234)


def initialize(num_sensors, num_barrier, num_individuals):
	return [Chromosome(num_sensors=num_sensors, K=num_barrier).optimize_path() for i in range(0, num_individuals)]


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
		child1, child2 = population[2 * i].cross_over_ox(population[2 * i + 1], pc=0.5)
		if child1 is not None:
			child1.optimize_path()
			population.append(child1)
		if child2 is not None:
			child2.optimize_path()
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
	prob_list = prob_list / np.sum(prob_list)
	return prob_list


def main():
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

	# print(wbg.o_adj_matrix[0][24], wbg.o_adj_matrix[24][12], wbg.o_adj_matrix[12][21], wbg.o_adj_matrix[21][17],
	# 		wbg.o_adj_matrix[17][6], wbg.o_adj_matrix[6][4], wbg.o_adj_matrix[4][14],wbg.o_adj_matrix[14][31])
	wbg.field_show()
	pass


if __name__ == '__main__':
	main()
