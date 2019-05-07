# 02373822922
import copy
from src.constants.my_constant import NUM_BARRIERS, NUM_SENSORS, POP_SIZE, ITER
import numpy as np

from src.ga.individuals_ox import Chromosome
from src.sensors.sensors_field import WBG

np.random.seed(1234)


def initialize(num_sensors, num_barrier, num_individuals, wbg):
	population = [Chromosome(num_sensors=num_sensors, K=num_barrier, wbg=wbg) for i in range(0, num_individuals)]
	for chromsome_ele in population:
		chromsome_ele.optimize_path()
	return population


def compute_fitness(chromosome, K, verbose=False):
	# print(chromosome.x)
	# print(K)
	paths = chromosome.convert_to_path()
	# print(paths)
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
	return population


def mutation(population, mutate_rate):
	# pop_size = len(population)
	for chromosome in population:
		rand1 = np.random.uniform(0, 1)
		if rand1 < mutate_rate:
			x_copy = copy.deepcopy(chromosome.x)
			mutate_chromosome = Chromosome(chromosome.num_sensors, chromosome.K, x=x_copy, wbg=chromosome.wbg)
			mutate_chromosome.local_mutation()
			if mutate_chromosome is not None:
				population.append(mutate_chromosome)
	return population


def selection(population, K, pop_size, verbose=False):
	prob_list = [compute_fitness(population[i], K) for i in range(len(population))]
	if verbose:
		# print(prob_list)
		arg_min = np.argmin(prob_list)
		print("best fitness: ", prob_list[arg_min])
		print(population[arg_min].x)
		print(population[arg_min].convert_to_path())
	rank_list = np.argsort(prob_list)
	new_population = [population[id] for id in rank_list[:pop_size-1]]
	new_population.append(Chromosome(num_sensors=population[rank_list[0]].num_sensors,
									K=population[rank_list[0]].K,
									x=copy.deepcopy(population[rank_list[0]].x).astype(int),
									wbg=population[rank_list[0]].wbg))
	# print("pop: ", len(new_population))
	return new_population



def main():

	global wbg
	wbg = WBG(lenght=500, height=100, mode='strong')
	wbg.create_sensors_randomly(num_sensor=NUM_SENSORS, r=20, alpha=np.pi * 1 / 4)
	wbg.build_WBG()
	wbg.o_adj_matrix[0][NUM_SENSORS+1] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.o_adj_matrix[NUM_SENSORS+1][0] = np.ceil(wbg.L/wbg.sensors_list[0].lr)
	wbg.show_matrix()
	population = initialize(num_sensors=NUM_SENSORS, num_barrier=NUM_BARRIERS, num_individuals=POP_SIZE, wbg=wbg)
	# khoi tao heuristic
	Pk, Nm = wbg.min_num_mobile_greedy(NUM_BARRIERS)
	print("######################################################")
	print(Pk)
	sen_set = set()
	for path in Pk:
		for sen in path:
			if sen > 0:
				sen_set.add(sen)
	# Pk.append([])
	for a in range(1, NUM_SENSORS + 1):
		if a not in sen_set:
			Pk[-1].append(a)
	x = np.concatenate(Pk)
	print(Nm)
	# print("######################################################")
	# population[-1].optimize_path()
	# population.append(Chromosome(num_sensors=num_sensors,
	# 							K=num_barrier,
	# 							wbg=wbg,
	# 							x=x))
	# population[-1].optimize_path()
	# print(population[-1].x)
	# x = compute_fitness(population[-1], num_barrier)
	print("x::: ^")
	for iter in range(ITER):
		if iter % 10 == 0:
			print("Iter: ", iter)
		population = crossover(population)
		# print("crossover", len(population))
		population = mutation(population, 0.1)
		# print("mutation", len(population))
		population = selection(population, K=NUM_BARRIERS, pop_size=POP_SIZE, verbose=iter % 10 == 0)
		# print("seletion", len(population))

	# x_test = [[0, 24, 12, 21, 17, 6, 4, 14, 31], [0, 22, 25, 16, 2, 30, 31], [0, 1, 26, 31]]
	# print(np.sum(wbg.length(path) for path in x_test))
	# print(wbg.o_adj_matrix[0][24], wbg.o_adj_matrix[24][12], wbg.o_adj_matrix[12][21], wbg.o_adj_matrix[21][17],
	# 		wbg.o_adj_matrix[17][6], wbg.o_adj_matrix[6][4], wbg.o_adj_matrix[4][14],wbg.o_adj_matrix[14][31])
	print("######################################################")
	print(Pk)
	print(Nm)
	wbg.field_show()
	pass


if __name__ == '__main__':
	main()
