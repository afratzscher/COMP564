import os
import subprocess
import sys
import math
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import random

T1 = '((..(((.......))))).....(((......)))....'
T2 = '(((((..(((((.(((((....))))))))))..))))).'
T3 = '((((((.((((....))))..((((....)))).))))))'
test = '((((.((((....)))).((((....)))).))))'
N = 100 # size of population # SHOULD BE 100
generations = 500 # number of generations SHOULD BE 500

# use class instead of list for easier management
class Individual:
	def __init__(self, seq, struct):
		self.seq = seq
		self.struct = struct

def getPop(targetLen):
	pop = []
	for i in range(N): # size of pop = N (100)
		seq = ''.join([random.choice("ACUG") for _ in range(targetLen)])
		p1 = subprocess.Popen(["printf", seq], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["RNAfold"], stdin=p1.stdout, stdout=subprocess.PIPE)
		to_parse = str(p2.communicate()[0])
		p1.stdout.close()
		p2.stdout.close()
		struct = to_parse[4+targetLen:to_parse.find(' ')]
		new_indi = Individual(seq, struct)
		new_indi.id = i
		pop.append(new_indi)
	os.remove('rna.ps')
	return pop

def getFitness(pop, target):
	beta = 2/len(target)
	dist_list = [] #total fitness
	for cell in pop:
		p1 = subprocess.Popen(["printf", cell.struct+'\\n'+target], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["RNAdistance", "--distance=P"], stdin=p1.stdout, stdout=subprocess.PIPE)
		p1.stdout.close()
		string_to_parse = str(p2.communicate()[0])
		p2.stdout.close()
		d = int(string_to_parse[string_to_parse.find(':')+2:string_to_parse.rfind(' ')-1])
		cell.fitness = math.exp((-1 * beta * d))
		cell.bp_dist = d
		dist_list.append(cell.fitness)
		cell.hamming_dist = sum(c1 != c2 for c1, c2 in zip(cell.struct, target)) #hamming distance
		
	Z = np.sum(dist_list)

	for cell in pop:
		cell.fitness = cell.fitness/ Z

	return pop

def mutate(seq, mutation_rate):
	new = ''
	mutated = False #check if mutated (to see if need to get struct again)
	for bp in seq:
		r = random.random()
		if r < mutation_rate:
			new = new + random.choice("ACUG")
			mutated=True
		else:
			new = new + bp # no change
	return(new, mutated)

def selection(pop, target, mutation_rate):
	parents = np.random.choice(pop, len(pop), p=[rna.fitness for rna in pop], replace=True)

	next_gen = []
	for i,p in enumerate(parents):
		new = Individual(p.seq, p.struct)
		new.id = i
		next_gen.append(new)

	#add mutations
	for rna in next_gen:
		mut_seq, mutated = mutate(rna.seq, mutation_rate)
		if mutated:
			rna.seq = mut_seq
			p1 = subprocess.Popen(["printf", mut_seq], stdout=subprocess.PIPE)
			p2 = subprocess.Popen(["RNAfold"], stdin=p1.stdout, stdout=subprocess.PIPE)
			to_parse = str(p2.communicate()[0])
			p1.stdout.close()
			p2.stdout.close()
			rna.struct = to_parse[4+len(target):to_parse.find(' ')]
		else:
			continue

	#update fitness
	getFitness(next_gen, target)
	return next_gen

def recordStats(pop, pop_stats):
	gen_bp_dist = [rna.bp_dist for rna in pop]
	gen_hamming_dist = [rna.hamming_dist for rna in pop]
	mean_bp_dist = np.mean(gen_bp_dist)
	mean_hamming_dist = np.mean(gen_hamming_dist)

	pop_stats.setdefault('mean_bp_dist', []).append(mean_bp_dist)
	pop_stats.setdefault('mean_hamming_dist', []).append(mean_hamming_dist)
	return None

def evolve(target, mutation_rate):
	pop_stats = {} #holds stats

	#for initial population
	initial_pop = getPop(len(target))
	getFitness(initial_pop, target)
	curr_gen = initial_pop

	#run for 500 generations
	for i in range(generations):
		if(i%20 == 0):
				print(i)
		#record stats
		recordStats(curr_gen, pop_stats)
		#get next generation
		new_gen = selection(curr_gen, target, mutation_rate=mutation_rate)
		curr_gen = new_gen
	return pop_stats # return stats

def main():
	target = test
	mutation_rate = [0.01, 0.02, 0.05]
	all_stats = {}
	for m in mutation_rate:
		stats = evolve(target, m)
		all_stats[m] = stats

	plt.title('target = ' + str(target))
	plt.xlabel('Generations')
	plt.ylabel('Average base pair distance')

	plt.plot(all_stats[0.01]['mean_bp_dist'], 'y.', label='u = 0.01')
	plt.plot(all_stats[0.02]['mean_bp_dist'], 'r.', label='u = 0.02')
	plt.plot(all_stats[0.05]['mean_bp_dist'], 'b.', label='u = 0.05')
	plt.legend(loc="upper right")
	plt.savefig('BP_distance.png')

	plt.clf()
	plt.title('target = ' + str(target))
	plt.xlabel('Generations')
	plt.ylabel('Average hamming distance')
	plt.plot(all_stats[0.01]['mean_hamming_dist'], 'y.', label='u = 0.01')
	plt.plot(all_stats[0.02]['mean_hamming_dist'], 'r.', label='u = 0.02')
	plt.plot(all_stats[0.05]['mean_hamming_dist'], 'b.', label='u = 0.05')
	plt.legend(loc="upper right")
	plt.savefig('Hamming_distance.png')

if __name__ == '__main__':
	main()
