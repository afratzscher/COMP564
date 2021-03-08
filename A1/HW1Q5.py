import os
import subprocess
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

T1 = '((..(((.......))))).....(((......)))....'
T2 = '(((((..(((((.(((((....))))))))))..))))).'
T3 = '((((((.((((....))))..((((....)))).))))))'
N = 100 # size of population
gens = 10 # number of generations

test = '((((.((((....)))).((((....)))).))))'

#CORRECT
def generateNSeq(targetLen):
	N_seqs = []
	i = 0
	while (i<N):
		# equal probability to select A,C,G,U, so can use random
		N_seqs.append(''.join(np.random.choice(list('ACUG')) for _ in range(targetLen)))
		i += 1

	return N_seqs

#CORRECT
def getStructures(N_seqs):
	structs = []
	for seq in N_seqs: # run RNAfold to get MFE for all sequences
		p1 = subprocess.Popen(["printf", seq], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["RNAfold"], stdin=p1.stdout, stdout=subprocess.PIPE)
		to_parse = str(p2.communicate()[0])
		p1.stdout.close()
		p2.stdout.close()
		structs.append(to_parse[4+len(N_seqs[0]):to_parse.find(' ')])
		del to_parse
			#first 2 chars = b', then have seq, then have \\n before start
			# use to_parse.find(' ') b/c formatted as (...)(..)_(value)  (space before '(value)', so can find ' ''
	os.remove('rna.ps')
	del N_seqs
	return structs

#CORRECT
def getDistance(structs, target):
	distances = []
	for struct in structs:
		p1 = subprocess.Popen(["printf", struct+'\\n'+target], stdout=subprocess.PIPE)
		p2 = subprocess.Popen(["RNAdistance", "--distance=P"], stdin=p1.stdout, stdout=subprocess.PIPE)
		p1.stdout.close()
		string_to_parse = str(p2.communicate()[0])
		p2.stdout.close()
		distances.append(int(string_to_parse[string_to_parse.find(':')+2:string_to_parse.rfind(' ')-1]))
		del string_to_parse
	del structs
	return distances

#CORRECT
def getReproductionRate(distances, L):
	beta = 2/L
	Z = sum([math.exp(-1 * beta * d) for d in distances])
	rate = [math.exp(-1*beta*d)/Z for d in distances]
		# rate high if distance low
	return rate 

def getNextPop(N_seqs, rate, error_rate):
	N_seqs = [x for _,x in sorted(zip(rate, N_seqs), reverse=True)]
	rate = sorted(rate, reverse=True)
	rounded = [round(x*N) for x in rate]
	temp = [x for x, number in zip(N_seqs, rounded) for _ in range(number)]
	N_seqs = []
	for seq in temp:
		while(len(N_seqs) < N):
			new = ''
			for i in range(len(seq)):
				if (np.random.random() < error_rate): # mutation occurs
					alphabet = ['A', 'C', 'G', 'U']
					alphabet.remove(seq[i])
					new += np.random.choice(alphabet)
				else:
					new += seq[i]
			N_seqs.append(new)
	del temp
	del rounded
	return N_seqs

def getHammingDistance(structs, target):
	distances = []
	for struct in structs:
		distances.append(sum(c1 != c2 for c1, c2 in zip(struct, target)))
	avg = sum(distances)/N
	del distances
	del structs
	return avg

def main():
	target = test #declare target (T1, T2, T3)

	dist_low_mut = []
	dist_mid_mut = []
	dist_high_mut = []
	hamm_low_mut = []
	hamm_mid_mut = []
	hamm_high_mut = []

	for error_rate in [0.01, 0.02, 0.05]:
		print('starting for ' + str(error_rate))
		gen = 0 #reset generation
		N_seqs = generateNSeq(len(target)) # first generation
		while (gen <= gens):
			# print(gen)
			if(gen%20 == 0):
				print(gen)
			structs = getStructures(N_seqs)
			distances = getDistance(structs, target)
			hamming = getHammingDistance(structs, target)
			if error_rate == 0.01:
				dist_low_mut.append(sum(distances)/N)
				hamm_low_mut.append(hamming)
			elif error_rate == 0.02:
				dist_mid_mut.append(sum(distances)/N)
				hamm_mid_mut.append(hamming)
			elif error_rate == 0.05:
				dist_high_mut.append(sum(distances)/N)
				hamm_high_mut.append(hamming)
			del structs #clearing memory
			del hamming
			rate = getReproductionRate(distances, len(target))
			del distances
			N_seqs = getNextPop(N_seqs, rate, error_rate)
			del rate
			gen+=1
		
	plt.title('target = ' + str(target))
	plt.xlabel('Generations')
	plt.ylabel('Average distance')
	plt.plot(dist_low_mut, 'y.', label='u = 0.01')
	plt.plot(dist_mid_mut, 'r.', label='u = 0.02')
	plt.plot(dist_high_mut, 'b.', label='u = 0.05')
	plt.legend(loc="upper right")
	plt.savefig('old_bp.png')
	plt.show()

	plt.clf()
	plt.title('target = ' + str(target))
	plt.xlabel('Generations')
	plt.ylabel('Average hamming distance')
	plt.plot(hamm_low_mut, 'y.', label='u = 0.01')
	plt.plot(hamm_mid_mut, 'r.', label='u = 0.02')
	plt.plot(hamm_high_mut, 'b.', label='u = 0.05')
	plt.legend(loc="upper right")
	plt.savefig('old_Hamming_distance.png')
	# plt.show()

if __name__ == '__main__':
	main()
