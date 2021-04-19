import numpy as np
import networkx as nx
import matplotlib.pyplot as pyplot
from scipy import sparse as sp
import collections

def createGraph(file):
	G = nx.Graph()
	f = open(file, "r")
	for line in f:
		line = line.replace('\n', '')
		line = line.split()
		G.add_edge(line[0], line[1])

	# get neighbours
	G_neigh = collections.OrderedDict()
	for i in G.nodes():
		G_neigh[i] = list(G.neighbors(i))
	return G, G_neigh

def greedy(RList, pairs, min_eigenvector, min_len):
	subgraph = nx.Graph()
	num = 0
	for i in range(len(RList)):
		if (min_len < 0): break
		max_idx = RList.index(max(RList))
		max_Rval = RList[max_idx]
		RList[max_idx] = min_eigenvector - 1
		edge = pairs[max_idx]

		if (not subgraph.has_node(edge[0]) and not subgraph.has_node(edge[1])):
			print(edge[0], edge[1], max_Rval)
			subgraph.add_edge(*edge)
			num+=1
			min_len -=1
	print("Number of nodes in common: ", num)
	
def align(fileM, fileN, threshold):
	m, m_neigh = createGraph(fileM)
	n, n_neigh = createGraph(fileN)
	size = len(m.nodes()) * len(n.nodes())

	#use sparse matrix to solve eigenvalue problem efficiently
	A = sp.dok_matrix((size, size), dtype = float)

	#do a(i,j)(u,v) = ... (equation 4 in assignment)
	# print(len(m.nodes()))
	for i_idx, i in enumerate(m.nodes()):
		for j_idx, j in enumerate(n.nodes()):
			for u_idx, u in enumerate(list(m_neigh.values())[i_idx]):
				for v_idx, v in enumerate(list(n_neigh.values())[j_idx]):
					N_u = len(list(m_neigh.values())[u_idx])
					N_v = len(list(n_neigh.values())[v_idx])
					x_idx = i_idx*len(n.nodes())+j_idx
					y_idx = u_idx*len(n.nodes())+v_idx
					A[x_idx, y_idx] = (1/(N_u*N_v))

	#linalg.eigs returns vals, vecs -> R = vectors
	R_lambda, R = sp.linalg.eigs(A, k=1)

	# list of eigenvectors
	RList = [abs(r.real) for r in R]
	RList = np.asarray(RList, dtype=float)
	RList = np.extract(RList > threshold, RList)
	RList = list(RList)

	# get pairs
	pairs = []
	m_nodes_sorted = np.sort(m.nodes())
	n_nodes_sorted = np.sort(n.nodes())
	for i in range(len(m_nodes_sorted)):
		for j in range(len(n_nodes_sorted)):
			pairs.append((m_nodes_sorted[i], n_nodes_sorted[j]))
	
	min_eigvector = min(RList)
	min_len = min(len(m.nodes()), len(n.nodes()))

	greedy(RList, pairs, min_eigvector, min_len)

def main():
	print("N1 and N2")
	align("N1.txt", "N2.txt", 0.0005)
	print("N1 and N3")
	align("N1.txt", "N3.txt", 0.0005)
	print("N2 and N3")
	align("N2.txt", "N3.txt", 0.0005)

if __name__ == '__main__':
	main()