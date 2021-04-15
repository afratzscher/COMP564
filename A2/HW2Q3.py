import numpy as numpy
import matplotlib.pyplot as plt
import os

x, y = [], []

filename = os.getcwd() + '/' + 'MD/workspace/molecule-bkbone-rmsd.xvg'
with open(filename) as file:
	for line in file:
		if not ('#' in line):
			if not ('@' in line):
				cols = line.split()
				x.append(float(cols[0]))
				y.append(float(cols[1]))

plt.title("RMSD vs time")
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (nm)")
plt.plot(x,y)
plt.savefig("RMSD.pdf")