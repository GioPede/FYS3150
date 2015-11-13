import commands,os
import matplotlib.pyplot as plt
import numpy as np

def getext(a):
	"""

	function analogous to the one in c++
	a=name of the program
	
	"""
	out = open(a,"r")
	
	c = out.read()
	b=c.split()
	i=0
	array= dict()
	j=0
	while (b[i]!="break"):
		array[j] = []
		while (b[i]!='change'):
			array[j].append(float(b[i]))
			i+=1
		i+=1
		j+=1
	for x in range(j):
		yield array[x]

a,b,c = getext("Data/test_t24_rand.txt")
for i in range(0,len(b)):
	b[i]/=400

plt.figure()
plt.plot(a, b, '-b')
plt.grid()
plt.xscale('log');
#plt.axis([0, 1000000, -2, -1.990])
plt.xlabel("Cycles")
plt.ylabel("E")
plt.savefig('Plots/ENERGY_T24_RAND_10000MEAN.png')

plt.figure()
plt.plot(a, c,'-b')
plt.xscale('log');
plt.axis([0, 1000000, 0, 0.03])
plt.grid()
plt.xlabel("Cycles")
plt.ylabel("M")
plt.savefig('Plots/MAGNETIZATION_T2_RAND_10000MEAN.png')

#plt.figure()
#plt.plot(a, d, '.b')
#plt.xscale('log');
#plt.grid()
#plt.xlabel("cycles")
#plt.ylabel("Susceptibility")
#plt.title("Susceptbility")
#plt.savefig('Susceptibility.png')
