""""
Monte Carlo simulation of the 2D Ising model with Umbrella Sampling
"""

import networkx as nx
from scipy import *
from scipy import weave
from pylab import *

Nitt = 100000  # total number of Monte Carlo steps
N = 20          # linear dimension of the lattice, lattice-size= N x N
warm = 50000     # Number of warmup steps
measure= 1     # How often to take a measurement


def RandomL(N):
    "Random lattice, corresponding to infinite temerature"
    latt = zeros((N,N), dtype=int)
    for i in range(N):
        for j in range(N):
            genrand = rand()
	    if genrand >= 0.5:
            	latt[i,j] = 1
            if genrand < 0.5:
		latt[i,j] = 0
    return latt

def CEnergy(latt):
    "Energy of a 2D Ising lattice at particular configuration"
    Energy = 0
    Chemical_potential_term = 0
    Interaction_energy = 0
    Particles = 0
    for i in range(len(latt)):
        for j in range(len(latt)):
            S = latt[i,j]
            Neighbor_interaction = -((latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1),j] + latt[i,(j-1)])*S)
	    Chemical_potential_term += S
	    Interaction_energy += Neighbor_interaction
            Energy = Chemical_potential_term + Interaction_energy
	    Particles += S
    return Energy/2

def SamplePython(Nitt, latt):
    "Monte Carlo sampling for the Ising model in Python"
    Energy = CEnergy(latt)  # Starting energy
    Particles  = sum(latt)         # Starting magnetization
    sample = open('%s_brute.dat' % (Target), 'w')
    N2 = float(N*N)
    Density = Particles/N2
    k = 0 #40000*Target**2 - 40000*Target + 10000
    W = 0.5*k*(Density-Target)**2

    Naver=0       # Measurements
    Eaver=0.0     # Additive energy
    Paver=0.0     # Additive particles
    Eaver2 = 0.0  # Energy^2
    Paver2 = 0.0  # Particles^2

    for itt in range(Nitt):
        t = int(rand()*N2)										
        (i,j) = (t % N, t/N)
        S = latt[i,j]         										  # Randomly selected lattice site
	Neighbor_sum = latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1),j] + latt[i,(j-1)]
        Neighbor_interaction = -((latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1),j] + latt[i,(j-1)])*S)  # Ineraction energy with site neighbors 
	Ham = 5*U*S + Neighbor_interaction + W 								  # Total Hamiltonian of Site with an Umbrella Weight	

	if S == int(1):     	 					# Particle is going
		Switch_dens = (Particles - 1)/N2 			# Density if particle is a removed		
		Umbrel_wt   = 0.5*k*(Switch_dens-Target)**2 		# Umbrella weight if particle is a removed
		E_diff = -5*U*S + -Neighbor_interaction + (Umbrel_wt - W)	# Energy difference of removing a particle

		if E_diff < 0: 
			latt[i,j] = int(0)					# Switch particle
			Energy += -5*U*S + -Neighbor_interaction + (Umbrel_wt - W)	# Change in total E
			Particles += -1						# Subtract 1 particle
			Density = Particles/N2					# Update density
			W = 0.5*k*(Density-Target)**2				# Update umbrella weight

		if E_diff >= 0:
			x = rand()
			if x < exp(Neighbor_interaction*T)*exp(U)*exp(W-Umbrel_wt):
				latt[i,j] = int(0)
				Energy += -5*U*S + -Neighbor_interaction + (Umbrel_wt - W)
				Particles += -1
				Density = Particles/N2
				W = 0.5*k*(Density-Target)**2

	if S == int(0):      # Particle is arriving
		Switch = 1
		Switch_dens = (Particles + 1)/N2 
		Umbrel_wt = 0.5*k*(Switch_dens-Target)**2
		Switch_interaction = -((latt[(i+1)%N, j] + latt[i,(j+1)%N] + latt[(i-1),j] + latt[i,(j-1)])*Switch)
		Switch_H = 5*U*Switch + Switch_interaction + Umbrel_wt
		E_diff = Switch_H - W

		if E_diff < 0:
			latt[i,j] = int(1)
			Energy += Switch_H - W
			Particles += 1
			Density = Particles/N2
			W = 0.5*k*(Density-Target)**2

		if E_diff >= 0:
			x = rand()
			Boltzmann = exp(-Switch_interaction*T)*exp(-U)*exp(Umbrel_wt-W)
			if x < exp(-Switch_interaction*T)*exp(-U)*exp(W-Umbrel_wt):
				latt[i,j] = int(1)
				Energy += Switch_H - W
				Particles += 1
				Density = Particles/N2
				W = 0.5*k*(Density-Target)**2

        if itt>warm and itt%measure==0:
	    sample.write('%s\n' % (Density)) 
            Naver += 1
            Eaver += Energy
            Paver += Particles
            Paver2 += Particles**2
            Eaver2 += Energy**2
    return (Paver/Naver, Eaver/Naver, Eaver2/Naver, Paver2/Naver)
if __name__ == '__main__':
    latt = RandomL(N)
    T = 0.2
    U = 0.4
    #(paver, eaver, eaver2, paver2) = SamplePython(Nitt, latt)
    #print paver/(N*N)
    wU = linspace(0,1,101)
    for Target in wU:
        (paver, eaver, eaver2, paver2) = SamplePython(Nitt, latt)
