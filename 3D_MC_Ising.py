# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 15:31:27 2021

@author: devoted4god

Slightly modified version of code by rajeshrinet from https://rajeshrinet.github.io/blog/2014/ising-model/#Monte-Carlo-simulation-of-2D-Ising-model
to a 3D Ising model
"""
import numpy as np

N = 10

""" input: size - width, length, depth = N
    output: 3D lattice with randomly distributed spins
"""
def initialState(N):
    state = 2*np.random.randint(2, size=(N,N,N))-1  # 3D lattice size used everywhere else
    return state
    
""" input: initial state for evolution, inverse temperature
    output: congiguration after one MC step
"""

def mcStep(config, beta):

    for i in range(N):
        for j in range(N):
            for k in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)
                c = np.random.randint(0, N)
                S = config[a,b,c]
                nb = config[(a+1)%N,b,c] + config[(a-1)%N,b,c] + config[a,(b+1)%N,c] + \
                    + config[a,(b-1)%N,c] + config[a,b,(c+1)%N] + config[a,b,(c-1)%N]
                cost = 2*S*nb
                if cost < 0:
                    S *= -1
                elif np.random.uniform(0,1) < np.exp(-cost*beta):
                    S *= -1
                config[a,b,c] = S
    return config

""" input: configuration of the 3D spin lattice
    output: energy value of the given configuration
    """

def calcEnergy(config):
    energy = 0
    for i in range(len(config)):
        for j in range(len(config)):
            for k in range(len(config)):
                S = config[i,j,k]
                nb = config[(i+1)%N,j,k] + config[(i-1)%N,j,k] + config[i,(j+1)%N,k] + \
                    + config[i,(j-1)%N,k] + config[i,j,(k+1)%N] + config[i,j,(k-1)%N]
                energy += -nb*S
    return energy/6.

"""
input: configuration of the 3D spin lattice
output: magnetization of the given configuration 
"""

def calcMag(config):
    mag = np.sum(config)
    return mag




## change these parameters for a smaller (faster) simulation 
nt      = 40         #  number of temperature points
N       = 10         #  size of the lattice, N x N
eqSteps = 1024       #  number of MC sweeps for equilibration
mcSteps = 1024       #  number of MC sweeps for calculation

T       = np.linspace(3.0, 7.0, nt); 
E,M,C,X = np.zeros(nt), np.zeros(nt), np.zeros(nt), np.zeros(nt)
n1, n2  = 1.0/(mcSteps*N*N*N), 1.0/(mcSteps*mcSteps*N*N*N) 
# divide by number of samples, and by system size to get intensive values
                
                
for tt in range(nt):
    E1 = M1 = E2 = M2 = 0
    config = initialState(N)
    iT = 1.0/T[tt]; iT2=iT*iT; #inverse temperature
    
    for i in range(eqSteps):  #equilibrate
        mcStep(config,iT)
    
    for i in range(mcSteps):
        mcStep(config, iT)
        Ene = calcEnergy(config)
        Mag = calcMag(config)
        
        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag 
        E2 = E2 + Ene*Ene
        
    E[tt] = n1*E1
    M[tt] = n1*M1
    C[tt] = (n1*E2 - n2*E1*E1)*iT2
    X[tt] = (n1*M2 - n2*M1*M1)*iT

import matplotlib.pyplot as plt    
                
f = plt.figure(figsize=(18, 10)); # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 );
plt.scatter(T, E, s=50, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);         plt.axis('tight');

sp =  f.add_subplot(2, 2, 2 );
plt.scatter(T, abs(M), s=50, marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight');

sp =  f.add_subplot(2, 2, 3 );
plt.scatter(T, C, s=50, marker='o', color='IndianRed')
plt.xlabel("Temperature (T)", fontsize=20);  
plt.ylabel("Specific Heat ", fontsize=20);   plt.axis('tight');   

sp =  f.add_subplot(2, 2, 4 );
plt.scatter(T, X, s=50, marker='o', color='RoyalBlue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Susceptibility", fontsize=20);   plt.axis('tight');                
                
                
                
                
                
                
                
                
                
                
                
    
                
                    
                
                