# -*- coding: utf-8 -*-
"""
Created on Sun Nov  5 14:04:37 2023

@author: Gilles Boonen
"""

" IMPORTS "
from IPython import get_ipython;
import numpy as np
import matplotlib.pyplot as plt
from ypstruct import structure
import fct
import codage as cd
# get_ipython().magic('reset -sf')

import ga_mono

" CODE "

#  Test functions
def f1_dejong(x):
    return sum(x**2)

def f2_dejong(x):
    n = len(x)
    result = 0
    for i in range(n - 1):
        result += 100 * (x[i + 1] - x[i]**2)**2 + (x[i] - 1)**2
    return result

def rastrigin(x):
    A = 10
    return A * len(x) + sum([(xi**2 - A * np.cos(2 * np.pi * xi)) for xi in x])

def easom(x):
    return -np.cos(x[0]) * np.cos(x[1]) * np.exp(-(x[0] - np.pi)**2 - (x[1] - np.pi)**2)


(t, date, P_el_total) = fct.GetProfile('groupe_5.csv', Shouse=1.)  
LDC_elec = fct.LoadDurationCurve(t, date, P_el_total, color='blue', name='Electricity Demand', steps=250)


(t, date, tau_pv) = fct.GetPVProduction('SolarData.csv', Pnom = 3000) #Pnom = puissance crête de l'installation PV en kWc
LDC_PV = fct.LoadDurationCurve(t, date, tau_pv, color='green', name='PV production', steps=250)


# Problem definition
problem = structure()
problem.costfunc = fct.CalculFacture(P_el_total, tau_pv, 150, 50) # Choix fonction test
problem.nvar = 1
problem.varmin = 1000
problem.varmax = 5000

# GA Parameters
params = structure()
params.maxit = 30
params.npop = 100
params.pc = 1
params.gamma = 0.1
params.mu = 0.01
params.sigma = 0.1
params.beta = 1

plt.clf()

# Run 
it = 3 # Nombre de run total
sol = np.zeros(it)
for i in range(it):
    out = ga_mono.run(problem, params)
    sol[i] = out.bestsol.cost

    # Results
    #plt.clf()
    
    plt.plot(out.bestcost)          # Plot linéaire
    #plt.semilogy(out.bestcost)     # Plot logarithmique
    
    plt.xlim(0, params.maxit)
    plt.xlabel('Iterations')
    plt.ylabel('Best Cost')
    plt.title('Genetic Algorithm (GA)')
    plt.grid(True)

plt.show()
for j in range(it):
    print("\nBest solution {} : ".format(j+1),sol[j])
