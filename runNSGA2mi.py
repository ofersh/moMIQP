# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:20:45 2022
@author: ofers
"""

# import numpy as np
# from pymoo.core.problem import ElementwiseProblem
# from pymoo.core.variable import Real, Integer #, Choice, Binary
from pymoo.visualization.scatter import Scatter
from pymoo.algorithms.moo.nsga2 import NSGA2 #, RankAndCrowdingSurvival
from pymoo.core.mixed import MixedVariableMating, MixedVariableSampling, MixedVariableDuplicateElimination #, MixedVariableGA
from pymoo.optimize import minimize
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
# from pymoo.visualization.pcp import PCP
import ellipsoidFunctions as Efunc
import numpy as np
NRUNS = 1
problemName = 'MIQP-HadEllipse'
N = Efunc.N;
c = 100.0
H = Efunc.genHadamardHellipse(2*N, c) #genRotatedHellipse(2*N, c) # genHellipse(N, c)
problem = Efunc.MixedVarsRotEllipsoid(H,c)
popSize = 15
MAX_ITER = 1000
algorithm = NSGA2(pop_size=popSize,
                  sampling=MixedVariableSampling(),
                  mating=MixedVariableMating(eliminate_duplicates=MixedVariableDuplicateElimination()),
                  eliminate_duplicates=MixedVariableDuplicateElimination(),
                  )
for k in range(NRUNS) :
    res = minimize(problem, algorithm, ('n_gen', MAX_ITER), seed=10+k, verbose=False)
    filename=problemName+str(2*N)+'_'+str(int(c))+'_'+str(k+1)+'_PF.dat'
    with open(filename, 'w') as fout :
        for i in range (popSize) :
           fout.write(str(res.F[i][0]) + ' ' + str(res.F[i][1]) + '\n')
    filename2=problemName+str(2*N)+'_'+str(int(c))+'_'+str(k+1)+'_PS.dat'
    with open(filename2, 'w') as fout2 :
        for j in range (popSize) :
            y = np.array([res.X[j][f"x{k:02}"] for k in range(1, 2*N+1)])
            fout2.write(str(y) + '\n')
    print(k)
#Scatter().add(res.F).show()
#
# algorithm = NSGA2(
#     pop_size=40,
#     n_offsprings=10,
#     sampling=FloatRandomSampling(),
#     crossover=SBX(prob=0.9, eta=15),
#     mutation=PM(eta=20),
#     eliminate_duplicates=True
# )