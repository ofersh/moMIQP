# moMIQP
Comparing White-Box and Black-Box Solvers for Multi-Objective Mixed-Integer Quadratic Optimization Problems

1. The quadratic models are defined in python by means of ellipsoidFunctions.py

2. The NSGA-II implementation relies on the pymoo package. To invoke it on the current MIQP setup, use runNSGA2mi.py 

3. The CPLEX-based DMA algorithm is implemented in OPL, using the following files - 
    modelR0.mod, modelR1.mod - underlying quadratic models
    dma_miRotEllipse_1.mod - the executation of the solver over MaxNumOfPoints iterations
    miRotEllipse.dat - data file with parameters + the defining Hessian matrix
    
4. The SMS-EMOA implementation is provided in MATLAB source files - *.m 

5. The obtained experimental datasets are located within the /datasets/ folder.

6. The attainment curves' calculation was done using the R script EAFplots.R
