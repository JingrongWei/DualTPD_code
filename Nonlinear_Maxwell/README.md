# Nonlinear Maxwell experiments code

This folder implements the nonlinear Maxwell experiments in **A Novel Preconditioning Framework for Solving Nonlinear PDEs Based on Fenchel-Rockafellar Duality and Transformed Primal-Dual Techniques**, by Long Chen, Ruchi Guo, Jingrong Wei, and Jun Zou. These scripts implement the numerical experiments in Section 5.4.

### Setup and first run

Install [iFEM package](https://github.com/lyc102/ifem) and add it to the MATLAB path. 


Main scripts for the nonlinear ferromagnetism and p-curl problems.

```matlab
Maxwell_main_ex1   % First experiment: nonlinear ferromagnetism
Maxwell_main_ex2   % Second experiment: p-curl problem
```

Each main script runs the following corresponding solver scripts sequentially. Mesh levels and solver parameters are set in the individual solver scripts.

### Experiment-to-file mapping

#### First experiment: nonlinear ferromagnetism

Run [Maxwell_main_ex1.m](Maxwell_main_ex1.m) to execute DualTPD, DualPD, PrimalTPD, and Newton-Krylov in the order listed below.

| File | Solver |
| --- | --- |
| `Dual_cubeMaxwellExp_DualTPD.m` | DualTPD |
| `Dual_cubeMaxwellExp_DualPD.m` | DualPD (Preconditioned Chambolle-Pock variant) |
| `cubeMaxwellExp_PrimalTPD.m` | TPD |
| `cubeMaxwellExp_Primal_NewtonKrylov.m` | Damped Newton-Krylov |

Use `maxIt = 5` in `Dual_cubeMaxwellExp_DualTPD.m` for all five rows of the paper's result

#### Second experiment: p-curl problem

Run [Maxwell_main_ex2.m](Maxwell_main_ex2.m) to execute DualTPD and DualPD for p = 1.5, followed by DualTPD and DualPD for p = 3, in the order listed below.

| File | p | Solver |
| --- | --- | --- |
| `Dual_cubeMaxwellpcurl_DualTPD_ex1.m` | p = 1.5 | DualTPD |
| `Dual_cubeMaxwellpcurl_DualPD_ex1.m` | p =  1.5 | DualPD |
| `Dual_cubeMaxwellpcurl_DualTPD_ex2.m` | p = 3 | DualTPD |
| `Dual_cubeMaxwellpcurl_DualPD_ex2.m` | p = 3 | DualPD |

Use `maxIt = 4` for all four rows of the paper's result. Reduce `maxIt` to run small size examples only. 

### Supporting files
- `MaxwellsaddleMat.m`: assembly of the Maxwell saddle-point matrices and boundary/source data.
- `mgMaxwell_1.m`: tensor-coefficient Maxwell solver with Hiptmair-Xu auxiliary-space preconditioning and multigrid.
- `mgHodgeLapELarge.m`: Hodge-Laplacian multigrid solver used by PrimalTPD.
- `data1.m`, `data2.m`, `data3.m`: nonlinear coefficients, derivatives, exact solutions, and source terms data for the experiments.
- `numerical_inverse.m`: numerical constitutive inversion used in the dual formulation.
