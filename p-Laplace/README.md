# p-Laplacian experiments code

This folder implements the p-Laplacian experiments in **A Novel Preconditioning Framework for Solving Nonlinear PDEs Based on Fenchel-Rockafellar Duality and Transformed Primal-Dual Techniques**, by Long Chen, Ruchi Guo, Jingrong Wei, and Jun Zou. These scripts implement the numerical experiments in Section 4.

### Setup and first run

Install [iFEM package](https://github.com/lyc102/ifem) and add it to the MATLAB path. Set the current folder to `p-Laplace`, then run `pLap_main_ex1` or `pLap_main_ex2`.

Set `coff_p` to choose p and edit `solverNames`  to select the solvers. For example:

```matlab
solverNames = {'DualTPDJ', 'DualTPDM'};
```


### Solvers


| Solver name | File | Method |
| --- | --- | --- |
| `DualTPDJ` | `pLapP0P1_DualTPDJ.m` | DualTPD-J with Jacobian preconditioning |
| `DualTPDM` | `pLapP0P1_DualTPDM.m` | DualTPD-M with mass preconditioning |
| `PGD` | `pLapP0P1_primal.m` (p > 2); `pLapP0P1_primalsmallp.m` (p <= 2) | Preconditioned gradient descent with line search |
| `PGDw/oLS` | `pLapP0P1_primal_nolinesearch.m` (p > 2); `pLapP0P1_primalsmallp_nolinesearch.m` (p <= 2) | Preconditioned gradient descent without line search |
| `Newton` | `pLapP0P1_Newton.m` | Newton method for the primal formulation |
| `PGD+Newton` | `pLapP0P1_PGDNewton.m` | PGD followed by Newton once the relative residual falls below 1e-2 |

### Experiment-to-file mapping

#### First experiment: on the unit square

| File | Role |
| --- | --- |
| `pLap_main_ex1.m` | Main file  |
| `pLapdata1.m` | Exact solution, source terms, and nonlinear coefficients |

The paper's solver-performance study uses both p = 1.5 and p = 4, with these step sizes:

| coff_p | DualTPD-J | DualTPD-M |
| --- | --- | --- |
| 1.5 | 1 | 0.6 |
| 4 | 0.6 | 1.2 |


#### Second experiment: constant source on the unit disk


| File | Role |
| --- | --- |
| `pLap_main_ex2.m` | Main file |
| `pLapdata2.m` | Exact solution, source terms, and nonlinear coefficients |



