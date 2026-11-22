%% First experiment: nonlinear ferromagnetism
% Add iFEM to the MATLAB path and run from the Nonlinear_Maxwell folder.

% Dual TPD
fprintf('\nSolver: DualTPD\n');
Dual_cubeMaxwellExp_DualTPD;


% Dual PD
fprintf('\nSolver: DualPD\n');
Dual_cubeMaxwellExp_DualPD;


% Primal TPD
fprintf('\nSolver: PrimalTPD\n');
cubeMaxwellExp_PrimalTPD;


% Newton-Krylov
fprintf('\nSolver: Newton-Krylov\n');
cubeMaxwellExp_Primal_NewtonKrylov;
