%% Second experiment: p-curl problem
% Add iFEM to the MATLAB path and run from the Nonlinear_Maxwell folder.

% Dual TPD, p = 1.5
fprintf('\nSolver: DualTPD, p = 1.5\n');
Dual_cubeMaxwellpcurl_DualTPD_ex1;


% Dual PD, p = 1.5
fprintf('\nSolver: DualPD, p = 1.5\n');
Dual_cubeMaxwellpcurl_DualPD_ex1;


% Dual TPD, p = 3
fprintf('\nSolver: DualTPD, p = 3\n');
Dual_cubeMaxwellpcurl_DualTPD_ex2;


% Dual PD, p = 3
fprintf('\nSolver: DualPD, p = 3\n');
Dual_cubeMaxwellpcurl_DualPD_ex2;
