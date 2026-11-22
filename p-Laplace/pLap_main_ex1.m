%% p-Laplacian - example 1

close all
clear variables
%% Settings
coff_p = 1.5;
pde = pLapdata1(coff_p);

% mg options
option.mg.solvermaxIt = 5;
option.mg.solver = 'Vcycle';
option.mg.tol = 1e-2;
option.mg.printlevel = 0;
option.rateflag = 1;
option = mfemoption(option);


%% Initialize err and time array
maxIt = 3;
maxN = 1e7;


% error
errsigmaL2 = zeros(maxIt,1);   
erruIuhL2 = zeros(maxIt,1); 
erruL2 = zeros(maxIt,1);   
uL2 = zeros(maxIt,1);   
relerruL2 = zeros(maxIt,1);   
errsigmaIsigmah = zeros(maxIt,1); 

% info
N = zeros(maxIt,1); 
h_lst = zeros(maxIt,1);
itStep = zeros(maxIt,1);
time = zeros(maxIt,1);
flag = zeros(maxIt,1);
Vcycle_num = zeros(maxIt,1);

%% Solvers
solverNames = {'DualTPDJ', 'DualTPDM'}; 
for i = 1:length(solverNames)
    % Generate mesh 
    h = 1/4;  % coarsest mesh 
    [node,elem] = squaremesh([0,1,0,1],h); 
    mesh = struct('node',node,'elem',elem);
    bdFlag = setboundary(node,elem,'Dirichlet');

   
    L0 = 4;
    for k = 1:L0
        [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        [~,~,isBdNode] = findboundary(elem,bdFlag);
        h = h/2;
    end

    solver = solverNames{i};
    fprintf('\nSolver: %s\n', solver);
    for k = 1:maxIt
                
        disp('solver begins')
        % solve the equation
        switch solver 

            case 'DualTPDJ'
                option.stepsize = 1;
                option.eps0 = 1e-4;
                [sigma,u,eqn,info] =  pLapP0P1_DualTPDJ(node,elem,bdFlag,pde,option);


            case 'DualTPDM'
                option.stepsize = 0.6;
                [sigma,u,eqn,info] =  pLapP0P1_DualTPDM(node,elem,bdFlag,pde,option);

            otherwise
                disp('Solver not specified.')
        end

        % compute error
        if isfield(pde,'exactsigma')
            errsigmaL2(k) = getL2error(node,elem,pde.exactsigma,sigma);
            sigmaI = Lagrangeinterpolate(pde.exactsigma,node,elem,'P0');
            sigmaIsigmah_diff = sqrt(sum((sigmaI-sigma).^2, 2));
            area = simplexvolume(node,elem);
            errsigmaIsigmah(k)=sqrt(dot((sigmaIsigmah_diff).^2,area));
        end
        if isfield(pde,'exactu')
            erruL2(k) = getL2error(node,elem,pde.exactu,u);
            uL2(k) = getL2error(node,elem,pde.exactu,zeros(size(u)));
            relerruL2(k) = erruL2(k)/uL2(k);
            % interpolation
            uI = Lagrangeinterpolate(pde.exactu,node,elem,'P1');        
            uIuh_diff = (uI(elem(:, 1)) - u(elem(:, 1)) + uI(elem(:, 2)) - u(elem(:, 2))...
                + uI(elem(:, 3)) - u(elem(:, 3)))/3;
            erruIuhL2(k) = sqrt(dot((uIuh_diff).^2,area));
        end
        
        %  solver information
        itStep(k) = info.itStep;
        flag(k) = info.flag;
        Vcycle_num(k) = info. Vcycle_num;
        time(k) = info.time;

        fprintf('\n ite stp: = %d \n',info.itStep);
        fprintf('\n time: = %g \n',info.time);
    
        % plot 
        N(k) = length(sigma) + length(u);
        h_lst(k) = h;
    
        if N(k) > maxN
            break;
        end
       
        % refine mesh
        [node,elem, bdFlag] = uniformrefine(node,elem,bdFlag);
        [~,~,isBdNode] = findboundary(elem,bdFlag);
        h = h/2;
    end

    % print iteration results 
    fprintf('\nIteration Results: %s\n', solver);
    fprintf('---------------------------------------------------------\n');
    iterationResults = table(h_lst, itStep, Vcycle_num, time, ...
        'VariableNames', {'h', 'Iteration', 'Vcycle', 'Time'});
    disp(iterationResults);


 end
