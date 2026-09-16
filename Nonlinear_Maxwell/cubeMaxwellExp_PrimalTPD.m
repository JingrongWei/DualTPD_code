%% CUBEMAXWELLSADDLE solves Maxwell type equations in a cube using lowest order element.
% This is a special case of div u = g being nozero.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

clear; close all;

%% Defacult setting
[node,elem] = cubemesh([-1,1,-1,1,-1,1],2/2);
bdFlag = setboundary3(node,elem,'Dirichlet');
h = 2/2;
for k = 1:1
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    h = h/2;
end

%% pde settings
pde = data0_dual;

%%
option.printlevel = 0;
% option.solver = 'mg';
% option.solver = 'diag';

%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
ite_lst = zeros(maxIt,1);
time_lst = zeros(maxIt,1);
Vcycle_lst = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);
stepsize = 1.5; 
stepsize_Q = 0.3;
stepsize_g = .5;
ToLmg = 0.2;
% set 0.99, 0.3 and 0.3 for a1=70

beta = 1;
NAA = 10; % every NAA steps to update the preconditioner IV and IQ;
Nstart = 10; % starting from the Nstart iteration, update every NAA steps
Nmg = 200; % the steps for HodgeLaplacian mg in every iteration
maxTPDIt = 1000;
tol = 10^(-5);

resultsTPDex.stepsize = [stepsize,stepsize_Q];
resultsTPDex.NAA = NAA;
resultsTPDex.Nstart = Nstart;
resultsTPDex.Nmg = Nmg;
resultsTPDex.tol = tol;
resultsTPDex.ToLmg = ToLmg;

%% Finite Element Method        
for k = 1:maxIt   
    rng(4);

   [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag); 

    %%%%%%%%%%%%%%%%%%%%%%%
    NT = size(elem,1);
    curlE = pde.curlu;
    [lambda,w] = quadpts3(2);
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ...
            + lambda(p,4)*node(elem(:,4),:);
        if isnumeric(curlE) % a constant vector
            curlEp = repmat(curlE,NT,1);
        else % function handel
            curlEp = curlE(pxy);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%

    elemSort = sort(elem,2);
    [elem2dof,edge] = dof3edge(elemSort);
    %curluh = ones(size(elem,1),1);

    uI = edgeinterpolate(pde.exactu,node,edge);
    uold = edgeinterpolate(@(p) ones(size(p)),node,edge); % ones(size(unew));
    [Dlambda,~] = gradbasis3(node,elemSort);
    % curl phi = 2*Dlambda_i cross Dlambda_j;
    clear curlPhi
    curlPhi(:,:,6) = 2*cross(Dlambda(:,:,3),Dlambda(:,:,4),2);
    curlPhi(:,:,1) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,2),2);
    curlPhi(:,:,2) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,3),2);
    curlPhi(:,:,3) = 2*cross(Dlambda(:,:,1),Dlambda(:,:,4),2);
    curlPhi(:,:,4) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,3),2);
    curlPhi(:,:,5) = 2*cross(Dlambda(:,:,2),Dlambda(:,:,4),2);
    curlEhp = zeros(NT,3);
    for kk = 1:6
        curlEhp = curlEhp + ...
            repmat(uold(elem2dof(:,kk)),1,3).*curlPhi(:,:,kk);
    end
    % assemble A^mu
    pde.mu = pde.nu(sum(curlEhp.^2,2).^(1/2)); % coefficient

    disp('generate eqn')
    [eqn,info] = MaxwellsaddleMat(node,elem,bdFlag,pde,option);
    NE = size(eqn.edge,1);
    NT = size(elemSort,1);% Ndof = max(elem2dof(:)); %N = size(node,1);

    ii = eqn.StifMat.ii;
    jj = eqn.StifMat.jj;
    index = 0;
    sA0 = eqn.StifMat.sA;
    sA = sA0;
    for i = 1:6
        for j = i:6
            sA(index+1:index+NT) = sA0(index+1:index+NT).*pde.mu;
            index = index + NT;
        end
    end
    diagIdx = eqn.StifMat.diagIdx;
    upperIdx = eqn.StifMat.upperIdx;
    isFreeNode = eqn.freeNode;
    isFreeEdge = eqn.freeEdge;

    A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),NE,NE);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),NE,NE);
    A0 = A + AU + AU';
    clear AU

    ub = zeros(size(uI));
    ub(~isFreeEdge) = uI(~isFreeEdge);
    f0 = eqn.f0;
    f = f0 - A0*ub;
    f = f(isFreeEdge);
    A  = A0(isFreeEdge,isFreeEdge);
    G  = eqn.G(isFreeEdge,isFreeNode);
    grad = eqn.grad(isFreeEdge,isFreeNode);
    pold = zeros(size(G,2),1);
    Gt = G';
    g = eqn.g - eqn.G'*ub;
    g0 = g(isFreeNode);
     
    clear eqn
    
    bigA = [A,G; Gt,sparse(size(G,2),size(G,2))];
    resi = bigA*[uold(isFreeEdge);zeros(size(G,2),1)] - [f;g0];
    residualErr1 = zeros(maxTPDIt,1);
    residualErr1(1) = norm(resi)/sqrt(length(resi));
    residualErr2 = residualErr1;
    residualErr3 = residualErr1;
    ite = 1;
    
    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation preserved
    volume = abs(volume);
    Mvlump = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);
    NumVcycle = 0;
    NumVcycle_ite = zeros(maxTPDIt,1);

    volumemu = volume.*(pde.mu);
    Mvlumpmu = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
            [volumemu;volumemu;volumemu;volumemu]/4,[max(elem(:)),1]);
    DMd = (Mvlump(isFreeNode).^2)./Mvlumpmu(isFreeNode);
    IQd = DMd;
    DMinv = spdiags(Mvlumpmu(isFreeNode)./(Mvlump(isFreeNode).^2),0,sum(isFreeNode),sum(isFreeNode));


   % Multigrid solver settings
    option.isFreeEdge = isFreeEdge;
    option.solvermaxit = Nmg;
    option.solver = 'CG'; % 'Vcycle', 'mg', diag'
    option.coarsegridsolver = 'direct';%'direct'; % 'direct'
    option.printlevel = 0;
    option.smoother = 'GSJAC';
    %option.smoother = 'GS';
    option.smoothingstep = 10; %8 for 10, 10 for 20, 12 for 30
    option.smoothingparameter = 0.5;
    option.smoothingratio = 3;
    option.tol = ToLmg; %0.03 best

    disp('start iteration')
    tic
    while abs(residualErr1(ite))>tol*residualErr1(1) && ite < maxTPDIt
    %while abs(residualErr1(ite))> 2.92*10^(-4)

        if mod(ite,NAA) == 1 || ite <= Nstart
            AA = A + G*DMinv*Gt;
            setup_index = 1;
        else
            setup_index = 0;
        end

        %stepsize = 1.5 - 1/ite;

        du_tmp = A*uold(isFreeEdge) - f + G*(DMinv*(Gt*uold(isFreeEdge) - g0)) + G*pold;

        % Direct solver
%         unew = uold(isFreeEdge) - stepsize*((A+beta*G*Gt)\(unew_tmp)); 

        % mg solver
        if setup_index == 1
            [du,info,Ai,Bi,BBi,Res,Pro] = mgHodgeLapELarge(AA,du_tmp,node,elem,bdFlag,option);
        elseif setup_index == 0
            [du,info] = mgHodgeLapELarge(AA,du_tmp,node,elem,bdFlag,option,Ai,Bi,BBi,Res,Pro);
        end

        NumVcycle = NumVcycle +info.itStep;
        NumVcycle_ite(ite) = NumVcycle_ite(ite)+info.itStep;
        Vcycle_lst(k) = info.itStep + Vcycle_lst(k);

        unew_tmp = uold(isFreeEdge) - du;

        IQd = (1-stepsize_Q*(1-stepsize_g))*IQd + stepsize_Q*(1-stepsize_g)*DMd;
        IQinv = spdiags(IQd.^(-1),0,sum(isFreeNode),sum(isFreeNode));
        pnew = pold - stepsize_Q*IQinv*(g0 - Gt*unew_tmp);

        unew = (1-stepsize)*uold(isFreeEdge) + stepsize*unew_tmp;

        unewAll = uI;   % boundary values
        unewAll(isFreeEdge) = unew;

        % compute curl phi
        % curl phi = 2*Dlambda_i cross Dlambda_j;
        curlEhp = zeros(NT,3);
        for kk = 1:6
            curlEhp = curlEhp + ...
                repmat(unewAll(elem2dof(:,kk)),1,3).*curlPhi(:,:,kk);
        end
        % assemble A^mu
        pde.mu = pde.nu(sum(curlEhp.^2,2).^(1/2)); % coefficient

        index = 0;
        sA = sA0;
        for i = 1:6
            for j = i:6
                sA(index+1:index+NT) = sA0(index+1:index+NT).*pde.mu;
                index = index + NT;
            end
        end
        A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),NE,NE);
        AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),NE,NE);
        A0 = A + AU + AU';
        A  = A0(isFreeEdge,isFreeEdge);
        clear AU

        f = f0 - A0*ub;
        f = f(isFreeEdge);

        volumemu = volume.*(pde.mu);
        Mvlumpmu = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
            [volumemu;volumemu;volumemu;volumemu]/4,[max(elem(:)),1]);
        DMd = (Mvlump(isFreeNode).^2)./Mvlumpmu(isFreeNode);
        DMinv = spdiags(Mvlumpmu(isFreeNode)./(Mvlump(isFreeNode).^2),0,sum(isFreeNode),sum(isFreeNode));

        uold = uI; 
        uold(isFreeEdge) = unew;
        pold = pnew;

        ite = ite+1;
        bigA = [A,G;G',sparse(size(G,2),size(G,2))];
        resi = bigA*[unew;pnew] - [f;g0];
%         resi1 = resi(1:sum(isFreeEdge));
%         resi2 = resi(sum(isFreeEdge)+1:end);
        %residualErr(ite) = norm(resi);
        %residualErr(ite) = max(abs(resi));
        residualErr1(ite) = norm(resi)/sqrt(length(resi));
        % max(norm(resi1)/sqrt(length(resi1)),norm(resi2)/sqrt(length(resi2)));
        residualErr2(ite) = sqrt(sum(du.*du_tmp)/length(du_tmp));
        residualErr3(ite) = max(abs(resi));
        disp(residualErr1(ite));
    end
    clear curlPhi;
    option = struct();
    u = unewAll;
    fprintf('\n ite stp: = %d \n',ite - 1);

    timeTotal = toc;
    time_lst(k) = timeTotal;
    fprintf('\n time: = %g \n',timeTotal);

    resultsTPDex.residualErr1 = residualErr1;
    resultsTPDex.residualErr2 = residualErr2;
    resultsTPDex.residualErr3 = residualErr3;
    resultsTPDex.u = uold;
    resultsTPDex.ite = ite;
    resultsTPDex.timeTotal = timeTotal;
    resultsTPDex.NumVcycle = NumVcycle;
    resultsTPDex.NumVcycle_ite = NumVcycle_ite;
    resultsTPDex.NormType = 'max';
    resultsTPDex.stp = 'residualErr1';
    %save('results','results');

    % compute error
    energyErr(k) = getHcurlerror3ND(node,elem,pde.curlu,u);
    L2Err(k) = getL2error3ND(node,elem,pde.exactu,u);
    uIuhErr(k) = sqrt((u-uI)'*A0*(u-uI));        
    fprintf('||curl(u-u_h)|| is %g \n',energyErr(k))
    N(k) = length(unewAll) + size(node,1);
    h(k) = 1./(size(node,1)^(1/3)-1);   
    fprintf('\n # of DoFs = %d \n',N(k));
    resultsTPDex.L2Err = L2Err(1);
    resultsTPDex.energyErr = energyErr(1);
    ite_lst(k) = ite-1;
end

% %% Plot convergence rates
% figure(1);
% showrateh3(h,energyErr,1,'k-+','|| curl (u-u_h) ||',...
%            h,uIuhErr,1,'r-+','|| curl (u_I-u_h) ||',...
%            h,L2Err,1,'b-+','|| u-u_h||');