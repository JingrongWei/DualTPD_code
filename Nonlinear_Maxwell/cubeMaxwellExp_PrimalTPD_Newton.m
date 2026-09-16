%% CUBEMAXWELLEXP_PRIMALTPD_NEWTON solves the primal nonlinear Maxwell saddle system.
% The Newton system is augmented by beta*G*D_p^{-1}*G' and solved by GMRES.

clear; close all;

%% Mesh and PDE data
[node,elem,HB] = cubemesh([-1,1,-1,1,-1,1],2/2);
bdFlag = setboundary3(node,elem,'Dirichlet');
for j = 1:1
    [node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB);
end
elemunSort = elem;
bdFlagunSort = bdFlag;
pde = data0_dual;

%% Parameters
maxRefine = 3;
maxNewtonIt = 100;
Newtonstep = 0.5;
tol = 1e-6;
beta = 1;
gmresRestart = 20;
gmresTol = 1e-4;
gmresMaxIt = 60;

N = zeros(maxRefine,1);
h_lst = zeros(maxRefine,1);
ite_lst = zeros(maxRefine,1);
time_lst = zeros(maxRefine,1);
energyErr = zeros(maxRefine,1);
L2Err = zeros(maxRefine,1);
pErr = zeros(maxRefine,1);
Vcycle_num = zeros(maxRefine,1);

for k = 1:maxRefine
    [node,elemunSort,bdFlagunSort,HB] = ...
        uniformrefine3(node,elemunSort,bdFlagunSort,HB);
    [elem,bdFlag] = sortelem3(elemunSort,bdFlagunSort);

    [elem2edge,edge] = dof3edge(elem);
    NT = size(elem,1);
    NE = size(edge,1);
    [Dlambda,~] = gradbasis3(node,elem);

    % Discrete curl from edge degrees of freedom to elementwise vectors.
    curlPhi = zeros(NT,3,6);
    curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
    curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
    curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
    curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
    curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
    curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);
    curlx = sparse(repmat((1:NT)',6,1),double(elem2edge(:)),...
        [curlPhi(:,1,1);curlPhi(:,1,2);curlPhi(:,1,3);curlPhi(:,1,4);...
         curlPhi(:,1,5);curlPhi(:,1,6)],NT,NE);
    curly = sparse(repmat((1:NT)',6,1),double(elem2edge(:)),...
        [curlPhi(:,2,1);curlPhi(:,2,2);curlPhi(:,2,3);curlPhi(:,2,4);...
         curlPhi(:,2,5);curlPhi(:,2,6)],NT,NE);
    curlz = sparse(repmat((1:NT)',6,1),double(elem2edge(:)),...
        [curlPhi(:,3,1);curlPhi(:,3,2);curlPhi(:,3,3);curlPhi(:,3,4);...
         curlPhi(:,3,5);curlPhi(:,3,6)],NT,NE);
    curl = [curlx;curly;curlz];

    uI = edgeinterpolate(pde.exactu,node,edge);
    u = edgeinterpolate(@(p) ones(size(p)),node,edge);
    curlu = curl*u;
    curlunorm = sqrt(sum(reshape(curlu,NT,3).^2,2));
    pde.mu = pde.nu(curlunorm);
    option = struct('printlevel',0);
    [eqn,~] = MaxwellsaddleMat(node,elem,bdFlag,pde,option);

    isFreeEdge = eqn.freeEdge;
    isFreeNode = eqn.freeNode;
    G = eqn.G(isFreeEdge,isFreeNode);
    Gt = G';
    u(~isFreeEdge) = uI(~isFreeEdge);
    p = zeros(sum(isFreeNode),1);
    g = eqn.g(isFreeNode) - eqn.G(~isFreeEdge,isFreeNode)'*u(~isFreeEdge);

    volume = abs(simplexvolume(node,elem));
    DpAll = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
        [volume;volume;volume;volume]/4,[size(node,1),1]);
    Dp = DpAll(isFreeNode);
    DpInv = 1./Dp;
    Mvol = spdiags([volume;volume;volume],0,3*NT,3*NT);

    mgOption = struct('printlevel',0,'N0',27,'tol',gmresTol,...
        'outsolver','CG','isBdEdge',~isFreeEdge);
    bdidx = double(~isFreeEdge);
    Tbd = spdiags(bdidx,0,NE,NE);
    Te = spdiags(1-bdidx,0,NE,NE);

    [rU,rP,Jaug,Jfree] = assembleNewtonData(u,p,curl,volume,pde,eqn,...
        isFreeEdge,isFreeNode,Gt,g,DpInv,beta,Mvol,Te,Tbd,NT);
    residual = [rU;rP];
    residual0 = norm(residual);
    if residual0 == 0, residual0 = 1; end
    residualErr = zeros(maxNewtonIt+1,1);
    residualErr(1) = residual0/sqrt(length(residual));
    totalVcycles = 0;

    tic
    for ite = 1:maxNewtonIt
        if norm(residual) < tol*residual0
            break;
        end

        rhs = [rU + beta*G*(DpInv.*rP); rP];
        preconditioner = @(r) blockDiagonalPreconditioner(r,Jaug,NE,...
            isFreeEdge,Dp,beta,node,elemunSort,edge,HB,mgOption);
        [delta,flag,~,iter] = gmres(@(x) augmentedNewtonOperator(x,Jfree,G,Gt),...
            rhs,gmresRestart,gmresTol,gmresMaxIt,preconditioner);
        if flag ~= 0
            warning('Newton GMRES did not reach its requested tolerance.');
        end
        totalVcycles = totalVcycles + max(0,(iter(1)-1)*gmresRestart + iter(2));

        u(isFreeEdge) = u(isFreeEdge) - Newtonstep*delta(1:sum(isFreeEdge));
        p = p - Newtonstep*delta(sum(isFreeEdge)+1:end);
        [rU,rP,Jaug,Jfree] = assembleNewtonData(u,p,curl,volume,pde,eqn,...
            isFreeEdge,isFreeNode,Gt,g,DpInv,beta,Mvol,Te,Tbd,NT);
        residual = [rU;rP];
        residualErr(ite+1) = norm(residual)/sqrt(length(residual));
        fprintf('Newton iteration %d: residual = %.4e\n',ite,residualErr(ite+1));
    end
    time_lst(k) = toc;
    ite_lst(k) = ite - (norm(residual) < tol*residual0);
    Vcycle_num(k) = totalVcycles;

    energyErr(k) = getHcurlerror3ND(node,elem,pde.curlu,u);
    L2Err(k) = getL2error3ND(node,elem,pde.exactu,u);
    pAll = zeros(size(node,1),1);
    pAll(isFreeNode) = p;
    pErr(k) = getL2error3(node,elem,pde.exactp,pAll);
    N(k) = length(u) + length(p);
    h_lst(k) = 1/(size(node,1)^(1/3)-1);

    fprintf('Refinement %d: Newton iterations = %d, time = %.3g s\n',...
        k,ite_lst(k),time_lst(k));
end

function [rU,rP,Jaug,Jfree] = assembleNewtonData(u,p,curl,volume,pde,eqn,...
    isFreeEdge,isFreeNode,Gt,g,DpInv,beta,Mvol,Te,Tbd,NT)
    curlCurrent = curl*u;
    curlNorm = sqrt(sum(reshape(curlCurrent,NT,3).^2,2));
    nuCurrent = pde.nu(curlNorm);
    flux = [volume;volume;volume].*[nuCurrent;nuCurrent;nuCurrent].*curlCurrent;
    nonlinearResidual = curl'*flux - eqn.f0 + eqn.G(:,isFreeNode)*p;
    rU = nonlinearResidual(isFreeEdge);
    rP = Gt*u(isFreeEdge) - g;
    
    % Ju is the elementwise derivative of nu(|curl u|) curl u.
    Ju = pde.Ju(curlCurrent,max(curlNorm,eps));
    Jcurl = curl'*(Mvol*Ju)*curl;
    Gp = eqn.G(:,isFreeNode);
    JaugFull = Jcurl + beta*Gp*spdiags(DpInv,0,length(DpInv),length(DpInv))*Gp';
    Jfree = JaugFull(isFreeEdge,isFreeEdge);
    Jaug = Te*JaugFull*Te + Tbd;
end

function y = augmentedNewtonOperator(x,Jfree,G,Gt)
    nu = size(G,1);
    y = [Jfree*x(1:nu) + G*x(nu+1:end);...
         Gt*x(1:nu)];
end

function z = blockDiagonalPreconditioner(r,Jaug,NE,isFreeEdge,Dp,beta,node,elem,edge,HB,option)
    ru = zeros(NE,1);
    ru(isFreeEdge) = r(1:sum(isFreeEdge));
    [zu,~] = mgMaxwell_1(Jaug,ru,node,elem,edge,HB,option);
    z = [zu(isFreeEdge); r(sum(isFreeEdge)+1:end)./(beta*Dp)];
end
