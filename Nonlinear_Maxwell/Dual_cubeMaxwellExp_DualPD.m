%% Primal-Dual Method for nonlinear ferromagnetism example

clear; close all;

%% generate mesh
[node,elem, HB] = cubemesh([-1,1,-1,1,-1,1],2/2);
bdFlag = setboundary3(node,elem,'Dirichlet');
h = 2/2;
for k = 1:1
    [node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB);
    h = h/2;
end
elemunSort = elem;
bdFlagunSort = bdFlag;

%% pde settings
pde = data1;
stepsize = 0.5;

%% Parameters
maxIt = 3; 
maxNewtonIt = 500;
tol = 10^(-6);
TolFactor = 10^(-1);

N = zeros(maxIt,1); 
h_lst = zeros(maxIt,1);
energyErr = zeros(maxIt,1);
pErr = zeros(maxIt,1);
sigmaErr = zeros(maxIt,1);
L2Err = zeros(maxIt,1);
uIuhErr = zeros(maxIt,1);
Vcycle_num = zeros(maxIt,1);
Vcycle_lst = zeros(maxIt,maxNewtonIt);
ite_lst=  zeros(maxIt,1);
time_lst=  zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt   
    rng(4)

    [node,elemunSort,bdFlagunSort,HB] = uniformrefine3(node,elemunSort,bdFlagunSort,HB);
    [elem,bdFlag] = sortelem3(elemunSort,bdFlagunSort);
    h = h/2;

    [elem2edge,edge] = dof3edge(elem);
    [Dlambda,~] =gradbasis3(node,elem);

    NT = size(elem,1);
    NE = size(edge,1);
    Nn = size(node,1);

    % Mass matrix in H(div) -- Msigma
    px = accumarray(repmat((1:NT)', 4, 1),...
        [node(elem(:,1),1);node(elem(:,2),1);node(elem(:,3),1);node(elem(:,4),1)]/4,[NT,1]);
    py = accumarray(repmat((1:NT)', 4, 1),...
        [node(elem(:,1),2);node(elem(:,2),2);node(elem(:,3),2);node(elem(:,4),2)]/4,[NT,1]);
    pz = accumarray(repmat((1:NT)', 4, 1),...
        [node(elem(:,1),3);node(elem(:,2),3);node(elem(:,3),3);node(elem(:,4),3)]/4,[NT,1]);


    % exact nodal values
    uI = edgeinterpolate(pde.exactu,node,edge);
    sigmaI = pde.exactsigma([px py pz]);
    uold = edgeinterpolate(@(p) ones(size(p)),node,edge); % random initial value

    % Element-wise basis
    % edge indices of 6 local bases: 
    % [1 2], [1 3], [1 4], [2 3], [2 4], [3 4]
    % phi = lambda_iDlambda_j - lambda_jDlambda_i;
    % curl phi = 2*Dlambda_i x Dlambda_j;
    curlPhi = zeros(NT, 3, 6);
    curlPhi(:,:,6) = 2*mycross(Dlambda(:,:,3),Dlambda(:,:,4),2);
    curlPhi(:,:,1) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,2),2);
    curlPhi(:,:,2) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,3),2);
    curlPhi(:,:,3) = 2*mycross(Dlambda(:,:,1),Dlambda(:,:,4),2);
    curlPhi(:,:,4) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,3),2);
    curlPhi(:,:,5) = 2*mycross(Dlambda(:,:,2),Dlambda(:,:,4),2);
    curlx = sparse(repmat((uint32(1:NT))',6,1), elem2edge(:) ,...
            [curlPhi(:,1,1); curlPhi(:,1,2); curlPhi(:,1,3);curlPhi(:,1,4); ...
            curlPhi(:,1,5); curlPhi(:,1,6)], NT, NE);
    curly = sparse(repmat((uint32(1:NT))',6,1), elem2edge(:) ,...
            [curlPhi(:,2,1); curlPhi(:,2,2); curlPhi(:,2,3);curlPhi(:,2,4); ...
            curlPhi(:,2,5); curlPhi(:,2,6)], NT, NE);
    curlz = sparse(repmat((uint32(1:NT))',6,1), elem2edge(:) ,...
            [curlPhi(:,3,1); curlPhi(:,3,2); curlPhi(:,3,3);curlPhi(:,3,4); ...
            curlPhi(:,3,5); curlPhi(:,3,6)], NT, NE);
    curl = [curlx; curly; curlz];    

    curlEhp = zeros(NT,3);
    for kk = 1:6
        curlEhp = curlEhp + ...
            repmat(uold(elem2edge(:,kk)),1,3).*curlPhi(:,:,kk);
    end
    pde.mu = pde.nu(sum(curlEhp.^2,2).^(1/2)); % coefficient

    disp('generate eqn')
    option = [];
    [eqn,~] = MaxwellsaddleMat(node,elem,bdFlag,pde,option);
    isFreeNode = eqn.freeNode;
    isFreeEdge = eqn.freeEdge;

    volume = abs(simplexvolume(node,elem)); % uniform refine in 3D is not orientation preserved
    Msigma = [volume; volume; volume];
    DMn = accumarray([elem(:,1);elem(:,2);elem(:,3);elem(:,4)],...
                        [volume;volume;volume;volume]/4,[max(elem(:)),1]);
    DMninv = 1./DMn;
    DMninvfree =  DMninv.*isFreeNode;

    % assemble A^0
    ii = eqn.StifMat.ii;
    jj = eqn.StifMat.jj;
    sA0 = eqn.StifMat.sA;
    diagIdx = eqn.StifMat.diagIdx;
    upperIdx = eqn.StifMat.upperIdx;

    A = sparse(ii(diagIdx),jj(diagIdx),sA0(diagIdx),NE,NE);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA0(upperIdx),NE,NE);
    A0 = A + AU + AU';
    clear AU
    
    % assemble A^mu
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

    A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),NE,NE);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),NE,NE);
    A = A + AU + AU';
    clear AU


    % incidence matrix of vertex-edge
    C = Msigma.*curl;
    Ct = C';

    % grad matrix
    grad = icdmat(double(edge),[-1,1]);
    gradft = grad(isFreeEdge,isFreeNode)';
    G = eqn.G;
    Gt = G';
    BBt = grad*Gt;




    % inital values and boundary values
    ub = zeros(size(uI));
    ub(~isFreeEdge) = uI(~isFreeEdge);
    uold = edgeinterpolate(@(p) ones(size(p)),node,edge);  
    uold(~isFreeEdge) = uI(~isFreeEdge);
    ubar = uold;


    sigma = ones(3*NT,1);

    p = zeros(Nn,1);
    p(~isFreeNode) = 0;

    curlu = curl*uold;
    curlunorm = sum(reshape(curlu, NT, 3).^2, 2).^(1/2);


    sigmanorm = zeros(NT,1);
    sigmamat = reshape(sigma, NT, 3);
    sigmal2norm = sum(sigmamat.^2, 2).^(1/2);
    z = pde.Phi_inv(sigmal2norm);
    K = pde.nu(z);
    volumemu = volume.*K;
    Jsigma = pde.Jsigma(sigma, K, z);
    Jsigmainv = pde.Jsigmainv(sigma, K, z);
    MJsigmainv = Jsigmainv./[volume;volume;volume];
     

    %  assemble residual

    ite = 1;
  
    resisigma = Msigma.*(sigma./[K; K; K]) - C*uold;
    resip = DMn.*p - Gt*uold + eqn.gall;
    resiu = Ct*sigma - eqn.f0 + G*p;
    resi = [resisigma; resip(isFreeNode); resiu(isFreeEdge)];
    residualErr1 = zeros(maxIt,1);
    residualErr1(1) = norm(resi)/sqrt(length(resi));
    

    % Dual PD iteraion
    disp('solver begins')
    tic;
    while residualErr1(ite)>tol*residualErr1(1) && ite <= maxNewtonIt

        % mg 
        rsigma = Msigma.*(sigma./[K; K; K]) - C*ubar;
        rp = DMn.*p - Gt*ubar+ eqn.gall;

        % % Multigrid solver settings
        option.printlevel = 0;
        option.N0 = 27;
        option.tol = 1e-1;

        sigma = sigma - stepsize*MJsigmainv*rsigma; 
        dp =  DMninvfree.*rp; 
        p(isFreeNode) = p(isFreeNode) - stepsize*dp(isFreeNode);


        % update Jacobian
        sigmanorm = zeros(NT,1);
        sigmamat = reshape(sigma, NT, 3);
        sigmal2norm = sum(sigmamat.^2, 2).^(1/2);
        z = pde.Phi_inv(sigmal2norm);
        K = pde.nu(z);
        volumemu = volume.*K;
        Jsigma = pde.Jsigma(sigma, K, z);
        Jsigmainv = pde.Jsigmainv(sigma, K, z);
        MJsigmainv = Jsigmainv./[volume;volume;volume];
        
        ru =  Ct*sigma - eqn.f0 + G*p;
        AJ = Ct*(Jsigmainv*curl);
        AA = (AJ + G*(DMninv.*Gt));
    
        bdidx = zeros(NE,1); 
        bdidx(~isFreeEdge) = 1;
        Tbd = spdiags(bdidx,0,NE,NE);
        Te = spdiags(1-bdidx,0,NE,NE);
        bigAA = Te*AA*Te + Tbd;
        ru(~isFreeEdge) = 0;
        
        [du,info] = mgMaxwell_1(bigAA,ru,node,elemunSort,edge,HB,option); 
        Vcycle_num(k) = Vcycle_num(k) + info.itStep;
        
        % relaxation 
        unew = zeros(size(uold));
        unew(isFreeEdge) = uold(isFreeEdge) - stepsize*du(isFreeEdge);
        ubar(isFreeEdge) = unew(isFreeEdge) + 0.8*(unew(isFreeEdge) - uold(isFreeEdge));
        uold(isFreeEdge) = unew(isFreeEdge);


        resisigma = Msigma.*(sigma./[K; K; K]) - C*uold;
        resip = DMn.*p - Gt*uold + eqn.gall;
        resiu = Ct*sigma - eqn.f0 + G*p;
        resi = [resisigma; resip(isFreeNode); resiu(isFreeEdge)];

        ite = ite+1;
        resi = resi(isFreeEdge);
        residualErr1(ite) = norm(resi)/sqrt(length(resi));       
        % disp(residualErr1(ite));
    end
    
    clear curlPhi;
    u = uold;
    fprintf('\n ite stp: = %d \n',ite);

    timeTotal = toc;
    fprintf('\n time: = %g \n',timeTotal);

  
    % compute error
    energyErr(k) = getHcurlerror3ND(node,elem,pde.curlu,u);
    sigmaErr(k) = getL2error3(node,elem,pde.exactsigmax,sigma(1:NT));
    L2Err(k) = getL2error3ND(node,elem,pde.exactu,u);
    uIuhErr(k) = sqrt((u-uI)'*(A0*(u-uI)));    
    pErr(k) = getL2error3(node,elem,pde.exactp,p);
    % fprintf('||curl(u-u_h)|| is %g \n',energyErr(k))
    N(k) = length(u) + length(sigma);
    % fprintf('\n # of DoFs = %d \n',N(k));

    h_lst(k) = h;   
    ite_lst(k) = ite;
    time_lst(k) = timeTotal;
end



%% Print iteration results
fprintf('\nIteration Results\n');
fprintf('---------------------------------------------------------\n');
iterationResults = table(h_lst, ite_lst, Vcycle_num, time_lst, ...
    'VariableNames', {'h', 'Iteration', 'Vcycle', 'Time'});
disp(iterationResults);
