function [sigma,u,eqn,info] = pLapP0P1_DualTPDM(node,elem,bdFlag,pde,option)


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1); 

maxIt = 3e4;
tol = 1e-6;
flag = 0;
stepsize = 1.5;

%% Data structure
[~,~,isBdNode] = findboundary(elem,bdFlag);
isfreeNode = ~isBdNode;

NT = size(elem,1);
Nu = N; Nsigma = 2*NT; Ndof = Nu + Nsigma;
Nufree = sum(isfreeNode);

% exact solutions
uI = pde.exactu(node);
sigmaI = Lagrangeinterpolate(pde.exactsigma,node,elem,'P0');
ub = zeros(Nu,1);
ub(~isfreeNode) = uI(~isfreeNode);

%% Initial values

% uoldAll = uI;
% sigmaoldAll = sigmaI;
% % 
% rng(4);
% uoldAll = rand(Nu,1);
% sigmaoldAll = rand(NT,2);

uoldAll = zeros(Nu,1);
sigmaoldAll = zeros(NT,2);

uoldAll(~isfreeNode) = uI(~isfreeNode);


%% fixed Diffusion coefficient
if ~isfield(pde,'K')
    pde.K = []; 
    K = [];
end
if ~isempty(pde.K) && isnumeric(pde.K)
    K = pde.K;                                 % d is an array
end
if ~isfield(option,'dquadorder'), option.dquadorder = 1; end

%% Assemble matrix for divergence operator
[Dlambda,area,~] = gradbasis(node,elem);
d1 = Dlambda(:,:,1).*[area, area];
d2 = Dlambda(:,:,2).*[area, area];
d3 = Dlambda(:,:,3).*[area, area];
Dx = sparse(elem(:), repmat((1:NT)',3,1),...
            [d1(:,1); d2(:,1); d3(:,1)], Nu, NT);
Dy = sparse(elem(:), repmat((1:NT)',3,1),...
            [d1(:,2); d2(:,2); d3(:,2)], Nu, NT);
B = [Dx Dy];  %this is grad transpose with mass
Bt = B';
clear d1 d2 d3

% [A,~,~] = assemblematrix(node,elem,1, []);
% A = A(isfreeNode, isfreeNode);

% Assemble right hand side.
b = zeros(N,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 3;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if isreal(pde.f) % f is a real number or vector and not a function
   switch length(pde.f)
       case NT  % f is piecewise constant
         bt = pde.f.*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[Ndof 1]);
       case N   % f is piecewise linear
         bt = zeros(NT,3);
         bt(:,1) = area.*(2*pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/12;
         bt(:,2) = area.*(2*pde.f(elem(:,2)) + pde.f(elem(:,3)) + pde.f(elem(:,1)))/12;
         bt(:,3) = area.*(2*pde.f(elem(:,3)) + pde.f(elem(:,1)) + pde.f(elem(:,2)))/12;
         b = accumarray(elem(:),bt(:),[N 1]);
       case 1   % f is a scalar e.g. f = 1
         bt = pde.f*area/3;
         b = accumarray(elem(:),[bt; bt; bt],[N 1]);
   end
end
if ~isempty(pde.f) && ~isreal(pde.f)  % f is a function 
    [lambda,weight] = quadpts(option.fquadorder);
    phi = lambda;                 % linear bases
	nQuad = size(lambda,1);
    bt = zeros(NT,3);
    for p = 1:nQuad
		% quadrature points in the x-y coordinate
		pxy = lambda(p,1)*node(elem(:,1),:) ...
			+ lambda(p,2)*node(elem(:,2),:) ...
			+ lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxy);
        for i = 1:3
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*fp;
        end
    end
    bt = bt.*repmat(area,1,3);
    b = accumarray(elem(:),bt(:),[Nu 1]);
end
clear pxy bt

F = zeros(2*NT + Nu, 1); % all righ hand side
F(2*NT+1:end, 1) = b;

g = zeros(Nsigma,1); % the right hand side of sigma
F(1:2*NT, 1) = g;


% Assemble Mass Matrix
if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
    sigmanorm = sqrt(sum(sigmaoldAll.^2, 2));
    K = pde.K(sigmanorm); % K is a scalar
    % Keps = pde.K(sigmanorm + keps); 
end

% M. Mass matrix for P0 element
M = [K.*area;K.*area];
% Meps = [Keps.*area;Keps.*area];

eqn = struct('M',M,'B',B,...
             'f',F(2*NT+1:end),'g',F(1:2*NT),...
             'Bff',B(isfreeNode,:),'Bfft',B(isfreeNode,:)', 'ff',b(isfreeNode));


% Jinv = pde.Juinv(uoldAll);
% MJinv = [1./area; 1./area].*Jinv;
% Minv = 1./M;
% S = eqn.Bff*( Minv.*(eqn.Bfft));
% % S = eqn.Bff*( MJinv*(eqn.Bfft));
% bigS = [S sparse(Nufree,N-Nufree); sparse(N-Nufree,Nufree) speye(N-Nufree)];
% 
% % set up multilevel structure for mg solver
% setupOption.solver = 'NO';
% setupOption.freeDof =  [true(Nufree, 1); false];
% [~,~,Ai,Bi,BBi,Res,Pro,isFreeDof] = mg(bigS,uoldAll,elem,setupOption);

% A = [M (B(1:Np-1, :))'; B(1:Np-1, :) sparse(Np-1, Np-1)]; % direct solver
% bigu = A\F(1:end-1);
% unewAll = bigu(1:Nu);
% pnewAll = bigu(Nu+1:end);
% pnewAll(Np) = 0;


%% Error
res = zeros(maxIt,1);
Btu =  Bt*uoldAll;
% % res0  = norm([M.*sigmaoldAll(:) - Btu - eqn.g; eqn.Bff*sigmaoldAll(:)  - eqn.ff]);
res0 = norm([eqn.g; eqn.ff]);
% res0 = norm(eqn.Bff*sigmaoldAll(:)  - eqn.ff);
%% Solve the saddle point problem
Vcycle_num = 0;

tic
for ite = 1:maxIt
%     fprintf('Iteration %i \n', ite);
     % Newton's method
    
     rsigma = M.*sigmaoldAll(:) - Btu - eqn.g;
     ru = eqn.Bff*sigmaoldAll(:)  - eqn.ff;

     % J = pde.Jsigma(sigmaoldAll);
     % MJ = [area; area].*J;
     % 

     % Minv = 1./M;

     % 
     Meps = M + 0.0001*[area;area];
     Minv = 1./Meps;

    % small p as preconditioner
    %  if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
    %     sigmanorm = sqrt(sum(sigmaoldAll.^2, 2));
    %     Ksmallp = pde.Ksmallp(sigmanorm); % K is a scalar
    %  end
    %  Msmall = [Ksmallp.*area; Ksmallp.*area];
    % Minv = 1./Msmall;

     % 
     % Jinv = pde.Jsigmainv(sigmaoldAll);
     % MJinv = [1./area; 1./area].*Jinv;


     % Schur complement
     S = eqn.B*( Minv.*(Bt));

     bdidx = zeros(N,1); 
     bdidx(~isfreeNode) = 1;
     Tbd = spdiags(bdidx,0,N,N);
     Te = spdiags(1-bdidx,0,N,N);
     bigS = Te*S*Te + Tbd;
     
     ruall = zeros(N,1);
     ru =  ru - eqn.Bff*( Minv.*rsigma);
     ruall(isfreeNode) = ru;

%      % multigrid
     % [Ai, Bi, BBi] = mgmat(S, Ai, Res, Pro);
     % [duAll, info] = mg(bigS, ruall, elem, option.mg,Ai,Bi,BBi,Res,Pro,isFreeDof);
     [duAll, info] = mg(bigS, ruall, elem, option.mg);
     Vcycle_num  = Vcycle_num  + info.itStep;

     dsigma = Minv.*(rsigma+eqn.Bfft*duAll(isfreeNode));
     sigmanewAll = sigmaoldAll - stepsize*reshape(dsigma,NT, 2);
     unewAll =  uoldAll;   % boundary values
     unewAll(isfreeNode) =  uoldAll(isfreeNode) - stepsize*duAll(isfreeNode);


   % M. Mass matrix for P0 element
    if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
         sigmanorm = sqrt(sum(sigmanewAll.^2, 2));
          K = pde.K(sigmanorm); % K is a scalar
          % Keps = pde.K(sigmanorm + keps); 
    end
    M = [K.*area;K.*area];
    % Meps = [Keps.*area; Keps.*area];
    eqn.M = M;

    % compute error
    % Btu =  eqn.Bfft*unewAll(isfreeNode);
    Btu =  Bt*unewAll;
    res(ite) = norm([M.*sigmanewAll(:) - Btu - eqn.g; eqn.Bff*sigmanewAll(:)  - eqn.ff]);
     % res(ite) = norm([eqn.Bff*sigmanewAll(:)  - eqn.ff]);
    disp(res(ite))
     
    if ite> 1 && res(ite) < tol*res0
        flag = 1;
        break;
    end

    sigmaoldAll = sigmanewAll;
    uoldAll = unewAll;


%% 
end
time = toc;



% Assemble Mass Matrix
if ~isempty(pde.K) && ~isnumeric(pde.K)       % d is a function   
    sigmanorm = sqrt(sum(sigmaI.^2, 2));
    K = pde.K(sigmanorm); % K is a scalar
end

% M. Mass matrix for P0 element
M = [K.*area;K.*area];

% x = [spdiags(M, 0, size(M,1), size(M,1)) Bt; eqn.Bff sparse(Nufree, Nu)]\[eqn.g; eqn.ff];
% sigmdirect = reshape(x(1:2*NT), NT, 2);

resexact = [spdiags(M, 0, size(M,1), size(M,1)) Bt; eqn.Bff sparse(Nufree, Nu)]*[sigmaI(:); uI] -[eqn.g; eqn.ff];
norm(resexact)


%% Output
info = struct('itStep', ite ,'res',res(1:ite), 'flag', flag, 'Vcycle_num', Vcycle_num, 'time', time);
% plot(1:ite, log10(res(1:ite)));
sigma = sigmanewAll;
u = unewAll;
end

