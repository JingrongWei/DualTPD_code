function [sigma,u,eqn,info] = pLapP0P1_primal(node,elem,bdFlag,pde,option)


if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1); 

maxIt = 3e4;
tol = 1e-6;
flag = 0;
stepsize = eps;
keps = 0.0001;

%% Data structure
[~,~,isBdNode] = findboundary(elem,bdFlag);
isfreeNode = ~isBdNode;

NT = size(elem,1);
Nu = N; Ndof = Nu;

% exact solutions
uI = pde.exactu(node);
ub = zeros(Nu,1);
ub(~isfreeNode) = uI(~isfreeNode);

%% Initial values

% uoldAll = uI;
% sigmaoldAll = sigmaI;

% uoldAll = zeros(Nu,1);

rng(4); 
uoldAll = rand(Nu,1);

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
grad = Bt./[area;area];
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

f = b;

% Assemble Mass Matrix
if ~isempty(pde.Ku) && ~isnumeric(pde.Ku)       % d is a function
    gradu = reshape(grad*uoldAll, NT, 2);
    gradunorm = sqrt(sum(gradu.^2, 2));
    K = pde.Ku(gradunorm ); % K is a scalar
end


% M = accumarray([elem(:,1);elem(:,2);elem(:,3)],[area;area;area]/3,[N,1]);
% A = B*(1./[area;area].*Bt) + sparse(1:Nu, 1:Nu, 0.0001*M);
% 
% bdidx = zeros(N,1); 
% bdidx(~isfreeNode) = 1;
% Tbd = spdiags(bdidx,0,N,N);
% Te = spdiags(1-bdidx,0,N,N);
% bigA = Te*A*Te + Tbd;


Ap = B*([K; K]./[area;area].*Bt);
Apeps = B*(([K; K]+keps)./[area;area].*Bt);
bdidx = zeros(N,1); 
bdidx(~isfreeNode) = 1;
Tbd = spdiags(bdidx,0,N,N);
Te = spdiags(1-bdidx,0,N,N);
bigAp = Te*Apeps*Te + Tbd;

eqn = struct('Ap',Ap,'B',B,...
             'f',f,...
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
% res0 = Ap*uoldAll - eqn.f;
% res0  = norm(res0(isfreeNode));
res0  = norm(eqn.f(isfreeNode));

energy0 =  sum(gradunorm.^(pde.coff_p).*area)/(pde.coff_p) - f'*uoldAll;
energy = zeros(maxIt,1);

stepsize_lst = zeros(maxIt,1);

%% Solve the saddle point problem
Vcycle_num = 0;


tic
for ite = 1:maxIt
%     fprintf('Iteration %i \n', ite);

     % Newton's method
     ruall = zeros(N,1);
     ru =  Ap*uoldAll - eqn.f;
     ruall(isfreeNode) = ru(isfreeNode);

%      % multigrid
     % [Ai, Bi, BBi] = mgmat(S, Ai, Res, Pro);
     % [duAll, info] = mg(bigS, ruall, elem, option.mg,Ai,Bi,BBi,Res,Pro,isFreeDof);
     [duAll, info] = mg(bigAp, ruall, elem, option.mg);
     Vcycle_num  = Vcycle_num  + info.itStep;


     du = zeros(N,1);
     du(isfreeNode) = -duAll(isfreeNode); % negative gradient, decent direction

     % if ite > 100
     %     stepsize = 1.05*stepsize;
     % end


     % --------
     % line search

     % bisection 
     a = eps;
     b = max(stepsize, eps);

     gradudu = reshape(grad*(uoldAll+b*du), NT, 2);
     gradudunorm = sqrt(sum(gradudu.^2, 2));
     K = pde.Ku(gradudunorm);
     Apb = B*([K; K]./[area;area].*Bt);
     Fb = du'*(Apb*(uoldAll+b*du)) - du'*f;

     while Fb < 0
         a = b;
         b = 2*b;

         gradudu = reshape(grad*(uoldAll+b*du), NT, 2);
         gradudunorm = sqrt(sum(gradudu.^2, 2));
         K = pde.Ku(gradudunorm); % K is a scalar
         Apb = B*([K; K]./[area;area].*Bt);
         Fb = du'*(Apb*(uoldAll+b*du)) - du'*f;     
     end


     gradudu = reshape(grad*(uoldAll+a*du), NT, 2);
     if ~isempty(pde.Ku) && ~isnumeric(pde.Ku)       % d is a function
        gradudunorm = sqrt(sum(gradudu.^2, 2));
        K = pde.Ku(gradudunorm); % K is a scalar
     end
     Apa = B*([K; K]./[area;area].*Bt);
     Fa =du'*(Apa*(uoldAll+a*du)) - du'*f;

     c = (b+a)/2;
     in_ite = 0;

     while b~= a && in_ite < 20

         gradudu = reshape(grad*(uoldAll+c*du), NT, 2);
         gradudunorm = sqrt(sum(gradudu.^2, 2));
         K = pde.Ku(gradudunorm); % K is a scalar
         Apc = B*([K; K]./[area;area].*Bt);
         Fc = du'*(Apc*(uoldAll+c*du)) - du'*f;

         if abs(Fc) < 10^(-12)
             break
         end

         if Fa*Fc > 0
             a = c;
             Fa = Fc;
         else
             b = c;
             Fb = Fc;
         end

         c = (b+a)/2;
         in_ite = in_ite + 1;

         if ite > 1 && in_ite > 8
             break;
         end
     end
     stepsize = a;
     stepsize_lst(ite) = stepsize;

     %  -----------------

     unewAll =  uoldAll;   % boundary values
     unewAll(isfreeNode) =  uoldAll(isfreeNode) + stepsize*du(isfreeNode);


    % Assemble Mass Matrix
    if ~isempty(pde.Ku) && ~isnumeric(pde.Ku)       % d is a function
        gradu = reshape(grad*unewAll, NT, 2);
        gradunorm = sqrt(sum(gradu.^2, 2));
        K = pde.Ku(gradunorm ); % K is a scalar
    end
    Ap = B*([K; K]./[area;area].*Bt);
    Apeps = B*(([K; K]+keps)./[area;area].*Bt);

    % bdidx = zeros(N,1); 
    bdidx(~isfreeNode) = 1;
    Tbd = spdiags(bdidx,0,N,N);
    Te = spdiags(1-bdidx,0,N,N);
    bigAp = Te*Apeps*Te + Tbd;

    eqn.Ap = Ap;

    % compute error
    % Btu =  eqn.Bfft*unewAll(isfreeNode);
    resAll =  Ap*unewAll - eqn.f;
    res(ite) = norm(resAll(isfreeNode));
    energy(ite) =  sum(gradunorm.^p.*area)/p - f'*unewAll;
    disp(res(ite))
     
    if ite> 1 && res(ite) < tol*res0
        flag = 1;
        break;
    end

    uoldAll = unewAll;


%% 
end
time = toc;

% Assemble Mass Matrix
% if ~isempty(pde.Ku) && ~isnumeric(pde.Ku)       % d is a function
%     gradu = reshape(grad*uI, NT, 2);
%     gradunorm = sqrt(sum(gradu.^2, 2));
%     K = pde.Ku(gradunorm ); % K is a scalar
% end
% Ap = B*([K; K]./[area;area].*Bt);
% resI =  Ap*uI - eqn.f;
% disp(norm(resI(isfreeNode)));

% u = unewAll;
% u(isfreeNode) = Ap(isfreeNode, isfreeNode)\eqn.f(isfreeNode);

%% Output
info = struct('itStep', ite ,'res',res(1:ite), 'flag', flag, ...
    'Vcycle_num', Vcycle_num, 'time', time, 'stepsize_lst', stepsize_lst);
% plot(1:ite, log10(res(1:ite)));
u = unewAll;
sigma = K.*gradu;
end

