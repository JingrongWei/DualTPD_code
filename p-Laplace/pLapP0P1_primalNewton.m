function [sigma,u,eqn,info] = pLapP0P1_primalNewton(node,elem,bdFlag,pde,option)
% Hybrid PGD-Newton solver for the P1-P0 discretization of p-Laplacian.

if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1);
maxIt = 3e4;
tol = 1e-6;
Phase_change_TOL = 1e-2;
pgdStepsize = 0.2;
newtonStepsize = 1;
keps = 1e-4;
flag = 0;

%% Data structure and initial value
[~,~,isBdNode] = findboundary(elem,bdFlag);
isfreeNode = ~isBdNode;
NT = size(elem,1);
Nu = N;

uI = pde.exactu(node);

rng(4);
uoldAll = rand(Nu,1);

% uoldAll = zeros(Nu,1);
uoldAll(~isfreeNode) = uI(~isfreeNode);

%% Assemble discrete gradient and load vector
[Dlambda,area,~] = gradbasis(node,elem);
d1 = Dlambda(:,:,1).*[area, area];
d2 = Dlambda(:,:,2).*[area, area];
d3 = Dlambda(:,:,3).*[area, area];
Dx = sparse(elem(:),repmat((1:NT)',3,1),...
    [d1(:,1); d2(:,1); d3(:,1)],Nu,NT);
Dy = sparse(elem(:),repmat((1:NT)',3,1),...
    [d1(:,2); d2(:,2); d3(:,2)],Nu,NT);
B = [Dx Dy];
Bt = B';
grad = spdiags(1./[area; area],0,2*NT,2*NT)*Bt;

b = zeros(Nu,1);
if ~isfield(option,'fquadorder'), option.fquadorder = 3; end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f == 0))
    pde.f = [];
end
if isreal(pde.f)
    switch length(pde.f)
        case NT
            bt = pde.f.*area/3;
            b = accumarray(elem(:),[bt; bt; bt],[Nu 1]);
        case Nu
            bt = zeros(NT,3);
            bt(:,1) = area.*(2*pde.f(elem(:,1)) + pde.f(elem(:,2)) + pde.f(elem(:,3)))/12;
            bt(:,2) = area.*(2*pde.f(elem(:,2)) + pde.f(elem(:,3)) + pde.f(elem(:,1)))/12;
            bt(:,3) = area.*(2*pde.f(elem(:,3)) + pde.f(elem(:,1)) + pde.f(elem(:,2)))/12;
            b = accumarray(elem(:),bt(:),[Nu 1]);
        case 1
            bt = pde.f*area/3;
            b = accumarray(elem(:),[bt; bt; bt],[Nu 1]);
    end
end
if ~isempty(pde.f) && ~isreal(pde.f)
    [lambda,weight] = quadpts(option.fquadorder);
    bt = zeros(NT,3);
    for q = 1:size(lambda,1)
        pxy = lambda(q,1)*node(elem(:,1),:) ...
            + lambda(q,2)*node(elem(:,2),:) ...
            + lambda(q,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for j = 1:3
            bt(:,j) = bt(:,j) + weight(q)*lambda(q,j)*fp;
        end
    end
    b = accumarray(elem(:),bt(:).*repmat(area,3,1),[Nu 1]);
end

bdidx = zeros(Nu,1);
bdidx(~isfreeNode) = 1;
Tbd = spdiags(bdidx,0,Nu,Nu);
Te = spdiags(1-bdidx,0,Nu,Nu);

[Ap,bigAp,gradu,K] = assemblePrimalOperators(uoldAll);
Jacobian = [];
bigJacobian = [];
eqn = struct('Ap',Ap,'Jacobian',Jacobian,'B',B,'f',b,...
    'Bff',B(isfreeNode,:),'Bfft',B(isfreeNode,:)','ff',b(isfreeNode));

%% PGD phase followed by Newton phase
res = zeros(maxIt,1);
stepsize_lst = zeros(maxIt,1);
res0 = norm(b(isfreeNode));
if res0 == 0, res0 = 1; end
Vcycle_num = 0;
phase = 'PGD';
phaseChangeIt = 0;
pgdIt = 0;
newtonIt = 0;

tic
for ite = 1:maxIt
    residual = Ap*uoldAll - b;
    rhs = zeros(Nu,1);
    rhs(isfreeNode) = residual(isfreeNode);

    if strcmp(phase,'PGD')
        [delta,mginfo] = mg(bigAp,rhs,elem,option.mg);
        stepsize = pgdStepsize;
        pgdIt = pgdIt + 1;
    else
        [delta,mginfo] = mg(bigJacobian,rhs,elem,option.mg);
        stepsize = newtonStepsize;
        newtonIt = newtonIt + 1;
    end
    Vcycle_num = Vcycle_num + mginfo.itStep;

    unewAll = uoldAll;
    unewAll(isfreeNode) = uoldAll(isfreeNode) - stepsize*delta(isfreeNode);
    [Ap,bigAp,gradu,K] = assemblePrimalOperators(unewAll);
    if strcmp(phase,'Newton')
        [Jacobian,bigJacobian] = assembleNewtonJacobian(gradu);
    end


    residual = Ap*unewAll - b;
    res(ite) = norm(residual(isfreeNode));
    stepsize_lst(ite) = stepsize;
    disp(res(ite))

    if ite > 1 && res(ite) < tol*res0
        flag = 1;
        break;
    end

    if strcmp(phase,'PGD') && res(ite)/res0 < Phase_change_TOL
        phase = 'Newton';
        phaseChangeIt = ite;
        [Jacobian,bigJacobian] = assembleNewtonJacobian(gradu);
    end

    eqn.Ap = Ap;
    eqn.Jacobian = Jacobian;

    uoldAll = unewAll;
end
time = toc;

info = struct('itStep',ite,'res',res(1:ite),'flag',flag,...
    'Vcycle_num',Vcycle_num,'time',time,'stepsize_lst',stepsize_lst(1:ite),...
    'Phase_change_TOL',Phase_change_TOL,'phaseChangeIt',phaseChangeIt,...
    'pgdIt',pgdIt,'newtonIt',newtonIt);
u = unewAll;
sigma = K.*gradu;

    function [A,bigA,graduLocal,KLocal] = assemblePrimalOperators(uLocal)
        graduLocal = reshape(grad*uLocal,NT,2);
        gradunorm = sqrt(sum(graduLocal.^2,2));
        KLocal = pde.Ku(gradunorm);

        A = B*spdiags([KLocal; KLocal]./[area; area],0,2*NT,2*NT)*Bt;
        Apeps = B*spdiags(([KLocal; KLocal] + keps)./[area; area],...
            0,2*NT,2*NT)*Bt;
        bigA = Te*Apeps*Te + Tbd;

    end

    function [J,bigJ] = assembleNewtonJacobian(graduLocal)
        Jlocal = pde.Ju(graduLocal);
        J = B*Jlocal*grad;
        bigJ = Te*J*Te + Tbd;
    end
end
