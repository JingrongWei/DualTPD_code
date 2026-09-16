function [sigma,u,eqn,info] = pLapP0P1_Newton(node,elem,bdFlag,pde,option)
% Newton solver for the P1-P0 discretization of the p-Laplacian equation.

if ~exist('option','var'), option = []; end
if ~exist('bdFlag','var'), bdFlag = []; end

N = size(node,1);
maxIt = 4e3;
tol = 1e-6;
flag = 0;
jacobianEps = 1e-4;
stepsize = 1;

%% Data structure
[~,~,isBdNode] = findboundary(elem,bdFlag);
isfreeNode = ~isBdNode;

NT = size(elem,1);
Nu = N;

% Exact solution is used only to impose the Dirichlet boundary data.
uI = pde.exactu(node);

% rng(4);
% uoldAll = rand(Nu,1);
% 
uoldAll = zeros(Nu,1);

uoldAll(~isfreeNode) = uI(~isfreeNode);


%% Assemble the discrete gradient and load vector
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

% [Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleNewtonSystem(uoldAll);
% [Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleJacobianFromJsigma(uoldAll);
[Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleJacobianFromJu(uoldAll);
eqn = struct('Ap',Ap,'Jacobian',Jacobian,'B',B,'f',b,...
    'Bff',B(isfreeNode,:),'Bfft',B(isfreeNode,:)','ff',b(isfreeNode));

%% Newton iteration
res = zeros(maxIt,1);
res0 = norm(b(isfreeNode));
if res0 == 0, res0 = 1; end
Vcycle_num = 0;


tic
for ite = 1:maxIt
    residual = Ap*uoldAll - b;
    rhs = zeros(Nu,1);
    rhs(isfreeNode) = residual(isfreeNode);

    % Solve J(u_k) delta = F(u_k) with the multigrid solver.
    [delta,mginfo] = mg(bigJacobian,rhs,elem,option.mg);
    Vcycle_num = Vcycle_num + mginfo.itStep;

    unewAll = uoldAll;
    unewAll(isfreeNode) = uoldAll(isfreeNode) - stepsize*delta(isfreeNode);
    % [Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleNewtonSystem(unewAll);
    % [Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleJacobianFromJsigma(unewAll);
    [Ap,Jacobian,bigJacobian,gradu,gradunorm,K] = assembleJacobianFromJu(unewAll);
    residual = Ap*unewAll - b;

    eqn.Ap = Ap;
    eqn.Jacobian = Jacobian;

    res(ite) = norm(residual(isfreeNode));
    disp(res(ite))

    if ite > 1 && res(ite) < tol*res0
        flag = 1;
        break;
    end

    uoldAll = unewAll;
end
time = toc;

info = struct('itStep',ite,'res',res(1:ite),'flag',flag,...
    'Vcycle_num',Vcycle_num,'time',time,'stepsize',stepsize);
u = unewAll;
sigma = K.*gradu;

    function [A,J,bigJ,graduLocal,gradunormLocal,KLocal] = assembleJacobianFromJu(uLocal)
        graduLocal = reshape(grad*uLocal,NT,2);
        gradunormLocal = sqrt(sum(graduLocal.^2,2));
        KLocal = pde.Ku(gradunormLocal);

        Jlocal = pde.Ju(graduLocal);
        A = B*spdiags([KLocal; KLocal]./[area; area],0,2*NT,2*NT)*Bt;
        J = B*Jlocal*grad;
        bigJ = Te*J*Te + Tbd;
    end

    function [A,J,bigJ,graduLocal,gradunormLocal,KLocal] = assembleJacobianFromJsigma(uLocal)
        graduLocal = reshape(grad*uLocal,NT,2);
        gradunormLocal = sqrt(sum(graduLocal.^2,2));
        KLocal = pde.Ku(gradunormLocal);
        sigmaLocal = KLocal.*graduLocal;

        % Jsigma is D_sigma(grad u); its inverse is the primal Newton
        % coefficient D_(grad u)(sigma).
        Jsigma = pde.Jsigma(sigmaLocal);
        Jlocal = Jsigma\speye(2*NT);

        A = B*spdiags([KLocal; KLocal]./[area; area],0,2*NT,2*NT)*Bt;
        J = B*Jlocal*grad;
        bigJ = Te*J*Te + Tbd;
    end

    function [A,J,bigJ,graduLocal,gradunormLocal,KLocal] = assembleNewtonSystem(uLocal)
        graduLocal = reshape(grad*uLocal,NT,2);
        gradunormLocal = sqrt(sum(graduLocal.^2,2));
        KLocal = pde.Ku(gradunormLocal);

        % A zero gradient makes the exact tangent singular.  The regularized
        % tangent is used only in the Newton solve; the residual remains Ap*u-f.
        tangentNorm = gradunormLocal;
        tangentNorm(tangentNorm < eps) = jacobianEps;
        if pde.coff_p == 2
            tangentFactor = zeros(NT,1);
            tangentK = ones(NT,1);
        else
            tangentFactor = (pde.coff_p - 2)*pde.Jdu(tangentNorm);
            tangentK = tangentNorm.^(pde.coff_p - 2);
        end

        gx = graduLocal(:,1);
        gy = graduLocal(:,2);
        j11 = tangentK + tangentFactor.*gx.^2;
        j22 = tangentK + tangentFactor.*gy.^2;
        j12 = tangentFactor.*gx.*gy;
        Jlocal = sparse([1:NT, NT+1:2*NT, 1:NT, NT+1:2*NT],...
            [1:NT, NT+1:2*NT, NT+1:2*NT, 1:NT],...
            [j11; j22; j12; j12],2*NT,2*NT);

        A = B*spdiags([KLocal; KLocal]./[area; area],0,2*NT,2*NT)*Bt;
        J = B*Jlocal*grad;
        bigJ = Te*J*Te + Tbd;
    end
end
