%% CONVERGENCE OF FAST SOLVERS FOR MIXED FINITE ELEMENT METHOD (P1-P0) FOR p-Laplacian EQUATIONS
%
% This example is to show the convergence of fast solvers of mixed finite
% element (P1-P0) approximation of the p-Laplacian equation.
%er

close all
clear variables
%% Setting
h = 1/4;
% [node,elem] = squaremesh([0,1,0,1],h); 
[node,elem] = circlemesh(0,0,1,h);
% [node,elem] = squaremesh([-1,1,-1,1],h);
mesh = struct('node',node,'elem',elem);
option.mg.solvermaxIt = 5;
option.mg.solver = 'Vcycle';
option.mg.tol = 1e-2;
option.mg.printlevel = 0;
option.rateflag = 1;
option = mfemoption(option);

L0 = 4;
maxIt = 5;
maxN = 1e7;

% pde = DFdata2;
pde = pLapdata1;
bdFlag = setboundary(node,elem,'Dirichlet');
% bdFlag = setboundary(node,elem,'Neumann');
% mfemDarcy(mesh,pde,option);

%% Generate an initial mesh 
for k = 1:L0
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    [~,~,isBdNode] = findboundary(elem,bdFlag);
    node(isBdNode, :) =  node(isBdNode, :)./(sqrt(sum(node(isBdNode, :).^2, 2)));
    % showmesh(node,elem);
%     if isfield(pde,'K') && isnumeric(pde.K)
%         pde.K = repmat(pde.K,4,1); % prolongate piecwise K to the fine grid
%     end    
end

%% Initialize err and time array
% error
errsigmaL2 = zeros(maxIt,1);   
erruIuhL2 = zeros(maxIt,1); 
erruL2 = zeros(maxIt,1);   
uL2 = zeros(maxIt,1);   
relerruL2 = zeros(maxIt,1);   
% erruHdiv = zeros(maxIt,1); 
errsigmaIsigmah = zeros(maxIt,1); 

% info
N = zeros(maxIt,1); 
h_lst = zeros(maxIt,1);
itStep = zeros(maxIt,1);
time = zeros(maxIt,1);
flag = zeros(maxIt,1);
Vcycle_num = zeros(maxIt,1);

%% Finite Element Method        
for k = 1:maxIt
    
    tic
    % solve the equation
  [sigma,u,eqn,info] =  pLapP0P1_DualTPDJ(node,elem,bdFlag,pde,option);
  % [sigma,u,eqn,info] =  pLapP0P1_DualTPDM(node,elem,bdFlag,pde,option);
     % [sigma, u,eqn,info] =  pLapP0P1_primal(node,elem,bdFlag,pde,option);
     % [sigma, u,eqn,info] =  pLapP0P1_primal_nolinesearch(node,elem,bdFlag,pde,option);
      % [sigma, u,eqn,info] =  pLapP0P1_primalsmallp(node,elem,bdFlag,pde,option);
        % [sigma, u,eqn,info] =  pLapP0P1_primalsmallp_nolinesearch(node,elem,bdFlag,pde,option);
   time(k) = toc;
    
    % compute error
    if isfield(pde,'exactsigma')
        errsigmaL2(k) = getL2error(node,elem,pde.exactsigma,sigma);
%         erruHdiv(k) = getHdiverrorRT0(node,elem,pde.exactudiv,u,[]);
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

    % plot 
    N(k) = length(sigma) + length(u);
    h_lst(k) = 2./(sqrt(size(node,1))-1);

    if N(k) > maxN
        break;
    end
   
    % refine mesh
    [node,elem, bdFlag] = uniformrefine(node,elem,bdFlag);
    [~,~,isBdNode] = findboundary(elem,bdFlag);
    node(isBdNode, :) =  node(isBdNode, :)./(sqrt(sum(node(isBdNode, :).^2, 2)));

    
end
%% Plot convergence rates
% if option.rateflag
%     figure;
%     set(gcf,'Units','normal'); 
%     set(gcf,'Position',[0.25,0.25,0.55,0.4]);
%     subplot(1,2,1)
%     showrateh2(h(1:k),errpIphL2(1:k),1,'-*','||p_I-p_h||_{\infty}',...
%                h(1:k),errpL2(1:k),1,'k-+','||p-p_h||');
%     subplot(1,2,2)
%     showrateh2(h(1:k),erruL2(1:k),1,'k-*','|| u - u_h||',...
%                h(1:k),erruIuh(1:k),1,'m-+','|| u_I - u_h||');
% %                h(1:k),erruHdiv(1:k),1,'-+','||div(u - u_h)||',... 
% end

%% Output
err = struct('h',h(1:k),'N',N(1:k),'uL2',erruL2(1:k),'uIuhL2',erruIuhL2(1:k),...
             'sigmaL2',errsigmaL2(1:k), 'sigmaIsigmaL2',errsigmaIsigmahL2(1:k));
            
%% Display error and CPU time
% disp('Table: Error')
% colname = {'#Dof','h','||p-p_h||','||p_I-p_h||','||u-u_h||'};
% disptable(colname,err.N,[],err.h,'%0.2e',err.pL2,'%0.5e',err.pIphL2,'%0.5e',...
%                      err.uL2,'%0.5e');

%% plot time growth
% N = [ 12417  49409 197121 787457 3147777];
% Time_DualTPDJ_p1d5 = [0.043 0.17 0.46 1.6 7.3];
% Time_DualTPDJ_p4 = [0.18 0.6 1.8 7.1 34 ];
% Time_DualTPDM_p1d5 = [0.11 0.38 1.1 3.9 16];
% Time_DualTPDM_p4 = [0.17 0.62 1.6 5.6 24];

% N = [ 49409 197121 787457 3147777];
% Time_DualTPDJ_p1d5 = [0.17 0.46 1.6 7.3];
% Time_DualTPDJ_p4 = [ 0.6 1.8 7.1 34 ];
% Time_DualTPDM_p1d5 = [0.38 1.1 3.9 16];
% Time_DualTPDM_p4 = [0.62 1.6 5.6 24];
% 
% r_DualTPDJ_p1d5 = showrate(N,Time_DualTPDJ_p1d5,1,'-*');
% hold on
% r_DualTPDM_p1d5 = showrate(N,Time_DualTPDM_p1d5,1,'k-x');
% r_DualTPDJ_p4 = showrate(N,Time_DualTPDJ_p4,1,'k-x');
% r_DualTPDM_p4 = showrate(N,Time_DualTPDM_p4,1,'g-x');
% 
% % title(['Rate of time growth'],'FontSize', 154);
% h_legend = legend('DualTPD-J, p = 1.5',['N^{' num2str(r_DualTPDJ_p1d5,2) '}'],...
%                      'DualTPD-M, p = 1.5',['N^{' num2str(r_DualTPDM_p1d5,2) '}'],...
%                   'DualTPD-J, p = 4',['N^{' num2str(r_DualTPDJ_p4,2) '}'],...
%                   'DualTPD-M, p = 4',['N^{' num2str(r_DualTPDM_p4,2) '}'],...
%                   'LOCATION','Best');
% set(h_legend,'FontSize',16);
% xlabel('N','FontSize',16);
% ylabel('CPU time','FontSize',16);
% 
% 
% 
% %% Plot convergence rates
% figure(1);
% h_lst = h;
% showrateh2(h_lst,erruL2,1,'k-+','|| u-u_h ||',...
%            h_lst,erruIuhL2,1,'r-+','|| u_I-u_h ||');
% 
% 
% figure(2);
% showrateh2(h_lst,errsigmaL2,1,'k-+','|| \sigma-\sigma_h ||',...
%            h_lst,errsigmaIsigmah,1,'r-+','||\sigma_I-\sigma_h||');
% 
