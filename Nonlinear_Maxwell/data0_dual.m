function pde = data0_dual
%% data for Maxwell equation
% exponential coefficient
% u = [0, 0, cos(omega*x).*cos(omega*y).*cos(omega*z)]


% parameters u7ed o define reluctivity
a0 = 10; %2.70675; %2.721 6.775 
a1 = 73.89; %20;  % 50
a2 = 1;  %10;  % 20
omega = 4;


pde = struct('nu',@nu,'Dnu',@Dnu, 'Phi_inv', @Phi_inv,'Jsigma', @Jsigma, ...
    'Jsigmainv', @Jsigmainv, 'J', @J, 'Ju', @Ju, 'exactu', @exactu,'g_D', @g_D,...
    'curlu', @curlu, 'g', @g,'exactsigma', @exactsigma,'exactp', @exactp, 'exactsigmax', @exactsigmax,'a1', a1);

    % coefficient
    function z = nu(s)
        z = a0 + a1*exp(-a2*s);
        % z = s.^3;
    end

    function z = Dnu(s)
        z = -a1*(a2)*exp(-a2*s);
        % z = 3*s.^2;
    end

    function z = Phi_inv(s)
        f = @(s) nu(s).*s;
        df = @(s) nu(s) + Dnu(s).*s;
        [res, z] = numerical_inverse(f, df, s, 20, 5);
        fprintf('numeric inv residual is %g \n',max(abs(res)));
   end

    % right hand side function
    function z = J(p)
        x = p(:,1); y = p(:,2); z=p(:, 3);
        l = ((cos(omega*x).*sin(omega*y).*cos(omega*z)).^2 ...
            + (sin(omega*x).*cos(omega*y).*cos(omega*z)).^2)*omega^2; % norm of curl(u)
        curlup = curlu(p);
        f1p = (-cos(omega*x).*sin(omega*y).*cos(omega*z)).*(sin(omega*x).*sin(omega*y).*cos(omega*z)) + ...
            (sin(omega*x).*cos(omega*y).*cos(omega*z)).*(cos(omega*x).*cos(omega*y).*cos(omega*z));
        f2p = (-cos(omega*x).*sin(omega*y).*cos(omega*z)).*(-cos(omega*x).*cos(omega*y).*cos(omega*z)) + ...
            (sin(omega*x).*cos(omega*y).*cos(omega*z)).*(-sin(omega*x).*sin(omega*y).*cos(omega*z));
        f3p = (-cos(omega*x).*sin(omega*y).*cos(omega*z)).*(cos(omega*x).*sin(omega*y).*sin(omega*z)) + ...
            (sin(omega*x).*cos(omega*y).*cos(omega*z)).*(-sin(omega*x).*cos(omega*y).*sin(omega*z));
        J1 = [sin(omega*x).*cos(omega*y).*sin(omega*z), ...
              cos(omega*x).*sin(omega*y).*sin(omega*z), ...
              2*cos(omega*x).*cos(omega*y).*cos(omega*z)]*omega^2;
        z = nu(l.^(1/2)).*J1 + ...
            a1*(-a2)*exp(-a2*l.^(1/2)).*(l.^(-1/2)).*cross([f1p,f2p,f3p],curlup)*omega^3;
    end

    function M = Jsigma(sigma, K, z)
        % K = 1/nu(z), z = Phi_inv(sigma)
        Nsigma = length(sigma);
        NT = Nsigma/3;
        sigmamat = reshape(sigma, NT, 3);
        sigmal2norm = sum(sigmamat.^2, 2).^(1/2);
        s = -Dnu(z)./((Dnu(z).*z + K).*(K.^2).*sigmal2norm);
        ii = [1:Nsigma, 1:Nsigma, 1:NT, 1:NT, NT+1:2*NT, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma];
        jj = [1:Nsigma, 1:Nsigma, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma, 1:NT, 1:NT, NT+1:2*NT];
        ss = [1./K; 1./K; 1./K; [s; s; s].*sigma.^2; ...
            s.*sigma(1:NT).*sigma(NT+1:2*NT); s.*sigma(1:NT).*sigma(2*NT+1:Nsigma);...
            s.*sigma(NT+1:2*NT).*sigma(2*NT+1:Nsigma); s.*sigma(NT+1:2*NT).*sigma(1:NT);...
            s.*sigma(2*NT+1:Nsigma).*sigma(1:NT); s.*sigma(2*NT+1:Nsigma).*sigma(NT+1:2*NT)];
        M = sparse(ii, jj, ss, Nsigma, Nsigma);
%         M = sparse(1:Nu, 1:Nu, [s; s], Nu, Nu);
%         M = M + sparse(1:Nu, 1:Nu, beta*u.^2./[unorm; unorm], Nu, Nu);
%         M = M + sparse(1:NT, NT+1:Nu,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
%         M = M + sparse( NT+1:Nu, 1:NT,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
    end

    function M = Ju(curlu, curlunorm)
        % K = nu(|curl u|)
        NT = length(curlu)/3;
        nuu = nu(curlunorm);
        Dnuu = Dnu(curlunorm);
        s = Dnuu./curlunorm;
        ii = [1:3*NT, 1:3*NT, 1:NT, 1:NT, NT+1:2*NT, NT+1:2*NT, 2*NT+1:3*NT, 2*NT+1:3*NT];
        jj = [1:3*NT, 1:3*NT, NT+1:2*NT, 2*NT+1:3*NT, 2*NT+1:3*NT, 1:NT, 1:NT, NT+1:2*NT];
        ss = [nuu; nuu; nuu; [s; s; s].*curlu.^2; ...
            s.*curlu(1:NT).*curlu(NT+1:2*NT); s.*curlu(1:NT).*curlu(2*NT+1:3*NT);...
            s.*curlu(NT+1:2*NT).*curlu(2*NT+1:3*NT); s.*curlu(NT+1:2*NT).*curlu(1:NT);...
            s.*curlu(2*NT+1:3*NT).*curlu(1:NT); s.*curlu(2*NT+1:3*NT).*curlu(NT+1:2*NT)];
        M = sparse(ii, jj, ss, 3*NT, 3*NT);
    %         M = sparse(1:Nu, 1:Nu, [s; s], Nu, Nu);
    %         M = M + sparse(1:Nu, 1:Nu, beta*u.^2./[unorm; unorm], Nu, Nu);
    %         M = M + sparse(1:NT, NT+1:Nu,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
    %         M = M + sparse( NT+1:Nu, 1:NT,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
    end

    function M = Jsigmainv(sigma, K, z)
        % K = 1/nu(z), z = Phi_inv(sigma)
        Nsigma = length(sigma);
        NT = Nsigma/3;
        sigmamat = reshape(sigma, NT, 3);
        sigmal2norm = sum(sigmamat.^2, 2).^(1/2);
        s = Dnu(z)./((Dnu(z).*z + K).*(K.^2).*sigmal2norm);
        sinv = s.*K.^2./(-1 + s.*K.*(sigmal2norm.^2));
        ii = [1:Nsigma, 1:Nsigma, 1:NT, 1:NT, NT+1:2*NT, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma];
        jj = [1:Nsigma, 1:Nsigma, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma, 1:NT, 1:NT, NT+1:2*NT];
        ss = [K; K; K; -[sinv; sinv; sinv].*sigma.^2; ...
            -sinv.*sigma(1:NT).*sigma(NT+1:2*NT); -sinv.*sigma(1:NT).*sigma(2*NT+1:Nsigma);...
            -sinv.*sigma(NT+1:2*NT).*sigma(2*NT+1:Nsigma); -sinv.*sigma(NT+1:2*NT).*sigma(1:NT);...
            -sinv.*sigma(2*NT+1:Nsigma).*sigma(1:NT); -sinv.*sigma(2*NT+1:Nsigma).*sigma(NT+1:2*NT)];
        M = sparse(ii, jj, ss, Nsigma, Nsigma);
%         M = sparse(1:Nu, 1:Nu, [s; s], Nu, Nu);
%         M = M + sparse(1:Nu, 1:Nu, beta*u.^2./[unorm; unorm], Nu, Nu);
%         M = M + sparse(1:NT, NT+1:Nu,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
%         M = M + sparse( NT+1:Nu, 1:NT,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
    end


    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2); z=p(:, 3);
        u = [0*x, 0*y, cos(omega*x).*cos(omega*y).*cos(omega*z)];
    end

    function u =  exactp(p)
        x = p(:,1);
        u = 0*x;
    end

    function sigma = exactsigma(p)
         x = p(:,1); y = p(:,2); z=p(:, 3);
         curlup =  curlu(p);
         l = ((cos(omega*x).*sin(omega*y).*cos(omega*z)).^2 ...
            + (sin(omega*x).*cos(omega*y).*cos(omega*z)).^2)*omega^2; % norm of curl(u)
         sigma = nu(l.^(1/2)).*curlup;
    end

    function sigma = exactsigmax(p)
         x = p(:,1); y = p(:,2); z=p(:, 3);
         l = ((cos(omega*x).*sin(omega*y).*cos(omega*z)).^2 ...
            + (sin(omega*x).*cos(omega*y).*cos(omega*z)).^2)*omega^2; % norm of curl(u)
         sigma = nu(l.^(1/2)).*-cos(omega*x).*sin(omega*y).*cos(omega*z)*omega;
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

   % curl of exact solution
   function u =  curlu(p)
       x = p(:,1); y = p(:,2); z=p(:, 3);
       u = [-cos(omega*x).*sin(omega*y).*cos(omega*z), ...
            sin(omega*x).*cos(omega*y).*cos(omega*z), ...
            0*z]*omega;
   end

    % divergence rhs
    function z =  g(p)
        x = p(:,1); y = p(:,2); z=p(:,3);
        z = -cos(omega*x).*cos(omega*y).*sin(omega*z)*omega;
    end
end