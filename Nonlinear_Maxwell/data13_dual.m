function pde = data13_dual
%% data for Maxwell equation
% p-curl example
% u = [0, 0, (x+1)(1-x) + (y + 1)(1- y) + (z + 1)(1- z)]


% parameters used to define reluctivity
coff_p = 1.5;
nueps = 0.0001;

pde = struct('nu',@nu,'Dnu',@Dnu, 'Ksigma', @Ksigma, 'Phi_inv', @Phi_inv,'Jsigma', @Jsigma, ...
    'Jsigmainv', @Jsigmainv, 'J', @J, 'Ju', @Ju, 'exactu', @exactu,'g_D', @g_D,...
    'curlu', @curlu, 'g', @g,'exactsigma', @exactsigma,'exactp', @exactp, 'exactsigmax', @exactsigmax);

    % coefficient
    function z = nu(s)
        if coff_p < 2
            s = s + (s < eps)*nueps;
        end

        z = s.^(coff_p-2);

        if coff_p > 2
            z = z + (z < eps)*nueps;
        end
        % z = s.^3;
    end

    function z = Dnu(s)
        z = (coff_p-2)*s.^(coff_p-3);
        % z = 3*s.^2;
    end

    function z = Phi_inv(s)
        z = s.^(1/(coff_p -1));
   end

    % right hand side function
    function z = J(p)
        x = p(:,1); y = p(:,2); z=p(:, 3);
        l = (-y).^2 + x.^2; % norm of curl(u)
        curlup = curlu(p);
        f1p = x;
        f2p = y;
        f3p = 0*z;
        J1 = [0*x, 0*y, 2+0*z];
        if coff_p >= 4
           z = nu(l.^(1/2)).*J1 + (coff_p-2)*l.^((coff_p - 4)/2).*cross([f1p,f2p,f3p],curlup);
        else
           z = nu(l.^(1/2)).*J1 + (coff_p-2)*(l + (l < nueps).*nueps).^((coff_p - 4)/2).*cross([f1p,f2p,f3p],curlup); 
        end
        
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

%     function M = Jsigmainv(sigma, K, z)
%         % K = 1/nu(z), z = Phi_inv(sigma)
%         Nsigma = length(sigma);
%         NT = Nsigma/3;
%         sigmamat = reshape(sigma, NT, 3);
%         sigmal2norm = sum(sigmamat.^2, 2).^(1/2);
%         s = Dnu(z)./((Dnu(z).*z + K).*(K.^2).*sigmal2norm);
%         sinv = s.*K.^2./(-1 + s.*K.*(sigmal2norm.^2));
%         ii = [1:Nsigma, 1:Nsigma, 1:NT, 1:NT, NT+1:2*NT, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma];
%         jj = [1:Nsigma, 1:Nsigma, NT+1:2*NT, 2*NT+1:Nsigma, 2*NT+1:Nsigma, 1:NT, 1:NT, NT+1:2*NT];
%         ss = [K; K; K; -[sinv; sinv; sinv].*sigma.^2; ...
%             -sinv.*sigma(1:NT).*sigma(NT+1:2*NT); -sinv.*sigma(1:NT).*sigma(2*NT+1:Nsigma);...
%             -sinv.*sigma(NT+1:2*NT).*sigma(2*NT+1:Nsigma); -sinv.*sigma(NT+1:2*NT).*sigma(1:NT);...
%             -sinv.*sigma(2*NT+1:Nsigma).*sigma(1:NT); -sinv.*sigma(2*NT+1:Nsigma).*sigma(NT+1:2*NT)];
%         M = sparse(ii, jj, ss, Nsigma, Nsigma);
% %         M = sparse(1:Nu, 1:Nu, [s; s], Nu, Nu);
% %         M = M + sparse(1:Nu, 1:Nu, beta*u.^2./[unorm; unorm], Nu, Nu);
% %         M = M + sparse(1:NT, NT+1:Nu,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
% %         M = M + sparse( NT+1:Nu, 1:NT,beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
%     end

    function M = Jsigmainv(sigma, teps, eps0)
         % sigma is dim*2 matrix
        sigma = reshape(sigma, length(sigma)/3, 3);
        Nsigma = size(sigma, 1);
        sigmanorm = sqrt(sum(sigma.^2, 2));
        if coff_p > 2 
            zerosnode = (sigmanorm < eps0);
            index = (1:Nsigma);
            index = index(~zerosnode);
            s = sigmanorm(index).^((coff_p-2)/(coff_p-1));
            % s = 1./(sigmanorm.^(-(coff_p-2)/(coff_p-1)));
            % t0 = -(coff_p -2)/(coff_p -1)*(sigmanorm(index).^(-(3*coff_p-4)/(coff_p-1)));
            % t = -s./(1./(t0.*s) + sigmanorm(index).^2);
            t = (coff_p - 2)*sigmanorm(index).^(coff_p/(1- coff_p));
            s1 = zeros(Nsigma,1);
            s1(index) = s;           
            s1(zerosnode) = teps^((coff_p-2)/(coff_p-1));
            ii = [1:3*Nsigma, index, index+Nsigma, index+2*Nsigma,...
                index, index+Nsigma, index, index+2*Nsigma, index+Nsigma,index+2*Nsigma];
            jj = [1:3*Nsigma, index, index+Nsigma, index+2*Nsigma, ...
                index + Nsigma, index, index+2*Nsigma, index, index+2*Nsigma, index+Nsigma];
            ss = [s1; s1; s1; t.*sigma(index, 1).^2; t.*sigma(index, 2).^2; t.*sigma(index, 3).^2;...
                  t.*sigma(index, 1).*sigma(index, 2); t.*sigma(index, 1).*sigma(index, 2);...
                 t.*sigma(index, 1).*sigma(index, 3);t.*sigma(index, 1).*sigma(index, 3);
                 t.*sigma(index, 3).*sigma(index, 2); t.*sigma(index, 3).*sigma(index, 2);];
            M = sparse(ii, jj, ss, 3*Nsigma, 3*Nsigma);
        elseif coff_p == 2
            M = sparse(1:3*Nsigma, 1:3*Nsigma, ones(3*Nsigma,1), 3*Nsigma, 3*Nsigma);
        else
            zerosnode = (sigmanorm.^((2-coff_p)/(coff_p-1)) < eps0);
            index = (1:Nsigma);
            index = index(~zerosnode);
            s = sigmanorm(index).^((coff_p-2)/(coff_p-1));
            % s = 1./(sigmanorm.^(-(coff_p-2)/(coff_p-1)));
            % t0 = -(coff_p -2)/(coff_p -1)*(sigmanorm(index).^(-(3*coff_p-4)/(coff_p-1)));
            % t = -s./(1./(t0.*s) + sigmanorm(index).^2);
            t = (coff_p - 2)*sigmanorm(index).^(coff_p/(1- coff_p));
            s1 = zeros(Nsigma,1);
            s1(index) = s;           
            s1(zerosnode) = 1./(sigmanorm(zerosnode).^((2-coff_p)/(coff_p-1)) + teps);
            ii = [1:3*Nsigma, index, index+Nsigma, index+2*Nsigma,...
                index, index+Nsigma, index, index+2*Nsigma, index+Nsigma,index+2*Nsigma];
            jj = [1:3*Nsigma, index, index+Nsigma, index+2*Nsigma, ...
                index + Nsigma, index, index+2*Nsigma, index, index+2*Nsigma, index+Nsigma];
            ss = [s1; s1; s1; t.*sigma(index, 1).^2; t.*sigma(index, 2).^2; t.*sigma(index, 3).^2;...
                  t.*sigma(index, 1).*sigma(index, 2); t.*sigma(index, 1).*sigma(index, 2);...
                 t.*sigma(index, 1).*sigma(index, 3);t.*sigma(index, 1).*sigma(index, 3);
                 t.*sigma(index, 3).*sigma(index, 2); t.*sigma(index, 3).*sigma(index, 2);];
            M = sparse(ii, jj, ss, 3*Nsigma, 3*Nsigma);

            % s = 1./(sigmanorm.^(-(coff_p-2)/(coff_p-1))+ teps);
            % if coff_p < 4/3
            %     t0 = -(coff_p - 2)/(coff_p -1)*(sigmanorm.^(-(3*coff_p-4)/(coff_p-1)));
            %     t = -(t0.*s.*s)./(1 + t0.*s.*sigmanorm.^2);
            %     % t = (coff_p - 2)*sigmanorm.^(coff_p/(1- coff_p));
            %     ii = [1:2*Nsigma, 1:2*Nsigma, 1:Nsigma, Nsigma+1:2*Nsigma];
            %     jj = [1:2*Nsigma, 1:2*Nsigma, Nsigma+1:2*Nsigma, 1:Nsigma];
            %     ss = [s; s; t.*sigma(:, 1).^2; t.*sigma(:, 2).^2; t.*sigma(:, 1).*sigma(:, 2); ...
            %          t.*sigma(:, 1).*sigma(:, 2)];
            %     M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);
            % else
            %     t0 = -(coff_p - 2)/(coff_p -1)*((sigmanorm+teps).^(-(3*coff_p-4)/(coff_p-1)));
            %     t = -(t0.*s.*s)./(1 + t0.*s.*sigmanorm.^2);
            %     ii = [1:2*Nsigma, 1:2*Nsigma, 1:Nsigma, Nsigma+1:2*Nsigma];
            %     jj = [1:2*Nsigma, 1:2*Nsigma, Nsigma+1:2*Nsigma, 1:Nsigma];
            %     ss = [s; s; t.*sigma(:, 1).^2; t.*sigma(:, 2).^2; t.*sigma(:, 1).*sigma(:, 2); ...
            %          t.*sigma(:, 1).*sigma(:, 2)];
            %     M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);
            % end
        end
       
        % detu = (beta*(u(1:NT).^2)./unorm+s).*(beta*(u(NT+1:end).^2)./unorm+s) ...
        %     - (beta^2)*(u(1:NT).^2).*(u(NT+1:end).^2)./(unorm.^2);
        % detuinv = 1./detu;
        % M = sparse(1:Nu, 1:Nu, [s + beta*(u(NT+1:end).^2./unorm) ; beta*(u(1:NT).^2)./unorm+s], Nu, Nu);
        % M = M + sparse(1:NT, NT+1:Nu,-beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
        % M = M + sparse( NT+1:Nu, 1:NT,-beta*u(1:NT).*u(NT+1:end)./unorm, Nu, Nu);
        % M = [ detuinv; detuinv].*M;
    end


 

    % exact solution
    function u =  exactu(p)
        x = p(:,1); y = p(:,2); z=p(:, 3);
        u = [0*x, 0*y, (x+1).*(1-x)/2 + (y+1).*(1-y)/2 + (z + 1).*(1- z)/2 ];
    end

    function u =  exactp(p)
        x = p(:,1);
        u = 0*x;
    end


    function sigma = exactsigma(p)
         x = p(:,1); y = p(:,2); z=p(:, 3);
         curlup =  curlu(p);
         l = (-y).^2 + x.^2; % norm of curl(u)
         sigma = nu(l.^(1/2)).*curlup;
    end

    function sigma = exactsigmax(p)
         x = p(:,1); y = p(:,2); z = p(:, 3);
         l = (-y).^2 + x.^2; % norm of curl(u)
         sigma = nu(l.^(1/2)).*(-y);
    end

    % Dirichlet boundary condition
    function u =  g_D(p)
        u =  exactu(p);
    end

   % curl of exact solution
   function u =  curlu(p)
       x = p(:,1); y = p(:,2); z=p(:, 3);
       u = [-y, x, 0*z];
   end

    % divergence rhs
    function z =  g(p)
        x = p(:,1); y = p(:,2); z=p(:,3);
        z = -z;
    end
end