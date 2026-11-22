function pde = pLapdata1(coff_p)
%% Data for p-Laplacian equation 

%     u = a*x*(x-1)y(y-1);

coff_a = 10;
nueps = 0.0001;
pde = struct('f',@f, 'exactsigma',@exactsigma,'exactu',@exactu,'g_D',@g_D, 'K', @K,...
          'Ku', @Ku,'Jsigma', @Jsigma, 'Jsigmainv', @Jsigmainv, 'coff_p', coff_p, 'exactgradu', @exactgradu);

    % load data (right hand side function)
    function rhs =  f(pt)
        x = pt(:,1); y = pt(:,2);
        gradu = exactgradu(pt);
        gradunorm = sqrt(sum(gradu.^2, 2));
        rhs = -2*coff_a*(gradunorm.^(coff_p -2)).*(x.*(x-1)+y.*(y-1))...
            -(coff_p-2)*(gradunorm.^(coff_p -4))...
            .*(((coff_a*y.*(y-1).*(2*x-1)).*(2*coff_a*y.*(y-1)) ...
            + (coff_a*x.*(x-1).*(2*y-1)).*(coff_a*(2*x-1).*(2*y-1))).*(coff_a*y.*(y-1).*(2*x-1)) +...
             ((coff_a*y.*(y-1).*(2*x-1)).*(coff_a*(2*y-1).*(2*x-1))...
             + (coff_a*x.*(x-1).*(2*y-1)).*(2*coff_a*x.*(x-1))).*(coff_a*x.*(x-1).*(2*y-1)));
    end


    function s = K(sigmanorm) % nonlinear coefficient
        if coff_p  > 2
            sigmanorm = sigmanorm + (sigmanorm < eps)*nueps;
        end
        s = sigmanorm.^(-(coff_p - 2)/(coff_p-1));
    end

    function s = Ku(gradunorm) % nonlinear coefficient
        if coff_p  < 2
           gradunorm = gradunorm + (gradunorm < eps)*nueps;
        end
        s = gradunorm.^(coff_p - 2);
    end
    
    function M = Jsigma(sigma)
        % sigma is dim*2 matrix
        Nsigma = size(sigma, 1);
        sigmanorm = sqrt(sum(sigma.^2, 2));
        s = sigmanorm.^(-(coff_p-2)/(coff_p-1));
        t = -(coff_p -2)/(coff_p -1)*(sigmanorm.^(-(3*coff_p-4)/(coff_p-1)));
        ii = [1:2*Nsigma, 1:2*Nsigma, 1:Nsigma, Nsigma+1:2*Nsigma];
        jj = [1:2*Nsigma, 1:2*Nsigma, Nsigma+1:2*Nsigma, 1:Nsigma];
        ss = [s; s; t.*sigma(:, 1).^2; t.*sigma(:, 2).^2; t.*sigma(:, 1).*sigma(:, 2); ...
             t.*sigma(:, 1).*sigma(:, 2)];
        M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);
    end

    function M = Jsigmainv(sigma, teps, eps0)
        Nsigma = size(sigma, 1);
        sigmanorm = sqrt(sum(sigma.^2, 2));
        if coff_p > 2 
            zerosnode = (sigmanorm < eps0);
            index = (1:Nsigma);
            index = index(~zerosnode);
            s = sigmanorm(index).^((coff_p-2)/(coff_p-1));
            t = (coff_p - 2)*sigmanorm(index).^(coff_p/(1- coff_p));
            s1 = zeros(Nsigma,1);
            s1(index) = s;           
            s1(zerosnode) = teps^((coff_p-2)/(coff_p-1));
            ii = [1:2*Nsigma, index, index+Nsigma, index, index+Nsigma];
            jj = [1:2*Nsigma, index, index+Nsigma, index + Nsigma, index];
            ss = [s1; s1; t.*sigma(index, 1).^2; t.*sigma(index, 2).^2; t.*sigma(index, 1).*sigma(index, 2); ...
                 t.*sigma(index, 1).*sigma(index, 2)];
            M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);
        elseif coff_p == 2
            M = sparse(1:2*Nsigma, 1:2*Nsigma, ones(2*Nsigma,1), 2*Nsigma, 2*Nsigma);
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
            ii = [1:2*Nsigma, index, index+Nsigma, index, index+Nsigma];
            jj = [1:2*Nsigma, index, index+Nsigma, index + Nsigma, index];
            ss = [s1; s1; t.*sigma(index, 1).^2; t.*sigma(index, 2).^2; t.*sigma(index, 1).*sigma(index, 2); ...
                 t.*sigma(index, 1).*sigma(index, 2)];
            M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);

         
        end
       
    end

    % exact solution
     function s = exactu(p)
        x = p(:,1); y = p(:,2);
        s = coff_a*x.*(x-1).*y.*(y-1);
    end


    function u = g_D(p)
        u = exactu(p);
    end

    function s =  exactsigma(pt)
        gradu = exactgradu(pt);
        gradunorm = sqrt(sum(gradu.^2, 2));
        s = (gradunorm.^(coff_p-2)).*gradu;
    end

    % the derivative of the exact solution
    function s = exactgradu(p)
        x = p(:,1); y = p(:,2);
        s(:, 1) = coff_a*y.*(y-1).*(2*x-1);
        s(:, 2)  = coff_a*x.*(x-1).*(2*y-1);
    end

   
end