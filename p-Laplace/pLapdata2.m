function pde = pLapdata2(coff_p)
%% p-Laplacian equation example 2

nueps = 0.0001;
if coff_p < 2
    Jeps0 = 0.0001;
else
    Jeps0 = eps;
end

pde = struct('f',@f, 'exactsigma',@exactsigma,'exactu',@exactu,'g_D',@g_D, 'K', @K,...
          'Ku', @Ku, 'Ju', @Ju, 'Jdu', @Jdu, 'Jsigma', @Jsigma, 'Jsigmainv', @Jsigmainv,...
          'coff_p', coff_p, 'exactgradu', @exactgradu);

    % load data (right hand side function)
    function rhs =  f(pt)
        x = pt(:,1);
        rhs = ones(length(x),1);
    end

    function s = Ku(gradunorm) % nonlinear coefficient
        if coff_p  < 2
           gradunorm = gradunorm + (gradunorm < eps)*nueps;
        end
        s = gradunorm.^(coff_p - 2);
    end

   
    function s = Jdu(gradunorm) % nonlinear coefficient
        s = gradunorm.^(coff_p - 4);
    end

    function M = Ju(gradu)
        % Regularized Jacobian of sigma = |grad u|^(p-2) grad u.
        NT = size(gradu,1);
        gradunorm = sqrt(sum(gradu.^2,2));
        
        if coff_p  > 2
            s = gradunorm.^(coff_p - 2) + nueps;
            zerosnode = gradunorm.^(coff_p - 2) < Jeps0;
        elseif coff_p < 2
            s = (gradunorm + (gradunorm < eps).*nueps).^(coff_p - 2);
            zerosnode = (gradunorm < Jeps0);
        end
        
        t = zeros(NT, 1);
        t(~zerosnode) = (coff_p - 2)*gradunorm(~zerosnode).^(coff_p - 4);

        j11 = s + t.*gradu(:,1).^2;
        j22 = s + t.*gradu(:,2).^2;
        j12 = t.*gradu(:,1).*gradu(:,2);
        ii = [1:NT, NT+1:2*NT, 1:NT, NT+1:2*NT];
        jj = [1:NT, NT+1:2*NT, NT+1:2*NT, 1:NT];
        M = sparse(ii,jj,[j11; j22; j12; j12],2*NT,2*NT);
    end



    function s = K(sigmanorm) % nonlinear coefficient
        if coff_p  > 2
            sigmanorm = sigmanorm + (sigmanorm < eps)*nueps;
        end
        s = sigmanorm.^(-(coff_p - 2)/(coff_p-1));
    end

    function M = Jsigma(sigma)
        % sigma is dim*2 matrix
        Nsigma = size(sigma, 1);
        sigmanorm = sqrt(sum(sigma.^2, 2));
        if coff_p > 2
            sigmanorm = sigmanorm + (sigmanorm < eps)*0.0001;
        end
        s = sigmanorm.^(-(coff_p-2)/(coff_p-1));
        t = -(coff_p -2)/(coff_p -1)*(sigmanorm.^(-(3*coff_p-4)/(coff_p-1)));
        ii = [1:2*Nsigma, 1:2*Nsigma, 1:Nsigma, Nsigma+1:2*Nsigma];
        jj = [1:2*Nsigma, 1:2*Nsigma, Nsigma+1:2*Nsigma, 1:Nsigma];
        ss = [s; s; t.*sigma(:, 1).^2; t.*sigma(:, 2).^2; t.*sigma(:, 1).*sigma(:, 2); ...
             t.*sigma(:, 1).*sigma(:, 2)];
        M = sparse(ii, jj, ss, 2*Nsigma, 2*Nsigma);
    end

    function M = Jsigmainv(sigma, teps, eps0)
         % sigma is dim*2 matrix
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
        s = (1/2)^(1/(coff_p-1))*((coff_p-1)/coff_p)*(1-((sqrt(x.^2 + y.^2)).^(coff_p/(coff_p-1))));
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
        beta = -(1/2)^(1/(coff_p-1))*(sqrt(x.^2+y.^2).^((2-coff_p)/(coff_p-1)));
        s(:, 1) = beta.*x;
        s(:, 2)  = beta.*y;
    end

   
end