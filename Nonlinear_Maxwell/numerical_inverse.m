function [fv,z] = numerical_inverse(f,fd,v,k_bisc,k_newt)
% solve the eqn f(x) = v where v is a vector

h = 5;
ss = 0:h:80;
fss = f(ss)-v;
[f_min,f_min_id] = min(abs(fss),[],2);
ss_est = ss(f_min_id)';

% bisection 
fv = f(ss_est)-v;
z = ss_est;

a = ss_est-h;
b = ss_est+h;
fa = f(a)-v;
fb = f(b)-v;

% find zero points
ida_eq = (abs(fa)<eps);
b(ida_eq) = a(ida_eq);

idb_eq = (abs(fb)<eps);
a(idb_eq) = b(idb_eq);

for k = 1:k_bisc

   c = (a+b)/2;
   fc = f(c)-v;
   idc_eq = (abs(fc)<eps);
   a(idc_eq) = c(idc_eq);
   b(idc_eq) = c(idc_eq);

   fac = fa.*fc;
   fbc = fb.*fc;
   ida_pos = (fac>0);
   a(ida_pos) = c(ida_pos);
   idb_pos= (fbc>0);
   b(idb_pos) = c(idb_pos);
   fa = f(a) - v;
   fb = f(b) - v;
end

sk = a;
fv = f(a)-v;
Newt_eps = 10^(-8);
while k <= k_newt && max(abs(fv)) > Newt_eps
    skp = sk - (f(sk)-v)./(fd(sk));
    sk = skp;
    fv = f(skp)-v;
    k = k+1;
end

z = sk;

end





