function d2V_dx2 = func_centralDiff2(V,h)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/6/2017

for n = 1:length(V)
    if n == 1
        % Use forward difference at the start:
        % df_dx = (f(n+2) - 2f(n+1) + f(n)) / h^2
        d2V_dx2(n) = (V(n+2) - 2*V(n+1) + V(n)) / h^2;
    elseif n == length(V)
        % Use backward difference at the end:
        % df_dx = (f(n) - 2f(n-1) + f(n-2)) / h^2
        d2V_dx2(n) = (V(n) - 2*V(n-1) + V(n-2)) / h^2;
    else
        % Use the central difference everywhere else:
        % df_dx = (f(n+1) - 2f(n) + f(n-1)) / h^2 
        d2V_dx2(n) = (V(n+1) - 2*V(n) + V(n-1)) / h^2;  
    end
end

