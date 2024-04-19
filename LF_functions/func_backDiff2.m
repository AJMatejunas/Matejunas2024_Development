function d2V_dx2 = func_backDiff2(V,h)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/6/2017

for n = 1:length(V)
    if n == 1 || n == 2
        % Use forward difference at the start:
        % df_dx = (f(n+2) - 2f(n+1) + f(n)) / h^2
        d2V_dx2(n) = (V(n+2) - 2*V(n+1) + V(n)) / h^2;
    else
        % Use the backward difference everywhere else:
        % df_dx = (f(n+1) - 2f(n) + f(n-1)) / h^2 
        d2V_dx2(n) = (V(n) - 2*V(n-1) + V(n-2)) / h^2;  
    end
end

