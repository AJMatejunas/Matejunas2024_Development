function dV_dx = func_centralDiff(V,h)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/6/2017

for n = 1:length(V)
    if n == 1
        % Use forward difference at the start:
        d2V_dx2(n) = (V(n+1) + V(n)) / h;
    elseif n == length(V)
        % Use backward difference at the end:
        d2V_dx2(n) = (V(n) - V(n-1)) / h;
    else
        % Use the central difference everywhere else:
        d2V_dx2(n) = (V(n+1) - V(n-1)) / h;  
    end
end

