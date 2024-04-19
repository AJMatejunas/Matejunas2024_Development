function d2V_dx2 = func_backDiff(V,h)
% Author: Lloyd Fletcher
% PhotoDyn Group, University of Southampton
% Date: 30/6/2017

for n = 1:length(V)
    if n == 1
        % Use forward difference at the start:
        d2V_dx2(n) = (V(n+1) - V(n)) / h;
    else
        % Use the backward difference everywhere else:
        d2V_dx2(n) = (V(n) - V(n-1)) / h;  
    end
end

