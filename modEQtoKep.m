function [y] = modEqToKep(x)
%   x is input as [p; f; g; h; k; L]
%   y is output as [a; e; i; Omega; omega; nu]
p = x(1,:);
f = x(2,:);
g = x(3,:);
h = x(4,:);
k = x(5,:);
L = x(6,:);

e = sqrt(f.^2 + g.^2);
a = p ./ (1 - e.^2);
i = atan2(2 * sqrt(h.^2 + k.^2), 1 - h.^2 - k.^2);
Omega = atan2(k, h);
omega = atan2(g, f) - Omega;
nu = L - Omega - omega;

y = [a; e; i; Omega; omega; nu];
end

