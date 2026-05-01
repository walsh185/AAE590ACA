function y = modEQtoCart(x)
%   x is input as [p; f; g; h; k; L]
mu = 1;
p = x(1,:);
f = x(2,:);
g = x(3,:);
h = x(4,:);
k = x(5,:);
L = x(6,:);

e = sqrt(f.^2 + g.^2);
a = p ./ (1 - e.^2);
i = 2 * atan( sqrt(h.^2 + k.^2) );
Omega = mod(atan2(k, h), 2*pi);
omega = mod(atan2(g, f) - Omega, 2*pi);
nu = mod(L - Omega - omega, 2*pi);

p = a .* (1 - e.^2);
r = p ./ (1 + e .* cos(nu));

y(1,:) = r .* (cos(Omega) .* cos(omega + nu) - sin(Omega) .* cos(i) .* sin(omega + nu));
y(2,:) = r .* (sin(Omega) .* cos(omega + nu) + cos(Omega) .* cos(i) .* sin(omega + nu));
y(3,:) = r .* sin(i) .* sin(omega + nu);

vr = sqrt(mu ./ p) .* e .* sin(nu);
vth = sqrt(mu ./ p) .* (1 + e .* cos(nu));

y(4,:) = vr .* (cos(Omega) .* cos(omega + nu) - sin(Omega) .* cos(i) .* sin(omega + nu)) - vth .* (cos(Omega) .* sin(omega + nu) + sin(Omega) .* cos(i) .* cos(omega + nu));
y(5,:) = vr .* (sin(Omega) .* cos(omega + nu) + cos(Omega) .* cos(i) .* sin(omega + nu)) - vth .* (sin(Omega) .* sin(omega + nu) - cos(Omega) .* cos(i) .* cos(omega + nu));
y(6,:) = vr .* sin(i) .* sin(omega + nu) + vth .* sin(i) .* cos(omega + nu);
end

