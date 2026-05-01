function y = odeSingle(t,x,u)
    % Assume mu = 1
    p = x(1,:);
    f = x(2,:);
    g = x(3,:);
    h = x(4,:);
    k = x(5,:);
    L = x(6,:);
    q = 1 + f * cos(L) + g * sin(L);
    s2 = 1 + h^2 + k^2;
    y = zeros(60,1);

    STM = reshape(x(7:42),6,6);
    STMb = reshape(x(43:60),6,3);

    y(1:6) = p^(1/2) * [0, 2*p/q, 0;...
        sin(L), ((q+1)*cos(L)+f)/q, -g*(h*sin(L)-k*cos(L))/q;...
        -cos(L), ((q+1)*sin(L)+g)/q, f*(h*sin(L)-k*cos(L))/q;...
        0, 0, s2*cos(L)/2/q;...
        0, 0, s2*sin(L)/2/q;...
        0, 0, (h*sin(L)-k*cos(L))/q] * u;
    y(6) = y(6) + q^2 / p^(3/2);
    y(7:42) = reshape(A_func(x(1:6), u) * STM, 36, 1);
    y(43:60) = reshape(A_func(x(1:6), u) * STMb + B_func(x(1:6), u), 18, 1);
end

