% Methods before were not great for iterating over time
% first tried to use solution for 60 mto propagate forward in time but peak
% stayed at 60. Then implemented trust region method. thena adaptable trust
% region method. then found going backwards in time is much better.

%close all;
mu_s = 4.28284e4; % [km^3/s^2]
Rmars = 3390; % [km]
l_s = 18000; % [km]
t_s = sqrt(l_s^3/mu_s); % [s]
v_s = l_s / t_s; % [km/s]
a_s = l_s / t_s^2; % [km/s^2]

N = 3; % number of deployed satallites
tol = 0.001;
Delta = 1000;
umin1 = 0.0001 / a_s;
umin2 = 0.00002 / a_s;


x0 = [20000 / l_s; 0.6; deg2rad(0); deg2rad(0); deg2rad(0); deg2rad(90)]; %keplerian
xf(:,1) = [10000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian
xf(:,2) = [11000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian
xf(:,3) = [12000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian

x0 = keptoModEQ(x0);
xf = keptoModEQ(xf);

t = linspace(0,1.5 * 2 * pi, 120);
dt = t(2) - t(1);
%ntc = 60;
n = length(t);

if 0
    %uguess = zeros(3,n-1,N);
    uguess = u;
    %uguess = U;
    %uguess = zeros(3,n-1,N);
    %uguess(:,1:n-2,:) = U(:,2:n-1,:);

    Usum = zeros(100);
    %U = zeros(3,n-1,N);
    %X = zeros(60, n, N);
    clear U;
    clear X;

    for ntc = 60:60
        nc = ntc;
        w = ones(1,n-1);
        w(1:nc) = 1/3;
        w(nc+1:end) = 1;

        [Usum(nc), U{nc}, X{nc}] = Optimize2(t,x0,xf,nc,uguess,umin1,umin2,w,tol,Delta);
        uguess = U{nc};
    end
end

ntc = 60;

Xc = zeros(6,n,N);
for j = 1:N
    Xc(:,:,j) = modEQtoCart(X{ntc}(1:6,:,j));
end



figure;
hold on;
axis equal;
for j = 1:N
    plot(Xc(1,:,j) * l_s,Xc(2,:,j) * l_s,'LineWidth',2);
end
plot(Xc(1,1:ntc,1) * l_s,Xc(2,1:ntc,1) * l_s,'k-','LineWidth',2);
legend('Agent 1', 'Agent 2', 'Agent 3', 'Carrier');
xlabel('r_1 [km]');
ylabel('r_2 [km]');

figure;
hold on;
for j = 1:N
    stairs(t(1:end-1) * t_s / 3600,vecnorm(U{ntc}(:,:,j) * a_s * 1000),'LineWidth',2)
end
stairs(t(1:ntc) * t_s / 3600, vecnorm(U{ntc}(:,1:ntc,1) * a_s * 1000),'k-','LineWidth',2)
legend('Agent 1', 'Agent 2', 'Agent 3', 'Carrier');
ylabel('||u||_2 [m/s^2]');
xlabel('t [hour]');




