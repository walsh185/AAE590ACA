close all;
mu_s = 4.28284e4; % [km^3/s^2]
Rmars = 3390; % [km]
l_s = 18000; % [km]
t_s = sqrt(l_s^3/mu_s); % [s]
v_s = l_s / t_s; % [km/s]
a_s = l_s / t_s^2; % [km/s^2]

N = 3; % number of deployed satallites
K = 5; % number of convex iterations
umin1 = 0.0001;
umin2 = 0.00002;

x0 = [20000 / l_s; 0.6; deg2rad(0); deg2rad(0); deg2rad(0); deg2rad(90)]; %keplerian
xf(:,1) = [10000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian
xf(:,2) = [11000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian
xf(:,3) = [12000 / l_s; 0.05; deg2rad(0); deg2rad(0); deg2rad(0); 0]; %keplerian

x0 = keptoModEQ(x0);
xf = keptoModEQ(xf);

t = linspace(0,1.5 * 2 * pi, 120);
n = length(t);


for nc = 60:60
    w = ones(1,n-1);
    w(1:nc-1) = 1/3;   % cheaper early control
    w(nc:end) = 1;     % normal cost later    

    X = zeros(60, n, N);
    Xc = zeros(6, n, N);
    A_k = zeros(6,6,n,N);
    B_k = zeros(6,3,n,N);
    c_k = zeros(6,n,N);
    u = zeros(3, N);
    for j = 1:N
        X(1:6, 1, j) = x0;
    end
    
    for j = 1:N
        for i = 2:n
            ic = [X(1:6,i-1,j); reshape(eye(6),36,1); zeros(18,1)];
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
           x_ode = ode45(@(t,x) odeSingle(t,x,u(:,j)), [t(i-1), t(i)], ic, opts);
           X(:,i,j) = deval(x_ode, t(i));
        end
    
       Xc(:,:,j) = modEQtoCart(X(1:6,:,j));
      A_k(:,:,:,j) = reshape(X(7:42,:,j),6,6,[]);
      B_k(:,:,:,j) = reshape(X(43:60,:,j),6,3,[]);
      for i = 2:n
          c_k(:,i,j) = X(1:6,i,j) - A_k(:,:,i,j) * X(1:6,i-1,j) - B_k(:,:,i,j) * u(:,j);
      end
   end


    cvx_begin
       variable x(6,n,N)
        variable u(3,n-1,N)
       cost = 0;
       for j = 1:N
           cost = cost + sum(norms(u(:,:,j),2,1));
       end
       minimize cost
       subject to 
            for j = 1:N
                x(:,1,j) == x0;
                x(1:5,n,j) == xf(1:5,j);
                    for i = 2:n
                    x(:,i,j) == A_k(:,:,i,j) * x(:,i-1,j) + B_k(:,:,i,j) * u(:,i-1,j) + c_k(:,i,j);
                    norm(u(:,i-1,j),2) * a_s <= 0.0001;
                end
            end
        

    cvx_end


    for k = 1:K
      X = zeros(60, n, N);
      Xc = zeros(6, n, N);
      A_k = zeros(6,6,n,N);
      B_k = zeros(6,3,n,N);
      c_k = zeros(6,n,N);

        for j = 1:N

            X(1:6, 1, j) = x0;

          for i = 2:n
                ic = [X(1:6,i-1,j); reshape(eye(6),36,1); zeros(18,1)];
                opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
                x_ode = ode45(@(t,x) odeSingle(t,x,u(:,i-1,j)), [t(i-1), t(i)], ic, opts);
                X(:,i,j) = deval(x_ode, t(i));
          end
    
            Xc(:,:,j) = modEQtoCart(X(1:6,:,j));
            A_k(:,:,:,j) = reshape(X(7:42,:,j),6,6,[]);
            B_k(:,:,:,j) = reshape(X(43:60,:,j),6,3,[]);
            for i = 2:n
                c_k(:,i,j) = X(1:6,i,j) - A_k(:,:,i,j) * X(1:6,i-1,j) - B_k(:,:,i,j) * u(:,i-1,j);
            end
        end


        cvx_begin
             variable x(6,n,N)
             variable u(3,n-1,N)
             cost = 0;
             for j = 1:N
                  cost = cost + sum(w .* norms(u(:,:,j),2,1));
             end
            minimize cost
            subject to 
                for j = 1:N
                    x(:,1,j) == x0;
                    x(1:5,n,j) == xf(1:5,j);
                    for i = 2:n
                        x(:,i,j) == A_k(:,:,i,j) * x(:,i-1,j) + B_k(:,:,i,j) * u(:,i-1,j) + c_k(:,i,j);
                    end
                    for i = 1:nc-1
                        u(:,i,j) == u(:,i,1);
                        norm(u(:,i,j),2) * a_s <= umin1;
                    end
                    for i = nc:n-1
                        norm(u(:,i,j),2) * a_s <= umin2;
                    end
                end
        
        cvx_end
    end
    
    costs(nc) = cost
end
