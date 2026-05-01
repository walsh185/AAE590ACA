function [Usum, U, X] = Optimize(t, x0, Xf, ntc, uguess, umin1, umin2, w, maxTol, Delta)
    % INPUT
    % t = timeseries for convex optimization
    % x0 = starting state
    % Xf = array of final states
    % ntc = index of t for deployment
    % uguess = guess solution for control
    % w = weights for cost function
    % tol = tolerance for solution accuracy

    Nt = length(t);
    dt = zeros(1,Nt-1);
    for i = 2:Nt
        dt(i-1) = t(i) - t(i-1);
    end
    [~,N] = size(Xf);

    tol = inf;
    u = uguess;

    X = zeros(60, Nt, N);
    Xc = zeros(6, Nt, N);
    A_k = zeros(6,6,Nt,N);
    B_k = zeros(6,3,Nt,N);
    c_k = zeros(6,Nt,N);
    for j = 1:N
        X(1:6, 1, j) = x0;
    end

    for j = 1:N
        for i = 2:Nt
            ic = [X(1:6,i-1,j); reshape(eye(6),36,1); zeros(18,1)];
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
            x_ode = ode45(@(t,x) odeSingle(t,x,u(:,i-1,j)), [t(i-1), t(i)], ic, opts);
            X(:,i,j) = deval(x_ode, t(i));
        end
    
        Xc(:,:,j) = modEQtoCart(X(1:6,:,j));
        A_k(:,:,:,j) = reshape(X(7:42,:,j),6,6,[]);
        B_k(:,:,:,j) = reshape(X(43:60,:,j),6,3,[]);
        for i = 2:Nt
            c_k(:,i,j) = X(1:6,i,j) - A_k(:,:,i,j) * X(1:6,i-1,j) - B_k(:,:,i,j) * u(:,i-1,j);
        end
    end


    while tol > maxTol

        u_prev = u;
                
        cvx_begin
            variable x(6,Nt,N)
            variable u(3,Nt-1,N)
            cost = 0;
            for j = 1:N
                cost = cost + sum(dt .* w .* norms(u(:,:,j),2,1));
            end
            minimize cost
            subject to 
                for j = 1:N
                    x(:,1,j) == x0;
                    x(1:5,Nt,j) == Xf(1:5,j);
                    for i = 2:Nt
                        x(:,i,j) == A_k(:,:,i,j) * x(:,i-1,j) + B_k(:,:,i,j) * u(:,i-1,j) + c_k(:,i,j);
                    end
                    for i = 1:ntc-1
                        u(:,i,j) == u(:,i,1);
                        norm(u(:,i,j),2) <= umin1;
                    end
                    for i = ntc:Nt-1
                        norm(u(:,i,j),2) <= umin2;
                        norm(u(:,i,j) - u_prev(:,i,j),2) <= Delta;
                    end
                end

                for i = 1:ntc-1
                    norm(u(:,i,1) - u_prev(:,i,1),2) <= Delta;
                end
        
        cvx_end
        
        Usum = cvx_optval;
        Delta = Delta / 10;
        
        for j = 1:N
            X(1:6, 1, j) = x0;
        end

        for j = 1:N
            for i = 2:Nt
                ic = [X(1:6,i-1,j); reshape(eye(6),36,1); zeros(18,1)];
                opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
                x_ode = ode45(@(t,x) odeSingle(t,x,u(:,i-1,j)), [t(i-1), t(i)], ic, opts);
                X(:,i,j) = deval(x_ode, t(i));
            end
    
            Xc(:,:,j) = modEQtoCart(X(1:6,:,j));
            A_k(:,:,:,j) = reshape(X(7:42,:,j),6,6,[]);
            B_k(:,:,:,j) = reshape(X(43:60,:,j),6,3,[]);
            for i = 2:Nt
                c_k(:,i,j) = X(1:6,i,j) - A_k(:,:,i,j) * X(1:6,i-1,j) - B_k(:,:,i,j) * u(:,i-1,j);
            end
        end

        tol = 0;

        

        for j = 1:N
            tol = tol + norm(X(1:5,end,j) - Xf(1:5,j));
        end
        tol
    end
    
    U = u;
end

