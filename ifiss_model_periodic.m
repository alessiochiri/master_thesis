function [Y_all, results] = ifiss_model_periodic(X)
    persistent A B Bx By G f g xy mv bound
    if isempty(A)
        load('C:\Users\aless\Documents\MATLAB\IFISS\ifiss3.7\ifiss3.7bis\datafiles\square_stokes_nobc.mat', ...
        'A','B','Bx','By','G','f','g','xy','mv','bound');
    end

    % Target points for A_max
    points = [
       -0.8, -0.8;   % corner bottom-left
        0.8, -0.8;   % corner bottom-right
       -0.8,  0.8;   % corner top-left
        0.8,  0.8;   % corner top-right
        0.0,  0.0;   % center cavity
        0.0,  0.5;   
        0.0, -0.5;   
       -0.5,  0.0;   
        0.5,  0.0;   
    ];
    n_points = size(points, 1);
    results = struct();
    N = size(X, 1);
    Y_all = zeros(N, 3 + n_points);
    warning('off', 'MATLAB:eigs:SingularA');
    
    for i = 1:N
        Re = X(i, 1);
        viscosity = 2 / Re;
        fprintf('Processing Re = %g\n', Re);
        try
            [flowsol, fst, gst, Jnst, Bst] = solve_navier_ifiss(viscosity, 3, 20, 20, 1e-8, 2, A, B, f, g, xy, mv, bound);
            nv = length(fst)/2;
            maxit_hb = 15;
            tol_hb = 1e-6;
            qmethod = 3;
            [w0, w1c, w1s, w2c, w2s, omega, res_hist, res_blocks, exit_flag, phi_c, phi_s] = harmonicbalanceHB2(viscosity, qmethod, maxit_hb, tol_hb, A, B, G, f, g, xy, mv, bound, flowsol);

            % Extract velocity components
            u0  = w0(1:2*nv);  u1c = w1c(1:2*nv);  u1s = w1s(1:2*nv);
            u2c = w2c(1:2*nv); u2s = w2s(1:2*nv);

            % Mean kinetic Energy(Parseval)
            E = 0.5*(u0'*G*u0) + 0.25*(u1c'*G*u1c + u1s'*G*u1s + u2c'*G*u2c + u2s'*G*u2s);

            % Mean enstrophy (Parseval)
            EN = enstrophy(w0,By,Bx,G,xy) + 0.5*(enstrophy(w1c,By,Bx,G,xy) + enstrophy(w1s,By,Bx,G,xy) + enstrophy(w2c,By,Bx,G,xy) + enstrophy(w2s,By,Bx,G,xy));
   
            idx_nodes = zeros(n_points, 1);
            for k = 1:n_points
                dist = (xy(:,1) - points(k,1)).^2 + (xy(:,2) - points(k,2)).^2;
                [~, idx_nodes(k)] = min(dist);
            end

            % Max amplitude
            n_t = 200;
            t = linspace(0, 2*pi/abs(omega), n_t);
            A_max_nodes = zeros(n_points, 1);
            for k = 1:n_points
                idx = idx_nodes(k);
                u_t = u0(idx)    + u1c(idx)   *cos(abs(omega)*t) - u1s(idx)   *sin(abs(omega)*t) + u2c(idx)   *cos(2*abs(omega)*t) - u2s(idx) *sin(2*abs(omega)*t);
                v_t = u0(idx+nv) + u1c(idx+nv)*cos(abs(omega)*t) - u1s(idx+nv)*sin(abs(omega)*t) + u2c(idx+nv)*cos(2*abs(omega)*t) - u2s(idx+nv)*sin(2*abs(omega)*t);
                A_max_nodes(k) = max(sqrt(u_t.^2 + v_t.^2));
            end

        
        catch ME
            fprintf('  Error: %s\n', ME.message);
            E = NaN; EN = NaN; omega = NaN;
            A_max_nodes = NaN(n_points, 1);
        end
        

        Y_all(i, 1) = E;
        Y_all(i, 2) = EN;
        Y_all(i, 3) = abs(omega);
        Y_all(i, 4:3+n_points) = A_max_nodes';
        results(i).Re = X(i,1);
        results(i).omega      = omega;
        results(i).res_hist   = res_hist;
        results(i).res_blocks = res_blocks;
        results(i).n_iter     = length(res_blocks);  
        results(i).exit_flag  = exit_flag;
    end
    warning('on', 'MATLAB:eigs:SingularA');
end