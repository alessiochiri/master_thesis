function Y_all = ifiss_model_parametric(X)
    persistent A B Bx By G f g xy mv bound
    if isempty(A)
        load('C:\Users\aless\Documents\MATLAB\IFISS\ifiss3.7\ifiss3.7bis\datafiles\square_stokes_nobc.mat', ...
            'A','B','Bx','By','G','f','g','xy','mv','bound');
    end
    N = size(X, 1); % Two inputs: first = viscosity, second = velocity;
    Y_all = zeros(N, 262);
    
    warning('off', 'MATLAB:eigs:SingularA');
    
    for i = 1:N
        viscosity = X(i, 1)/X(i, 2); % viscosity_IFISS = viscosity/V_lid
        vel = X(i, 2); 
        try
            [flowsol, fst, gst, Jnst, Bst] = solve_navier_ifiss_parametric(viscosity, 3, 20, 20, 1e-8, 2, A, B, f, g, xy, mv, bound, vel);
            
            try
                [real_hopf, imag_hopf] = autovalori(xy, bound, Jnst, G, Bst, gst);
                
                if isempty(real_hopf) || isempty(imag_hopf)
                    real_hopf = NaN;
                    imag_hopf = NaN;
                end
            catch ME
                fprintf('  Warning: eigenvalues failure (%s)\n', ME.message);
                real_hopf = NaN;
                imag_hopf = NaN;
            end
            
            profiles = centerline_profiles(xy, flowsol, fst, gst);
            nv = length(fst)/2;
            u_sol = flowsol(1:2*nv);
            K = 0.5 * u_sol' * G * u_sol;
            EN = enstrophy(flowsol, By, Bx, G, xy);
            
        catch ME
            fprintf('  Error: %s\n', ME.message);
            profiles.u_centerline = NaN(129, 1);
            profiles.v_centerline = NaN(129, 1);
            real_hopf = NaN;
            imag_hopf = NaN;
            K = NaN;
            EN = NaN;
        end
        
        % Save
        Y_all(i, 1:129)   = profiles.u_centerline(:)';
        Y_all(i, 130:258) = profiles.v_centerline(:)';
        Y_all(i, 259) = real_hopf;
        Y_all(i, 260) = imag_hopf;
        Y_all(i, 261) = K;
        Y_all(i, 262) = EN;
    end
    
    warning('on', 'MATLAB:eigs:SingularA');
end