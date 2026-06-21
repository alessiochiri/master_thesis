function [w0, w1c, w1s, w2c, w2s, omega, res_hist, res_blocks, exit_flag, phi_c, phi_s] = harmonicbalanceHB2(viscosity, qmethod, ...
    maxit_hb, tol_hb, A, B, G, f, g, xy, mv, bound, flowsol)
% HB Harmonic Balance N_H=2 for Navier-Stokes (IFISS)
%
% Equation: M*dw/dt + L*w + N(w,w) = 0 
% Ansatz:  w_hat_1 = (1/2)*(w1c + i*w1s); w_hat_2=(1/2)*(w2c + i*w2s)
% N(w,w) = w*grad(w)
%
% Residuals:
%   R0  = [nu*A + N(w0)]*w0 - fext + 0.5*[N(w1c,w1c)+N(w1s,w1s)+N(w2c,w2c)+N(w2s,w2s)]
%   R1c = Lh*w1c + N(w0,w1c) + N(w1c,w0) + 0.5*(N(w1c,w2c) + N(w2c,w1c) + N(w1s,w2s) + N(w2s, w1s)) - om*M*w1s
%   R1s = Lh*w1s + N(w0,w1s) + N(w1s,w0) + 0.5*(N(w1c,w2s) + N(w2s,w1c) - N(w1s,w2c) - N(w2c, w1s)) + om*M*w1c
%   R2c = Lh*w2c + N(w0,w2c) + N(w2c,w0) + 0.5*(N(w1c,w1c) - N(w1s,w1s)) - 2*om*M*w2s
%   R2s = Lh*w2s + N(w0,w2s) + N(w2s,w0) + 0.5*(N(w1c,w1s) + N(w1s,w1c)) + 2*om*M*w2c
%   Rphi = phi_c(1:2*nv)' * G * w1c(1:2*nv) - eps  
% Augmented vector state:
%   X_aug = [w0; lam0; w1c; lam1c; w1s; lam1s; w2c; lam2c; w2s; lam2s; omega]
%  where lam* are the multipliers for the pressure pinning

fprintf('\n=== Harmonic Balance N_H=2 (IFISS) ===\n')

%% --- Size ---
nv  = size(xy,1);
np  = length(g);
nw  = 2*nv + np;
nwa = nw + 1;  % size augmented block

%% --- Matrices ---
Mnst = mass_q2_bc(G, xy, bound);
Mext = [Mnst,            sparse(2*nv, np);
         sparse(np, 2*nv), sparse(np,   np)];

[Ast_h, Bst_h, ~, ~] = flowbc(A, B, f*0, g*0, xy, bound);

[~, Bst, ~, gst] = flowbc(A, B, f, g, xy, bound);
Lmat_h = [viscosity*Ast_h, Bst_h'; Bst_h, sparse(np,np)];

%% --- Phase 1: Hopf eigenvector ---
fprintf('\nFase 1: Computation of Hopf eigenvetor ...\n')
[Nxx, Nxy, Nyx, Nyy] = newton_q2(xy, mv, flowsol);
Nst  = navier_q2(xy, mv, flowsol);
J    = viscosity*A + [Nst+Nxx, Nxy; Nyx, Nst+Nyy];
Jnst = newtonbc(J, xy, bound);
[real_hopf, imag_hopf, phi_c, phi_s] = eigenvalvec(xy, bound, Jnst, G, Bst, gst);
omega0 = abs(imag_hopf);
fprintf('omega - initial guess: %.6f\n', omega0)
fprintf('Re(lambda) : %.6f\n', real_hopf)

% Normalization
phi_vel_c = phi_c(1:2*nv);
phi_vel_s = phi_s(1:2*nv);

E_phi = sqrt(phi_vel_c'*G*phi_vel_c);
phi_c = phi_c / E_phi;
phi_s = phi_s / E_phi;

w0    = flowsol;
eps   = 0.02;
w1c   = eps * phi_c;
w1s   = eps * phi_s;
w2c   = 0.01 * eps * phi_c;
w2s   = 0.01 * eps * phi_s;
omega = omega0;

X_aug = [w0; 0; w1c; 0; w1s; 0; w2c; 0; w2s; 0; omega];

% Extraction Indices
idx_w0  = 1:nw;
idx_l0  = nw+1;
idx_w1c = nw+2 : 2*nw+1;
idx_l1c = 2*nw+2;
idx_w1s = 2*nw+3 : 3*nw+2;
idx_l1s = 3*nw+3;
idx_w2c = 3*nw+4 : 4*nw+3;
idx_l2c = 4*nw+4;
idx_w2s = 4*nw+5 : 5*nw+4;
idx_l2s = 5*nw+5;
idx_om  = 5*nw+6;

%% Newton HB 
fprintf('\nPhase 3: Newton HB...\n\n')
res_hist   = zeros(maxit_hb, 1);
res_blocks = struct('R0',{},'R1c',{},'R1s',{},'R2c',{},'R2s',{},'Rphi',{},'Rtot',{});
exit_flag  = 'maxit';  % default: esaurisce iterazioni senza convergere

for it = 1:maxit_hb

    w0_k  = X_aug(idx_w0);
    w1c_k = X_aug(idx_w1c);
    w1s_k = X_aug(idx_w1s);
    w2c_k = X_aug(idx_w2c);
    w2s_k = X_aug(idx_w2s);
    om_k  = X_aug(idx_om);

    %% Residuals
    [R0, R1c, R1s, R2c, R2s, Rphi] = compute_residual_blocks(w0_k, w1c_k, w1s_k, w2c_k, w2s_k, om_k, ...
        viscosity, A, B, Lmat_h, Mext, f, g, nv, np, xy, mv, bound, Mnst, phi_c, phi_s, eps, G);

    % Augmented residual
    R_aug    = [R0; 0; R1c; 0; R1s; 0; R2c; 0; R2s; 0; Rphi];
    res_norm = norm(R_aug);

    % 
    res_blocks(it).R0   = norm(R0);
    res_blocks(it).R1c  = norm(R1c);
    res_blocks(it).R1s  = norm(R1s);
    res_blocks(it).R2c  = norm(R2c);
    res_blocks(it).R2s  = norm(R2s);
    res_blocks(it).Rphi = abs(Rphi);
    res_blocks(it).Rtot = res_norm;
    res_blocks(it).omega = om_k;        
    res_blocks(it).condJ = NaN;         
    res_blocks(it).dX_norm = NaN;       
    fprintf('Newton HB iter %3d — Rtot=%.3e  R0=%.3e  R1c=%.3e  R1s=%.3e  R2c=%.3e  R2s=%.3e  Rphi=%.3e  omega=%.4f\n', it, res_norm, norm(R0), norm(R1c), norm(R1s), norm(R2c), norm(R2s), abs(Rphi), om_k);

    % Stagnation
    if it > 2 && res_norm > 0.99 * res_hist(it-1)
        exit_flag = 'stagnation';
        fprintf('  [STAGNATION] iter %d — res=%.3e not >1%% — exit\n', it, res_norm);
        break
    end

    res_hist(it) = res_norm;

    % Convergence criterion
    if res_norm < tol_hb
        exit_flag = 'residual_tol';
        fprintf('  [CONVERGED: residual_tol] in %d iterations\n', it)
        break
    end

    %% Jacobiano 
    Nw0  = navier_q2(xy, mv, w0_k);
    Nw1c = navier_q2(xy, mv, w1c_k);
    Nw1s = navier_q2(xy, mv, w1s_k);
    Nw2c = navier_q2(xy, mv, w2c_k);
    Nw2s = navier_q2(xy, mv, w2s_k);

    [J00_uu, J00_uv, J00_vu, J00_vv] = newton_q2(xy, mv, w0_k);
    [J1c_uu, J1c_uv, J1c_vu, J1c_vv] = newton_q2(xy, mv, w1c_k);
    [J1s_uu, J1s_uv, J1s_vu, J1s_vv] = newton_q2(xy, mv, w1s_k);
    [J2c_uu, J2c_uv, J2c_vu, J2c_vv] = newton_q2(xy, mv, w2c_k);
    [J2s_uu, J2s_uv, J2s_vu, J2s_vv] = newton_q2(xy, mv, w2s_k);

    Mext_aug = [Mnst,           sparse(2*nv,np), sparse(2*nv,1);
                 sparse(np,2*nv), sparse(np,np),   sparse(np,1);
                 sparse(1,2*nv),  sparse(1,np),    0];

    % dR0/dw0
    J00_vel  = newtonbc(viscosity*A + [Nw0+J00_uu, J00_uv; J00_vu, Nw0+J00_vv], xy, bound);
    dR0dw0   = [J00_vel,        Bst',           zeros(2*nv,1);
                Bst,            sparse(np,np),   ones(np,1)/np;
                zeros(1,2*nv),  ones(1,np)/np,   0];

    J1c_vel  = offdiagbc([Nw1c+J1c_uu, J1c_uv; J1c_vu, Nw1c+J1c_vv], nv, bound);
    J1s_vel  = offdiagbc([Nw1s+J1s_uu, J1s_uv; J1s_vu, Nw1s+J1s_vv], nv, bound);
    J2c_vel  = offdiagbc([Nw2c+J2c_uu, J2c_uv; J2c_vu, Nw2c+J2c_vv], nv, bound);
    J2s_vel  = offdiagbc([Nw2s+J2s_uu, J2s_uv; J2s_vu, Nw2s+J2s_vv], nv, bound);

    dR0dw1c = 0.5*[J1c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR0dw1s = 0.5*[J1s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR0dw2c = 0.5*[J2c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR0dw2s = 0.5*[J2s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR0dom  = sparse(nwa, 1);

    dR1cdw0 = [offdiagbc([Nw1c+J1c_uu, J1c_uv; J1c_vu, Nw1c+J1c_vv], nv, bound), sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    J1C = viscosity*A + [Nw0+J00_uu, J00_uv; J00_vu, Nw0+J00_vv] + 0.5*[Nw2c+J2c_uu, J2c_uv; J2c_vu, Nw2c+J2c_vv];
    J1C_bc = newtonbc(J1C, xy, bound);
    dR1cdw1c = [J1C_bc, Bst_h', zeros(2*nv,1); Bst_h, sparse(np,np), ones(np,1)/np; zeros(1,2*nv), ones(1,np)/np, 0];
   
    J2s_aug  = [J2s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR1cdw1s = -om_k * Mext_aug + 0.5*J2s_aug;

    dR1cdom  = [-Mnst*w1s_k(1:2*nv); sparse(np,1); 0];

    dR1cdw2c = 0.5*[J1c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR1cdw2s = 0.5*[J1s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    dR1sdw0 = [offdiagbc([Nw1s+J1s_uu, J1s_uv; J1s_vu, Nw1s+J1s_vv], nv, bound), sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    dR1sdw1c = om_k * Mext_aug + 0.5*J2s_aug;
    dR1sdom  = [Mnst*w1c_k(1:2*nv); sparse(np,1); 0];

    J1S = viscosity*A + [Nw0+J00_uu, J00_uv; J00_vu, Nw0+J00_vv] - 0.5*[Nw2c+J2c_uu, J2c_uv; J2c_vu, Nw2c+J2c_vv];
    J1S_bc = newtonbc(J1S, xy, bound);
    dR1sdw1s = [J1S_bc, Bst_h', zeros(2*nv,1); Bst_h, sparse(np,np), ones(np,1)/np; zeros(1,2*nv), ones(1,np)/np, 0];
    dR1sdw2c = -0.5*[J1s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR1sdw2s =  0.5*[J1c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    dR2cdw0 = [offdiagbc([Nw2c+J2c_uu, J2c_uv; J2c_vu, Nw2c+J2c_vv], nv, bound), sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    dR2cdw1c =  0.5*[J1c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR2cdw1s = -0.5*[J1s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];

    J2C = viscosity*A + [Nw0+J00_uu, J00_uv; J00_vu, Nw0+J00_vv];
    J2C_bc = newtonbc(J2C, xy, bound);
    dR2cdw2c = [J2C_bc, Bst_h', zeros(2*nv,1); Bst_h, sparse(np,np), ones(np,1)/np; zeros(1,2*nv), ones(1,np)/np, 0];
    dR2cdw2s = -2*om_k * Mext_aug;
    dR2cdom  = [-2*Mnst*w2s_k(1:2*nv); sparse(np,1); 0];

    dR2sdw0 = [offdiagbc([Nw2s+J2s_uu, J2s_uv; J2s_vu, Nw2s+J2s_vv], nv, bound), sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR2sdw1c =  0.5*[J1s_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR2sdw1s =  0.5*[J1c_vel, sparse(2*nv,np), sparse(2*nv,1); sparse(np,2*nv), sparse(np,np), sparse(np,1); sparse(1,2*nv), sparse(1,np), 0];
    dR2sdw2c = 2*om_k * Mext_aug;

    J2S = viscosity*A + [Nw0+J00_uu, J00_uv; J00_vu, Nw0+J00_vv];
    J2S_bc = newtonbc(J2S, xy, bound);
    dR2sdw2s = [J2S_bc, Bst_h', zeros(2*nv,1); Bst_h, sparse(np,np), ones(np,1)/np; zeros(1,2*nv), ones(1,np)/np, 0];
    dR2sdom  = [2*Mnst*w2c_k(1:2*nv); sparse(np,1); 0];

    dRphi_dw0  = sparse(1, nwa);
    dRphi_dw1c = [phi_s(1:2*nv)'*G, sparse(1,np), 0];
    dRphi_dw1s = [phi_c(1:2*nv)'*G, sparse(1,np), 0];
    dRphi_dw2c = sparse(1, nwa);
    dRphi_dw2s = sparse(1, nwa);
    dRphi_dom  = 0;
    dRphi = [dRphi_dw0, dRphi_dw1c, dRphi_dw1s, dRphi_dw2c, dRphi_dw2s, dRphi_dom];

    J_HB = [dR0dw0,  dR0dw1c,  dR0dw1s,  dR0dw2c,  dR0dw2s,  dR0dom;
            dR1cdw0, dR1cdw1c, dR1cdw1s, dR1cdw2c, dR1cdw2s, dR1cdom;
            dR1sdw0, dR1sdw1c, dR1sdw1s, dR1sdw2c, dR1sdw2s, dR1sdom;
            dR2cdw0, dR2cdw1c, dR2cdw1s, dR2cdw2c, dR2cdw2s, dR2cdom;
            dR2sdw0, dR2sdw1c, dR2sdw1s, dR2sdw2c, dR2sdw2s, dR2sdom;
            dRphi];

    %% Newton step with backtracking
    t_pre  = tic;
    dX_aug = -J_HB \ R_aug;
    fprintf('  [solve] %.2f s\n', toc(t_pre))
    
    dX_n = norm(dX_aug);
    cJ   = dX_n / (norm(R_aug) + 1e-14);
    res_blocks(it).condJ   = cJ;
    res_blocks(it).dX_norm = dX_n;
    fprintf('  proxy_cond    = %.3e\n', cJ);
    fprintf('  norm(dX)      = %.3e\n', dX_n);
    alpha = 1.0;
    X_new = X_aug + dX_aug;
    for ls = 1:20
        w0_n  = X_new(idx_w0);
        w1c_n = X_new(idx_w1c); w1c_n([bound; bound+nv]) = 0;
        w1s_n = X_new(idx_w1s); w1s_n([bound; bound+nv]) = 0;
        w2c_n = X_new(idx_w2c); w2c_n([bound; bound+nv]) = 0;
        w2s_n = X_new(idx_w2s); w2s_n([bound; bound+nv]) = 0;
        om_n  = X_new(idx_om);
        [R0n,R1cn,R1sn,R2cn,R2sn,Rphin] = compute_residual_blocks(w0_n,w1c_n,w1s_n,w2c_n,w2s_n,om_n, viscosity,A,B,Lmat_h,Mext,f,g,nv,np,xy,mv,bound,Mnst,phi_c,phi_s,eps,G);
        R_new = [R0n;0;R1cn;0;R1sn;0;R2cn;0;R2sn;0;Rphin];
        if norm(R_new) < res_norm && om_n > 0.5 * omega0; break; end
        alpha = alpha / 2;
        X_new = X_aug + alpha * dX_aug;
        fprintf('  [backtrack] ls=%d  alpha=%.5f  res=%.3e  omega=%.4f\n', ls, alpha, norm(R_new), om_n)
    end
    X_aug = X_new;
    dX_norm = norm(dX_aug) / (norm(X_aug) + 1e-14);
    if dX_norm < tol_hb
        exit_flag = 'increment_tol';
        fprintf('  [CONVERGED: increment_tol] in %d iter\n', it)
        break
    end

    fprintf('  norm(w1c)=%.3e  norm(w1s)=%.3e  norm(w2c)=%.3e  norm(w2s)=%.3e  omega=%.6f\n', ...
        norm(X_aug(idx_w1c)), norm(X_aug(idx_w1s)), norm(X_aug(idx_w2c)), norm(X_aug(idx_w2s)), X_aug(idx_om))
end

%% Output 
res_hist   = res_hist(1:it);
res_blocks = res_blocks(1:it);

w0    = X_aug(idx_w0);
w1c   = X_aug(idx_w1c);
w1s   = X_aug(idx_w1s);
w2c   = X_aug(idx_w2c);
w2s   = X_aug(idx_w2s);
omega = X_aug(idx_om);

fprintf('\n=== HB completed [%s] ===\n', exit_flag)
fprintf('N. iter : %d\n',   it)
fprintf('omega   : %.6f\n', omega)
fprintf('Amplitude w1c  : %.6e\n', norm(w1c))
fprintf('Amplitude w1s  : %.6e\n', norm(w1s))
fprintf('Amplitude w2c  : %.6e\n', norm(w2c))
fprintf('Amplitude w2s  : %.6e\n', norm(w2s))
end

%% RESIDUALS

function [R0, R1c, R1s, R2c, R2s, Rphi] = compute_residual_blocks(w0_k, w1c_k, w1s_k, w2c_k, w2s_k, om_k, viscosity, A, B, Lmat_h, Mext, f, g, nv, np, xy, mv, bound, Mnst, phi_c, phi_s, eps, G)

    Nw0   = navier_q2(xy, mv, w0_k);
    Anst  = viscosity*A + [Nw0, sparse(nv,nv); sparse(nv,nv), Nw0];
    [Anst_bc, Bst_k, fst_k, gst_k] = flowbc(Anst, B, f, g, xy, bound);
    Lmat_0 = [Anst_bc, Bst_k'; Bst_k, sparse(np,np)];

    N11cc = assemble_N(w1c_k, w1c_k, nv, np, xy, mv);
    N11ss = assemble_N(w1s_k, w1s_k, nv, np, xy, mv);
    N22cc = assemble_N(w2c_k, w2c_k, nv, np, xy, mv);
    N22ss = assemble_N(w2s_k, w2s_k, nv, np, xy, mv);
    N_harmonics = 0.5*(N11cc + N11ss + N22cc + N22ss);
    N_harmonics([bound; bound+nv]) = 0;
    R0 = Lmat_0*w0_k + N_harmonics - [fst_k; gst_k];

    N01c  = assemble_N(w0_k,  w1c_k, nv, np, xy, mv) + assemble_N(w1c_k, w0_k,  nv, np, xy, mv);
    N01s  = assemble_N(w0_k,  w1s_k, nv, np, xy, mv) + assemble_N(w1s_k, w0_k,  nv, np, xy, mv);
    N1c2c = assemble_N(w1c_k, w2c_k, nv, np, xy, mv) + assemble_N(w2c_k, w1c_k, nv, np, xy, mv);
    N1s2s = assemble_N(w1s_k, w2s_k, nv, np, xy, mv) + assemble_N(w2s_k, w1s_k, nv, np, xy, mv);
    N1c2s = assemble_N(w1c_k, w2s_k, nv, np, xy, mv) + assemble_N(w2s_k, w1c_k, nv, np, xy, mv);
    N1s2c = assemble_N(w1s_k, w2c_k, nv, np, xy, mv) + assemble_N(w2c_k, w1s_k, nv, np, xy, mv);

    R1c = Lmat_h*w1c_k + N01c + 0.5*N1c2c + 0.5*N1s2s - om_k*(Mext*w1s_k);
    R1s = Lmat_h*w1s_k + N01s + 0.5*N1c2s - 0.5*N1s2c + om_k*(Mext*w1c_k);

    N02c  = assemble_N(w0_k,  w2c_k, nv, np, xy, mv) + assemble_N(w2c_k, w0_k,  nv, np, xy, mv);
    N02s  = assemble_N(w0_k,  w2s_k, nv, np, xy, mv) + assemble_N(w2s_k, w0_k,  nv, np, xy, mv);
    N1c1s = assemble_N(w1c_k, w1s_k, nv, np, xy, mv) + assemble_N(w1s_k, w1c_k, nv, np, xy, mv);

    R2c = Lmat_h*w2c_k + N02c + 0.5*(N11cc - N11ss) - 2*om_k*(Mext*w2s_k);
    R2s = Lmat_h*w2s_k + N02s + 0.5*N1c1s           + 2*om_k*(Mext*w2c_k);

    Rphi = phi_s(1:2*nv)' * G * w1c_k(1:2*nv) + phi_c(1:2*nv)' * G * w1s_k(1:2*nv);

    R1c([bound; bound+nv]) = 0;
    R1s([bound; bound+nv]) = 0;
    R2c([bound; bound+nv]) = 0;
    R2s([bound; bound+nv]) = 0;

    % fprintf('norm(R1c bordo) = %.3e\n', norm(R1c([bound; bound+nv])));
    % fprintf('norm(R1s bordo) = %.3e\n', norm(R1s([bound; bound+nv])));
    % fprintf('norm(R2c bordo) = %.3e\n', norm(R2c([bound; bound+nv])));
    % fprintf('norm(R2s bordo) = %.3e\n', norm(R2s([bound; bound+nv])));
    % fprintf('norm(R0 bordo)  = %.3e\n', norm(R0([bound; bound+nv])));
    % fprintf('norm(w1c bordo) = %.3e\n', norm(w1c_k([bound; bound+nv])));
    % fprintf('norm(w1s bordo) = %.3e\n', norm(w1s_k([bound; bound+nv])));
end

function Nab = assemble_N(wa, wb, nv, np, xy, mv)
    Na  = navier_q2(xy, mv, wa);
    ub  = wb(1:nv);
    vb  = wb(nv+1:2*nv);
    Nab = [Na*ub; Na*vb; sparse(np,1)];
end

function Jbc = offdiagbc(J, nv, bound)
    Jt = J';
    Jt(:, bound) = 0;
    Jt(:, bound+nv) = 0;
    J = Jt';
    J(:, bound) = 0;
    J(:, bound+nv) = 0;
    Jbc = J;
end