% Bifurcation Identification via bisection on Re(lambda*)
load('square_stokes_nobc.mat', 'A','B','Bx','By','G','f','g','xy','mv','bound');
Re_a = 7800; 
Re_b = 8200;  

[f_a, om_a] = eval_rightmost_eig(Re_a, A, B, G, f, g, xy, mv, bound);
[f_b, om_b] = eval_rightmost_eig(Re_b, A, B, G, f, g, xy, mv, bound);
fprintf('Re=%.1f: Re(lambda*)=%.4e\n', Re_a, f_a);
fprintf('Re=%.1f: Re(lambda*)=%.4e\n', Re_b, f_b);

if sign(f_a) == sign(f_b)
    error('Same sign at the extrama: increase range');
end

tol_bisect = 1e-2;
max_iter   = 10;

Re_hist = zeros(max_iter,1);
f_hist  = zeros(max_iter,1);
om_hist = zeros(max_iter,1);

fprintf('%-6s %-10s %-14s %-10s\n', 'Iter', 'Re', 'Re(lambda*)', 'Im(lambda*)');

for k = 1:max_iter
    Re_mid = 0.5*(Re_a + Re_b);
    [f_mid, om_mid] = eval_rightmost_eig(Re_mid, A, B, G, f, g, xy, mv, bound);

    Re_hist(k) = Re_mid;
    f_hist(k)  = f_mid;
    om_hist(k) = om_mid;

    fprintf('%-6d %-10.4f %-14.4e %-10.4f\n', k, Re_mid, f_mid, om_mid);

    if sign(f_mid) == sign(f_a)
        Re_a = Re_mid; f_a = f_mid;
    else
        Re_b = Re_mid; f_b = f_mid;
    end

    if abs(Re_b - Re_a) < tol_bisect
        break
    end
end

Re_c   = 0.5*(Re_a + Re_b);
om_c   = om_mid;
n_iter = k;

fprintf('\n=== Re_c = %.4f  (Im(lambda*) at Re_c = %.4f) ===\n', Re_c, om_c);

Re_hist = Re_hist(1:n_iter);
f_hist  = f_hist(1:n_iter);
om_hist = om_hist(1:n_iter);

save('bisection_Rec.mat', 'Re_hist', 'f_hist', 'om_hist', 'Re_c', 'om_c');

% --- Plot ---
figure('Position',[100 100 900 400]);
subplot(1,2,1)
plot(Re_hist, f_hist, 'bo-', 'LineWidth',1.5, 'MarkerFaceColor','b'); hold on
yline(0,'k--'); xline(Re_c,'r--','LineWidth',1.5);
xlabel('Re'); ylabel('Re(\lambda^*)','Interpreter','tex')
title('Bisection convergence'); grid on; box on;

subplot(1,2,2)
plot(1:n_iter, Re_hist, 'bo-','LineWidth',1.5); hold on
yline(Re_c,'r--','LineWidth',1.5);
xlabel('Iteration'); ylabel('Re')
title('Re_k \rightarrow Re_c'); grid on; box on;

function [f_val, om_val] = eval_rightmost_eig(Re, A, B, G, f, g, xy, mv, bound)
    viscosity = 2/Re;
    [flowsol, fst, gst, Jnst, Bst] = solve_navier_ifiss(viscosity, 3, 20, 20, 1e-8, 2, A, B, f, g, xy, mv, bound);

    [Nxx, Nxy, Nyx, Nyy] = newton_q2(xy, mv, flowsol);
    Nst  = navier_q2(xy, mv, flowsol);
    J    = viscosity*A + [Nst+Nxx, Nxy; Nyx, Nst+Nyy];
    Jnst_full = newtonbc(J, xy, bound);

    [real_hopf, imag_hopf] = eigenvalueproblem(xy, bound, Jnst_full, G, Bst, gst);

    f_val  = real_hopf;
    om_val = imag_hopf;
end