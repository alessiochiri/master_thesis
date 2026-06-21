function [real_hopf, imag_hopf, phi_c, phi_s] = eigenvalvec(xy, bound, Jnst, G, Bst, gst)

np = length(gst);
Gst = mass_q2_bc(G, xy, bound);
delta = -1e-3;

Mfull = [-Gst,        delta*Bst';
          delta*Bst,  sparse(np,np)];
Jfull = [Jnst, Bst'; Bst, sparse(np,np)];

k = 350;
tic;
[V, Dr] = eigs(Jfull, Mfull, k, 'smallestabs');
fprintf('Tempo di calcolo eigs: %.2f secondi\n', toc);

lambda = diag(Dr);

% seleziona autovalore complesso con parte reale massima
idx_hopf = abs(imag(lambda)) > 1e-6;
lambda_hopf = lambda(idx_hopf);
V_hopf = V(:, idx_hopf);

[~, idx] = max(real(lambda_hopf));
lambda_crit = lambda_hopf(idx);
phi_complex = V_hopf(:, idx);

% separa modi coseno e seno con normalizzazione coerente
scale = norm(phi_complex);
phi_c = real(phi_complex) / scale;
phi_s = -imag(phi_complex) / scale;

fprintf('Autovalore critico: %.6f + %.6fi\n', real(lambda_crit), imag(lambda_crit));

if real(lambda_crit) < 0
    disp('Flusso stabile');
else
    disp('Flusso instabile - Biforcazione di Hopf');
end

real_hopf = real(lambda_crit);
imag_hopf = imag(lambda_crit);

end