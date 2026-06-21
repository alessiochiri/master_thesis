function [real_hopf, imag_hopf] = eigenvalueproblem(xy, bound, Jnst, G, Bst, gst)

% Solve the generalized eigenvalue problem Jnstx=Gx
%
% INPUT:
%   xy       : node coordinates [nvtx x 2]
%   bound    : node coordinates [boundary]
%   gst      : pressure rhs vector (for np)
%   Jnst     : Jacobian Navier Stokes
%   Bst      : diffusion matrix
%   G        : mass matrix
%
% OUTPUT:
% real_hopf     : real part max eigenvalue
% imag_hopf     : imaginary part max eigenvalue

np = length(gst);
Gst = mass_q2_bc(G, xy, bound);

delta = -1e-3;
Mfull = [-Gst,       delta*Bst';
          delta*Bst, sparse(np,np)];
Jfull = [Jnst, Bst'; Bst, sparse(np,np)];

k = 350;
tic;
[~, Dr] = eigs(Jfull, Mfull, k, 'smallestabs');
fprintf('Tempo di calcolo eigs: %.2f secondi\n', toc);
lambda = diag(Dr);

idx_hopf = abs(imag(lambda)) > 1e-6;
lambda_hopf = lambda(idx_hopf);
[~, idx] = max(real(lambda_hopf));
lambda_crit = lambda_hopf(idx);

fprintf('%.5f + %.5fi\n', real(lambda_crit), imag(lambda_crit));

if real(lambda_crit) < 0
    disp('Stable Flow');
elseif real(lambda_crit) > 0
    disp('Unstable Flow');
end

real_hopf = real(lambda_crit);
imag_hopf = imag(lambda_crit);