function Mnst = mass_q2_bc(M, xy, bound)

% MASS_Q2_BC Apply Dirichlet BC to Q2 mass matrix
%
%   Mnst = mass_q2_bc(M, xy, bound)
%
%   M      : 2*nvtx x 2*nvtx mass matrix
%   xy     : nodal coordinates [nvtx x 2]
%   bound  : vector of boundary node indices
%   Mnst   : mass matrix with Dirichlet BC applied (same size as M)

nvtx = size(xy,1);
nbd = length(bound);

% build diagonal vector for constrained DOFs
dA = zeros(nvtx,1);
dA(bound) = ones(nbd,1);

% prepare blocks for u and v velocity components
Iu = 1:nvtx;
Iv = nvtx+1:2*nvtx;

% copy original matrix
Mnst = M;

% zero out rows and columns corresponding to constrained DOFs
null_col = sparse(2*nvtx, nbd);

% u block
Mnst(:,bound) = null_col;
Mnst(bound,:) = null_col';

% v block
Mnst(:,nvtx+bound) = null_col;
Mnst(nvtx+bound,:) = null_col';

% add unit diagonal on constrained DOFs
Mnst = Mnst + [spdiags(dA,0,nvtx,nvtx), sparse(nvtx,nvtx);
               sparse(nvtx,nvtx), spdiags(dA,0,nvtx,nvtx)];
end