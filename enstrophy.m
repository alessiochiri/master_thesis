function EN =enstrophy(flowsol,By,Bx,G,xy)
% Compute Enstrophy of the flow
%
%   input
%          flowsol    flow solution vector
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          G          scalar velocity mass matrix
%          xy         velocity nodal coordinate vector  
%   output
%          EN          enstrophy vector   
nvtx = length(xy); 
nu = 2*nvtx;
u = flowsol(1:nu);
f = [-By, Bx] * u;
w = G(1:nvtx, 1:nvtx) \ f;
EN = 0.5 * w' * G(1:nvtx, 1:nvtx) * w;
end


