function [flowsol, fst, gst, Jnst, Bst] = solve_navier_ifiss(viscosity, qmethod, maxit_p, maxit_n, tol_nl, spc, A, B, f, g, xy, mv, bound)
%SOLVE_NAVIER_IFISS Solves the 2D enclosed Navier-Stokes problem using IFISS
%
%   INPUT:
%       viscosity  - kinematic viscosity
%       qmethod    - finite element method: 0=Q1-P0, 1=Q1-Q1, >1=Q2-Q1
%       maxit_p    - maximum number of Picard iterations
%       maxit_n    - maximum number of Newton iterations
%       tol_nl     - nonlinear tolerance
%       spc        - streamline tracing type (1=uniform, 2=non-uniform)
%       A, B, f, g, xy, mv, bound    - mesh
%
%   OUTPUT:
%       flowsol    - solution vector [u; v; p]
%       xy, ev, mv - mesh data
%       bound      - boundary nodes
%       A, B       - finite element matrices

fprintf('Enclosed flow problem (IFISS)...\n')

%% Load assembled matrices
%gohome
%cd datafiles
%load square_stokes_nobc.mat   % loads A, B, f, g, xy, ev, mv, bound

%% Solve the Stokes problem
%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2; np=length(gst); 
if qmethod>1
   beta=0;
% xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
%----------------------------------- stabilized version
xstz=[Ast,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
       zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
elseif qmethod==1
	  beta=1/4;     % default parameter
%  xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
%----------------------------------- stabilized version
xstz=[Ast,Bst',zeros(2*nv,1);Bst,-beta*C,ones(np,1)/np; ...
      zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
elseif qmethod==0
   fprintf('computing pressure stabilized solution...\n')
   beta=1;
%  xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
%----------------------------------- stabilized version
xstz=[Ast,Bst',zeros(2*nv,1);Bst,-beta*C,ones(np,1)/np; ...
      zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
end

%% Compute the Stokes solution residual
if qmethod>1
   N = navier_q2(xy,mv,xst);
elseif qmethod<=1,
   nubeta=beta/viscosity;
   N = navier_q1(xy,ev,xst);
end
Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
if     qmethod>1, nlres = [Anst,Bst';Bst,sparse(np,np)]*xst-[fst;gst];
elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xst-[fst;gst];
end
nlres_norm  = norm(nlres);

%%% plot solution
%contourn = default('number of contour lines (default 50)',50);

%flowplot(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,spc,33);

%flowplot09(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,133)

pause(1) 
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
%
pde=4;
it_p = 0;
%

% nonlinear iteration 
%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p,
  %nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
   it_p = it_p+1;
   fprintf('\nPicard iteration number %g \n',it_p),
% compute Picard correction and update solution
   if     qmethod>1,
%  dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Anst,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
                zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   elseif qmethod<=1,
%  dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Anst,Bst',zeros(2*nv,1);Bst,-nubeta*C,ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,xns);
   end
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
% plot solution
   %%flowplot(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,spc,66); drawnow;
 %  flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,166);drawnow
   pause(1)
   flowsol = xns;
%% end of Picard iteration loop
end
%%
%
it_nl = it_p;
it_n = 0;
%% Newton iteration loop
while (nlres_norm > nlres0_norm*tol_nl) && (it_nl < maxit_p + maxit_n),
  %nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
   it_n = it_n+1;
   it_nl = it_nl+1;
   fprintf('\nNewton iteration number %g \n',it_n),
% compute Jacobian of current solution
   if     qmethod>1,  [Nxx,Nxy,Nyx,Nyy] = newton_q2(xy,mv,flowsol);
   elseif qmethod<=1, [Nxx,Nxy,Nyx,Nyy] = newton_q1(xy,ev,flowsol);
   end
   J = viscosity*A + [N + Nxx, Nxy; Nyx, N + Nyy];
   Jnst = newtonbc(J,xy,bound); 
% compute Newton correction and update solution
   if qmethod>1,
%  dxns = -[Jnst,Bst';Bst,sparse(np,np)]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Jnst,Bst',zeros(2*nv,1);Bst,sparse(np,np),ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   elseif qmethod<=1,
%  dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
%----------------------------------- stabilized version
   dxnsz= -[Jnst,Bst',zeros(2*nv,1);Bst,-nubeta*C,ones(np,1)/np; ...
            zeros(1,2*nv),ones(1,np)/np,zeros(1,1)]\[nlres;0];
   dxns=dxnsz(1:end-1); multiplier=dxnsz(end);
   end
   xns = flowsol + dxns;
% compute residual of new solution
   if     qmethod>1,  N = navier_q2(xy,mv,xns);
   elseif qmethod<=1, N = navier_q1(xy,ev,xns);
   end
   Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
   [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
   if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
   elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
   end
   nlres_norm = norm(nlres);
   nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
   fprintf('nonlinear residual is %e',nlres_norm)
   fprintf('\n   velocity change is %e\n',soldiff)
% plot solution
 %
 %   flowplot(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,spc,66); drawnow;
%flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,166);drawnow
   pause(1)
   flowsol = xns;
%% end of Newton iteration loop 
end
if nlres_norm <= nlres0_norm * tol_nl, 
   fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
% explicitly reorthogonalize to remove hydrostatic pressure component   
   nlres = nlres - [zeros(2*nv,1);(sum(nlres(2*nv+1:2*nv+np))/np)*ones(np,1)];
else
   fprintf('\nfinished, stopped on iteration counts\n\n');
end
%
%%% estimate errors
%if qmethod==1
 %  [jmpx,jmpy,els] = stressjmps_q1p0(viscosity,flowsol,xy,ev,ebound);
  % [error_x,error_y,fex,fey,ae] = navierpost_q1p0_p(viscosity,flowsol,jmpx,jmpy,els,xy,ev);
   %[error_x,error_y] = navierpost_q1p0_bc(viscosity,ae,fex,fey,...
    %                                    error_x,error_y,xy,ev,ebound);
   %error_div = q1div(xy,ev,flowsol);
   %errorest=sqrt(sum(error_x.^2 + error_y.^2 + error_div.^2));
   %fprintf('estimated overall error is %10.6e \n',errorest)
   %ee_error=sqrt((error_x.^2 + error_y.^2 + error_div.^2));
%% plot element errors
   %eplot(ee_error,ev,xy,x,y,67);
%% plot macroelement errors
   %mplot(ee_error,ev,xy,x,y,67);  title('Estimated error')
   %pause(5), figure(66)
%elseif qmethod==0
   %error_div = q1div(xy,ev,flowsol);	
%elseif qmethod>1, navierpost,

% Ensure Jnst is always defined for output
if ~exist('Jnst','var') || isempty(Jnst)
    fprintf('Computing Jnst from final solution (Newton not executed or converged early)...\n')
    [Nxx,Nxy,Nyx,Nyy] = newton_q2(xy,mv,flowsol);
    N = navier_q2(xy,mv,flowsol);
    J = viscosity*A + [N + Nxx, Nxy; Nyx, N + Nyy];
    Jnst = newtonbc(J,xy,bound);
end
end