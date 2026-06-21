function profiles = centerline_profiles(xy, flowsol, fst, gst)
% Extract velocity profiles on cavity centerlines (2x2 cavity)
%
% INPUT:
%   xy       : node coordinates [nvtx x 2]
%   flowsol  : solution vector [u; v; p]
%   fst      : velocity rhs vector (for nv)
%   gst      : pressure rhs vector (for np)
%
% OUTPUT:
%   profiles : struct containing centerline profiles and coordinates

nv = length(fst) / 2;
np = length(gst);
u = flowsol(1:nv);
v = flowsol(nv+1:2*nv);
p = flowsol(2*nv+1:2*nv+np);

x = xy(:, 1);
y = xy(:, 2);

x_center = (max(x) + min(x)) / 2;
y_center = (max(y) + min(y)) / 2;

tol = 1e-6;

%% Vertical line (x = x_center)
idx_vert = find(abs(x - x_center) < tol);
if isempty(idx_vert)
    tol = 1e-3;
    idx_vert = find(abs(x - x_center) < tol);
end

y_vert = y(idx_vert);
u_vert = u(idx_vert);

[y_vert, sort_idx] = sort(y_vert);
u_vert = u_vert(sort_idx);

%% Horizontal line (y = y_center)
idx_horiz = find(abs(y - y_center) < tol);
if isempty(idx_horiz)
    tol = 1e-3;
    idx_horiz = find(abs(y - y_center) < tol);
end

x_horiz = x(idx_horiz);
v_horiz = v(idx_horiz);

[x_horiz, sort_idx] = sort(x_horiz);
v_horiz = v_horiz(sort_idx);

%% Output struct
profiles.y_centerline = y_vert;
profiles.u_centerline = u_vert;
profiles.x_centerline = x_horiz;
profiles.v_centerline = v_horiz;
profiles.x_center = x_center;
profiles.y_center = y_center;
