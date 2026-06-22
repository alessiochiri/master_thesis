function plot_monitor_point(w0, w1c, w1s, w2c, w2s, omega, xy, nv)
% PLOT_MONITOR_POINT
% Reconstructs u and v time series at the Ghia et al. monitor point
%
% Cavity domain: [-1,1]² with moving lid at u_1 = +1 (top boundary)
% Ghia et al. reference point: (14/16, 13/16) in [0,1]²
% Mapped to [-1,1]² domain: (-0.75, 0.625)
% Symmetry transformation: u_1 changes sign, u_2 unchanged

% MONITOR POINT
x_mon = -0.75;
y_mon =  0.625;

dist = sqrt((xy(:,1) - x_mon).^2 + (xy(:,2) - y_mon).^2);
[~, idx] = min(dist);

% Extract velocity
u0_n  = w0(idx);    v0_n  = w0(idx+nv);
u1c_n = w1c(idx);   v1c_n = w1c(idx+nv);
u1s_n = w1s(idx);   v1s_n = w1s(idx+nv);
u2c_n = w2c(idx);   v2c_n = w2c(idx+nv);
u2s_n = w2s(idx);   v2s_n = w2s(idx+nv);

% Time reconstruction
T = 2*pi / abs(omega);
t = linspace(0, 2*T, 500);
t_norm = t / T;

% u_1 changes sign due to symmetry transformation
u = -(u0_n + u1c_n*cos(omega*t) - u1s_n*sin(omega*t) + u2c_n*cos(2*omega*t) - u2s_n*sin(2*omega*t));

% u_2 remains unchanged
v =   (v0_n + v1c_n*cos(omega*t) - v1s_n*sin(omega*t) + v2c_n*cos(2*omega*t) - v2s_n*sin(2*omega*t));

% Harmonic amplitudes
u1_mean = -u0_n;
u1_amp_1H = sqrt(u1c_n^2 + u1s_n^2);
u1_amp_2H = sqrt(u2c_n^2 + u2s_n^2);

u2_mean = v0_n;
u2_amp_1H = sqrt(v1c_n^2 + v1s_n^2);
u2_amp_2H = sqrt(v2c_n^2 + v2s_n^2);

%% --- Plot time series ---
figure('Name', 'Monitor Point Time Series', 'Position', [100 100 1000 700]);

% u - Streamwise velocity
subplot(2,1,1);
plot(t_norm, u, 'b-', 'LineWidth', 2);
hold on
yline(u1_mean, 'k--', 'LineWidth', 1, 'Label', 'Mean flow');
hold off
xlabel('t / T', 'FontSize', 11);
ylabel('u_1 [m/s]', 'FontSize', 11);
title(sprintf('Streamwise velocity u(t) — Monitor point (%.4f, %.4f)', xy(idx,1), xy(idx,2)), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('HB2 solution', 'Mean flow', 'FontSize', 10, 'Location', 'best');
xlim([0 2]);

% v - Transverse velocity
subplot(2,1,2);
plot(t_norm, v, 'r-', 'LineWidth', 2);
hold on
yline(u2_mean, 'k--', 'LineWidth', 1, 'Label', 'Mean flow');
hold off
xlabel('t / T', 'FontSize', 11);
ylabel('u_2 [m/s]', 'FontSize', 11);
title(sprintf('Transverse velocity v(t) — Monitor point (%.4f, %.4f)', xy(idx,1), xy(idx,2)), ...
    'FontSize', 12, 'FontWeight', 'bold');
grid on;
legend('HB2 solution', 'Mean flow', 'FontSize', 10, 'Location', 'best');
xlim([0 2]);

% Print statistics
fprintf('\n HARMONIC COEFFICIENTS AND AMPLITUDES\n');

fprintf('\n--- Horizontal velocity (u_1) ---\n');
fprintf('Mean (u_1)^0            : %.6f m/s\n', u1_mean);
fprintf('Amplitude 1st harmonic  : %.6e m/s\n', u1_amp_1H);
fprintf('Amplitude 2nd harmonic  : %.6e m/s\n', u1_amp_2H);
fprintf('Ratio A_2 / A_1         : %.4f%%\n', 100*u1_amp_2H/(u1_amp_1H+1e-16));

fprintf('\n--- Vertical velocity (u_2) ---\n');
fprintf('Mean (u_2)^0            : %.6f m/s\n', u2_mean);
fprintf('Amplitude 1st harmonic  : %.6e m/s\n', u2_amp_1H);
fprintf('Amplitude 2nd harmonic  : %.6e m/s\n', u2_amp_2H);
fprintf('Ratio A_2 / A_1         : %.4f%%\n', 100*u2_amp_2H/(u2_amp_1H+1e-16));

fprintf('\n--- Periodicity ---\n');
fprintf('Period T                : %.6f s\n', T);
fprintf('Angular frequency ω     : %.6f rad/s\n', abs(omega));
fprintf('Frequency f             : %.6f Hz\n', abs(omega)/(2*pi));
fprintf('=======================================\n\n');

end