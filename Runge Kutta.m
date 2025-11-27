% rk4_rl_script.m
clc; clear; close all;

% ===== Nhap tham so =====
t0 = input('\nNhap t0 = ');
R  = input('Nhap R (Ohm) = ');
L  = input('Nhap L (H)   = ');
Vs = input('Nhap Vs (V)  = ');
i0 = input('Nhap i(t0) (A) = ');
tf = input('Nhap thoi gian ket thuc tf = ');
h  = input('Nhap buoc thoi gian h (vd: 0.1) = ');

% ===== Luoi thoi gian =====
t = t0:h:tf;                      % day thoi gian
N = numel(t);
i = zeros(1,N);                   % nghiem RK4
i(1) = i0;

% ===== Giai RK4 (di/dt = (Vs - R*i)/L) =====
for n = 1:N-1
    in = i(n);

    k1 = h * (Vs - R*in)/L;
    k2 = h * (Vs - R*(in + k1/2))/L;
    k3 = h * (Vs - R*(in + k2/2))/L;
    k4 = h * (Vs - R*(in + k3))/L;

    i(n+1) = in + (k1 + 2*k2 + 2*k3 + k4)/6;
end

% ===== Nghiem dung va sai so =====
alpha = R/L;
iex = i0*exp(-alpha*(t - t0)) + (Vs/R)*(1 - exp(-alpha*(t - t0)));
err_abs = abs(iex - i);
err_rel = err_abs ./ max(abs(iex), eps);

% ===== In bang ket qua =====
fprintf('\nBANG KET QUA Runge-Kutta-4 (h = %.4g)\n', h);
fprintf('--------------------------------------------------------\n');
fprintf(' Thoi gian (s)\tRunge-Kutta-4 (A)\t\tNghiem dung (A)\t|Sai so|\n');
fprintf('--------------------------------------------------------\n');
for k = 1:N
    fprintf('   %.4f\t\t%.6f\t%.6f\t%.6f\n', t(k), i(k), iex(k), err_abs(k));
end
fprintf('--------------------------------------------------------\n');

% ===== Ve do thi nghiem =====

figure('Name','i(t): Runge-Kutta-4 vs Exact');
plot(t, iex, 'b-', 'LineWidth', 1.6); hold on;     % Blue nghiệm đúng
plot(t, i, 'ro--', 'LineWidth', 1.2, 'MarkerSize', 5); % Red RK4
grid on; xlabel('t (s)'); ylabel('i (A)');
title('So sanh Ngiem Runge-Kutta-4 va Nghiem dung');
legend('Nghiem dung','RK4','Location','best');
% ===== Ve do thi sai so =====
figure('Name','Sai so tuyet doi');
plot(t, err_abs, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on; xlabel('t (s)'); ylabel('|sai so| (A)');
title('Sai so tuyet doi');


figure('Name','Sai so tuong doi');
plot(t, err_rel, '-o', 'LineWidth', 1.2, 'MarkerSize', 5);
grid on; xlabel('t (s)'); ylabel('Sai so tuong doi');
title('Sai so tuong doi');
