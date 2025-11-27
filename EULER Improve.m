% heun_rl_script.m
clc; clear; close all;

% ===== Nhap tham so =====
t0 = input('\nNhap t0 = ');
R  = input('Nhap R (Ohm) = ');
L  = input('Nhap L (H)   = ');
Vs = input('Nhap Vs (V)  = ');
i0 = input('Nhap i(t0) (A) = ');
tf = input('Nhap thoi gian ket thuc tf = ');
h  = input('Nhap buoc thoi gian h (vd: 0.01) = ');

% ===== Luoi thoi gian =====
t = t0:h:tf;
N = numel(t);
i = zeros(1,N);            % nghiem Heun
i(1) = i0;

% ===== Giai Euler cai tien (Heun) =====
for n = 1:N-1
    in = i(n);
    % Predictor (Euler)
    k1 = (Vs - R*in)/L;
    i_pred = in + h*k1;
    
    % Corrector
    k2 = (Vs - R*i_pred)/L;
    
    % Heun update
    i(n+1) = in + h*(k1 + k2)/2;
end

% ===== Nghiem dung va sai so =====
alpha = R/L;
iex = i0*exp(-alpha*(t - t0)) + (Vs/R)*(1 - exp(-alpha*(t - t0)));
err_abs = abs(iex - i);

% ===== In bang ket qua =====
fprintf('\nBANG KET QUA Euler-cai-tien (Heun) (h = %.4g)\n', h);
fprintf('-------------------------------------------------------------\n');
fprintf(' Thoi gian (s)   Heun (A)        Nghiem dung (A) |Sai so|\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:N
    fprintf('   %.4f          %.6f        %.6f        %.6f\n', ...
            t(k), i(k), iex(k), err_abs(k));
end
fprintf('-------------------------------------------------------------\n');

% ===== Ve do thi nghiem =====
figure('Name','i(t): Heun vs Exact');
plot(t, iex, 'b-', 'LineWidth', 1.6); hold on;
plot(t, i, 'ro--', 'LineWidth',
