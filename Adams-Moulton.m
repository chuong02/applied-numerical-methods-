%% Adam-Moulton Bậc 2 (AM2) cho PTPV Mạch RL - Phiên bản cuối cùng
% Giai PTPV: di/dt = Us/L - (R/L)*i
% Khoi dong bang Runge-Kutta Bậc 4 (RK4)

clc;        % Xóa cửa sổ Command
clear all;  % Xóa tất cả các biến

disp('==================================================');
disp('   GIAI PTPV MACH RL BANG PP ADAM-MOULTON BAC 2    ');
disp('==================================================');

% --- 1. Nhập Tham số từ người dùng ---
t0 = input('Nhap t0 (thoi gian bat dau) = ');
R = input('Nhap R (Ohm) = ');
L = input('Nhap L (H) = ');
Us = input('Nhap Us (V) = '); % ĐÃ SỬA TÊN BIẾN TỪ Vs -> Us
i0 = input('Nhap i(t0) (A) = ');
tf = input('Nhap thoi gian ket thuc tf = ');
h = input('Nhap buoc thoi gian h (vd: 0.01) = ');

% --- 2. Thiết lập Động ---
tau_inv = R / L;          
steady_state = Us / R;    % ĐÃ SỬA: Us / R
C = (i0 - steady_state) * exp(tau_inv * t0); 

% Định nghĩa hàm f(i, t) và nghiệm chính xác
f = @(i_val) steady_state * tau_inv - tau_inv * i_val;
i_exact_func = @(t_val) steady_state + C * exp(-tau_inv * t_val);

t = t0:h:tf;
N = length(t);
i_approx_am2 = zeros(1, N);  
i_approx_ab2 = zeros(1, N);  
f_val = zeros(1, N);      

i_approx_am2(1) = i0;
i_approx_ab2(1) = i0;
f_val(1) = f(i0);

if N < 3
    error('So buoc qua it. Vui long tang tf hoac giam h.');
end

% --- 3. Khởi động bằng RK4 (Tính i(2) và i(3)) ---
for n = 1:2
    k1 = h * f(i_approx_am2(n));
    k2 = h * f(i_approx_am2(n) + k1/2);
    k3 = h * f(i_approx_am2(n) + k2/2);
    k4 = h * f(i_approx_am2(n) + k3);
    i_approx_am2(n+1) = i_approx_am2(n) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    i_approx_ab2(n+1) = i_approx_am2(n+1); 
    f_val(n+1) = f(i_approx_am2(n+1));
end

%% --- 4. Giải bằng Adam-Moulton Bậc 2 (Dự đoán-Hiệu chỉnh) ---

for n = 3:N-1
    f_n_minus_1 = f_val(n-1);
    f_n = f_val(n);
    
    % A. Dự đoán (Predictor - AB2)
    i_star = i_approx_am2(n) + (h/2) * (3*f_n - f_n_minus_1);
    i_approx_ab2(n+1) = i_star; % LƯU giá trị dự đoán AB2 (i*)
    
    % B. Tính f dự đoán
    f_star = f(i_star);
    
    % C. Hiệu chỉnh (Corrector - AM2)
    i_approx_am2(n+1) = i_approx_am2(n) + (h/2) * (f_n + f_star);
    f_val(n+1) = f(i_approx_am2(n+1));
end

% --- 5. Trình bày Bảng Giá trị Chi tiết ---

i_exact = i_exact_func(t);

% Tính Sai số Tuyệt đối (Sử dụng cho cả Bảng và Đồ thị)
Error_AM2 = abs(i_approx_am2 - i_exact);
Error_AB2 = abs(i_approx_ab2 - i_exact); 

fprintf('\n==========================================================================================================\n');
fprintf('                             BẢNG SO SÁNH GIÁ TRỊ VÀ SAI SỐ (h = %.4f)\n', h);
fprintf('==========================================================================================================\n');
fprintf('| %-6s | %-12s | %-12s | %-12s | %-12s | %-12s |\n', ...
    't (s)', 'i_exact', 'i_AM2', 'i_AB2*', '|AM2-i_chính xác|', '|AB2-i_chính xác|'); 
fprintf('----------------------------------------------------------------------------------------------------------\n');

for k = 1:N
    % Định dạng chuỗi cho giá trị i_AB2 (Hiển thị RK4 cho 2 điểm đầu)
    if k <= 2
        i_ab2_str = '   RK4   ';
    else
        i_ab2_str = sprintf('%12.6f', i_approx_ab2(k));
    end
    
    % In hàng dữ liệu
    fprintf('| %-6.2f | %-12.6f | %-12.6f | %-12s | %-12.6f | %-12.6f |\n', ...
        t(k), i_exact(k), i_approx_am2(k), i_ab2_str, Error_AM2(k), Error_AB2(k));
end

fprintf('==========================================================================================================\n');

% --- 6. Vẽ Đồ thị ---

figure('Position', [100, 100, 1000, 400]);

% Đồ thị 1: So sánh Dòng điện (Giữ nguyên)
subplot(1, 2, 1);
plot(t, i_exact, 'b-', 'LineWidth', 2, 'DisplayName', 'i_{chính xác}');
hold on;
plot(t, i_approx_am2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'i_{AM2}');
plot(t, i_approx_ab2, 'g:', 'LineWidth', 1.5, 'DisplayName', 'i_{AB2}*');
title('So sánh Dòng điện');
xlabel('Thời gian t (s)');
ylabel('Dòng điện i(t) (A)');
legend('show', 'Location', 'southeast');
grid on;

% Đồ thị 2: Sai số Tuyệt đối (Giữ nguyên)
subplot(1, 2, 2);
plot(t, Error_AM2, 'r-', 'LineWidth', 2, 'DisplayName', '|AM2 - i_{chính xác}|');
hold on;
plot(t, Error_AB2, 'b--', 'LineWidth', 1.5, 'DisplayName', '|AB2* - i_{chính xác}|'); 
title('Sai số Tuyệt đối');
xlabel('Thời gian t (s)');
ylabel('Sai số Tuyệt đối');
legend('show');
grid on;
