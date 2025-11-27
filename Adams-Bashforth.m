%% 1. Khai báo Tham số và Nhập liệu

% Nhập giá trị h từ người dùng
h_str = input('Nhập giá trị bước thời gian h (ví dụ: 0.01): ', 's');
h = str2double(h_str);

if isnan(h) || h <= 0
    error('Giá trị h không hợp lệ. Vui lòng nhập một số dương.');
end

Us = 40;
R = 10;
L = 1;

T_final = 0.5;      % Thời gian mô phỏng cuối cùng (cố định)
N = round(T_final / h); % Số bước lặp tổng cộng

if N < 3
    error('Bước thời gian h quá lớn (h > T_final/3). Cần ít nhất 3 bước để khởi tạo AB3.');
end

t = 0:h:T_final;    % Vector thời gian

% Hàm vi phân f(i) = di/dt = 40 - 10*i
f = @(i) Us/L - R/L * i;

%% 2. Khởi tạo Giá trị Ban đầu và Giải tích

% Khởi tạo các vector kết quả
i_exact = zeros(1, N+1);
i_rk4 = zeros(1, N+1);
i_ab2 = zeros(1, N+1);
i_ab3 = zeros(1, N+1);

% Điều kiện ban đầu i(0) = 0
i_exact(1) = 0;
i_rk4(1) = 0;
i_ab2(1) = 0;
i_ab3(1) = 0;

% Hàm Giải tích (Chính xác): i(t) = 4 - 4*exp(-10*t)
i_exact = 4 - 4*exp(-10*t);

%% 3. Khởi tạo RK4 cho i1 và i2

% --- Tính i_1 (dùng RK4) ---
k1 = h * f(i_rk4(1));
k2 = h * f(i_rk4(1) + k1/2);
k3 = h * f(i_rk4(1) + k2/2);
k4 = h * f(i_rk4(1) + k3);
i_rk4(2) = i_rk4(1) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
i_ab2(2) = i_rk4(2);
i_ab3(2) = i_rk4(2);

% --- Tính i_2 (dùng RK4) ---
k1 = h * f(i_rk4(2));
k2 = h * f(i_rk4(2) + k1/2);
k3 = h * f(i_rk4(2) + k2/2);
k4 = h * f(i_rk4(2) + k3);
i_rk4(3) = i_rk4(2) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
i_ab2(3) = i_rk4(3);
i_ab3(3) = i_rk4(3);

%% 4. Lặp bằng Phương pháp Adams-Bashforth (AB2 và AB3)

% Bắt đầu lặp từ k=3 (tương đương n=2 trong công thức AB)
for k = 3:N
    % --- Công thức AB2 (Tính i_k+1) ---
    f_k = f(i_ab2(k));
    f_k_minus_1 = f(i_ab2(k-1));
    i_ab2(k+1) = i_ab2(k) + h * ( (3/2)*f_k - (1/2)*f_k_minus_1 );

    % --- Công thức AB3 (Tính i_k+1) ---
    f_k_ab3 = f(i_ab3(k));
    f_k_minus_1_ab3 = f(i_ab3(k-1));
    f_k_minus_2_ab3 = f(i_ab3(k-2));
    i_ab3(k+1) = i_ab3(k) + h * ( (23/12)*f_k_ab3 - (16/12)*f_k_minus_1_ab3 + (5/12)*f_k_minus_2_ab3 );
end

%% 5. Tính Sai số Tuyệt đối
Error_AB2 = abs(i_exact - i_ab2);
Error_AB3 = abs(i_exact - i_ab3);

%% 6. Xuất Bảng Giá trị Lặp (Giới hạn đến i5)
% Số lượng điểm cần hiển thị (0.05 / h) + 1. Đảm bảo T_final_table >= 0.05
T_final_table = 0.05;
N_table = round(T_final_table / h) + 1; 

fprintf('\n========================================================================\n');
fprintf('                BẢNG SO SÁNH GIÁ TRỊ DÒNG ĐIỆN i(t)    \n');
fprintf('  (Bước thời gian h = %.4f)  \n', h);
fprintf('  (Hiển thị giá trị từ i0 đến i5) \n');
fprintf('========================================================================\n');
fprintf('  t (s)  | i_ChinhXac (A) | i_AB2 (A) | i_AB3 (A) | Sai số AB2 | Sai số AB3 \n');
fprintf('------------------------------------------------------------------------\n');
% Chỉ lặp từ k=1 đến N_table
for k = 1:N_table
    fprintf(' %.4f |    %.6f    | %.6f  | %.6f  |  %.6f  |  %.6f \n', ...
        t(k), i_exact(k), i_ab2(k), i_ab3(k), Error_AB2(k), Error_AB3(k));
end
fprintf('========================================================================\n');

%% 7. Vẽ Đồ thị

% --- Đồ thị 1: So sánh các giá trị i(t) ---
figure(1);
% Đổi màu đường Giải tích thành xanh lá cây đậm ('g-')
plot(t, i_exact, 'g-', 'LineWidth', 2, 'DisplayName', 'Giải tích (Chính xác)');
hold on;
plot(t, i_ab2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Adams-Bashforth Bậc 2 (AB2)');
plot(t, i_ab3, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Adams-Bashforth Bậc 3 (AB3)');
hold off;

title(['So Sánh Dòng Điện i(t) (h = ', num2str(h), ')']);
xlabel('Thời gian t (s)');
ylabel('Dòng điện i (A)');
legend('Location', 'southeast');
grid on;

% --- Đồ thị 2: Đồ thị Sai số Tuyệt đối ---
figure(2);
plot(t, Error_AB2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Sai số AB2');
hold on;
plot(t, Error_AB3, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Sai số AB3');
hold off;

title(['Sai Số Tuyệt Đối (h = ', num2str(h), ')']);
xlabel('Thời gian t (s)');
ylabel('|i_{exact} - i_{approx}| (A)');
legend('Location', 'best');
grid on;

% --- Nhận xét tổng kết ---
fprintf('\n\n------------------------------------------------------------------------\n');
fprintf('NHẬN XÉT: Đánh giá độ chính xác với h = %.4f\n', h);
fprintf('------------------------------------------------------------------------\n');
fprintf('Sai số Tuyệt đối LỚN NHẤT (trên toàn bộ thời gian) của AB2: %.6f A\n', max(Error_AB2));
fprintf('Sai số Tuyệt đối LỚN NHẤT (trên toàn bộ thời gian) của AB3: %.6f A\n', max(Error_AB3));
fprintf('\nQuan sát đồ thị Sai số Tuyệt đối (Figure 2) để thấy rõ sự khác biệt:\n');
fprintf('- Phương pháp AB3 duy trì độ chính xác cao hơn AB2, thể hiện qua sai số tuyệt đối nhỏ hơn.\n');
fprintf('- Sai số của các phương pháp đa bước (AB2, AB3) phụ thuộc rất lớn vào độ chính xác của các giá trị khởi tạo (i1, i2).\n');
fprintf('------------------------------------------------------------------------\n');
