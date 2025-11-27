%EULER

clc; clear; close all;


R1 = input('Nhập giá trị R (Ohm): ');
h = input('Nhập bước thời gian h (s): ');


L1 = 1;
E1 = 440;
T = 0.85;          
t = 0:h:T;
N = length(t);


i_nghiem = (E1/R1) * (1 - exp(-R1/L1 * t));

% Euler
i_euler = zeros(1,N);
i_euler(1) = 0;   % điều kiện đầu
for k = 1:N-1
    di = (E1 - R1*i_euler(k)) / L1;
    i_euler(k+1) = i_euler(k) + h * di;
end


error_percent = abs((i_nghiem - i_euler) ./ i_nghiem) * 100;
error_percent(i_nghiem == 0) = 0; % tránh chia cho 0 ở t=0

% In bảng kết quả
fprintf('-------------------------------------------------------------\n');
fprintf(' Thoi gian (s)     Euler (A)     PTVP (A)     Sai so (%%)\n');
fprintf('-------------------------------------------------------------\n');
for k = 1:N
    fprintf('   %6.4f         %7.4f      %7.4f       %7.3f\n', ...
        t(k), i_euler(k), i_nghiem(k), error_percent(k));
end
fprintf('-------------------------------------------------------------\n');

% Tính sai số trung bình
mean_error = mean(error_percent(2:end)); % bỏ t=0
fprintf('Sai số trung bình: %.3f%%\n', mean_error);

% Vẽ đồ thị so sánh
figure;
subplot(2,1,1);
plot(t, i_nghiem, 'r-', 'LineWidth', 2); hold on;
plot(t, i_euler, 'bo--', 'LineWidth', 1);
legend('PTVP (Chính xác)', 'Euler', 'Location', 'southeast');
xlabel('Thời gian t (s)');
ylabel('Dòng điện i (A)');
title('So sánh nghiệm i(t) giữa Euler và nghiệm chính xác');
grid on;

% Vẽ đồ thị sai số phần trăm
subplot(2,1,2);
plot(t, error_percent, 'm-', 'LineWidth', 2);
xlabel('Thời gian t (s)');
ylabel('Sai số (%)');
title('Đồ thị sai số phần trăm giữa Euler và nghiệm chính xác');
grid on;
