%% Duffing方程非线性振动特性的计算与分析————王海波
%% 2.1Duffing的自由振动,图2.2
%% 奈奎斯特采样定理:只要采样频率超过信号最高频率的两倍，原始信号就可以从采样值中完美恢复，避免混叠效应
clear;clc;close all;
m = 1;c = 0;k = 1;beta = 0.1;   % 质量m、阻尼c、线性刚度k、非线性刚度beta
F0 = 0;omega = 0;                 % 无外力

% 定义Duffing方程的状态微分方程
duff = @(t, x) [x(2); (F0*cos(omega*t) - c*x(2) - k*x(1) - beta*x(1)^3)/m];
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);   % 设置求解精度（提高精度以减小长期积累误差）
tspan = linspace(0, 10000, 500000); % 模拟时间0-10000s
X0 = [1; 0];                         % 初始条件 [x(1); x(2)]，例如初始位移1，初始速度0
[t1, X1] = ode45(duff, tspan, X0, opts); % 时间区间 [0,100]

%% 绘图设置
subplot(2, 2, 1)
plot(t1, X1(:,1)); xlabel('Time (s)'); ylabel('Displacement (m)');
xlim([9950 10000]);ylim([min(X1(:,1)) max(X1(:,1))]);
title('Displacement-time');

subplot(2, 2, 2)
plot(t1, X1(:,2)); xlabel('Time (s)'); ylabel('Acceleration');
xlim([9950 10000]);
title('Acceleration-time');

subplot(2, 2, 3)
plot(X1(:,1), X1(:,2)); xlabel('Displacement (m)'); ylabel('Acceleration');
xlim([min(X1(:,1)) max(X1(:,1))]);
title('Phase trajectory');

% FFT频谱分析
x_resp = X1(7000:10000,1);            % 时域位移信号,样本要足够大
Fs = 1/mean(diff(t1));       % 采样频率(HZ)，diff(t2)计算时间点之间的间隔
N = length(x_resp);
f = (0:N-1)*(Fs/N);          % 除N是为了归一化        
Y = fft(x_resp);
subplot(2, 2, 4)
% 转换为单边谱的实际振幅
P2 = abs(Y)/N;               % 双边谱归一化
P1 = P2(1:floor(N/2)+1);     % 取前半部分（含Nyquist点）
P1(2:end-1) = 2*P1(2:end-1); % 非直流分量×2（得到实际振幅）
plot(f(1:floor(N/2)+1), P1);xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0 3]);
title('Spectrogram');
