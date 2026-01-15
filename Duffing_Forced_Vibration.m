%% Duffing方程非线性振动特性的计算与分析————王海波
%% 2.2Duffing的受迫振动（主共振），图2.6（情形1）
%% 2.3Duffing的受迫振动的分数谐波共振,图2.11（情形2）
%% 2.4Duffing的受迫振动的混沌振动,图2.13（情形3）
%% 相比情形1将阻尼增大至5倍（情形4）
%% 最后绘制了正反扫频（幅频曲线图），情形1发生分岔跳跃现象，情形4正常
clear;clc;close all;
m = 1;c = 0.04;k = 1;beta = 0.15;F0 = 0.2;omega = 1.25;    % 情形1
%m = 1;c = 0.01;k = 1;beta = 0.2;F0 = 3;omega = 3.1;        % 情形2 
%m = 1;c = 0.05;k = 0.1;beta = 1.0;F0 = 8;omega = 1;       % 情形3
%m = 1;c = 0.20;k = 1;beta = 0.15;F0 = 0.2;omega = 1.25;        % 情形4                
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);  % 设置求解精度（提高精度以减小长期积累误差）
tspan = linspace(0, 10000, 50000);             % 模拟时间0-10000s
X0 = [2.9; 0];                                 % 假定初始位移和速度为零
[t2, X1] = ode45(@(t,x) duffing(t, x, m, c, k, beta, F0, omega), tspan, X0, opts);

%% 绘图设置
subplot(2, 2, 1)
plot(t2, X1(:,1)); xlabel('时间 (s)'); ylabel('位移 (m)');
xlim([900 1000]);ylim([-2.5 2.5]);
title('位移-时间图');

subplot(2, 2, 2)
plot(t2, X1(:,2)); xlabel('时间 (s)'); ylabel('加速度');
xlim([900 1000]);
title('加速度-时间图');

subplot(2, 2, 3)
plot(X1(900:1000,1), X1(900:1000,2)); xlabel('位移 (m)'); ylabel('加速度');
xlim([-2.5 2.5]);
title('相轨迹图');

% FFT频谱分析
x_resp = X1(7000:10000,1);            % 时域位移信号,样本要足够大
Fs = 1/mean(diff(t2));                % 采样频率(HZ)，diff(t2)计算时间点之间的间隔
N = length(x_resp);
f = (0:N-1)*(Fs/N);                   % 除N是为了归一化        
Y = fft(x_resp);
% 转换为单边谱的实际振幅
P2 = abs(Y)/N;                        % 双边谱归一化
P1 = P2(1:floor(N/2)+1);              % 取前半部分（含Nyquist点）
P1(2:end-1) = 2*P1(2:end-1);          % 非直流分量×2（得到实际振幅）
subplot(2, 2, 4)
plot(f(1:floor(N/2)+1), P1);xlabel('频率 (Hz)'); ylabel('幅值');
xlim([0 3]);
title('频谱图');

% 扫频
omegas = linspace(0.5, 1.5, 400);
amps_up = zeros(size(omegas));
amps_down = zeros(size(omegas));
tspan2 = linspace(0, 400, 2000);       

x_init = X0;
for i = 1:length(omegas)
    omega = omegas(i);
    [t, x] = ode45(@(t,x) duffing(t,x,m,c,k,beta,F0,omega), tspan2, x_init,opts);
    x_end = x(end-500:end,1);
    amps_up(i) = (max(x_end) - min(x_end))/2;
    x_init = [x(end,1); x(end,2)];
end

x_init = [x(end,1); x(end,2)];
for i = length(omegas):-1:1
    omega = omegas(i);
    [t, x] = ode45(@(t,x) duffing(t,x,m,c,k,beta,F0,omega), tspan2, x_init, opts);
    x_end = x(end-500:end,1);
    amps_down(i) = (max(x_end) - min(x_end))/2;
    x_init = [x(end,1); x(end,2)];
end

figure;
plot(omegas, amps_up, '-ro', 'DisplayName','正扫描');
hold on;
plot(omegas, amps_down, '-bs', 'DisplayName','反扫描');
xlabel('激励频率 ω'); ylabel('稳态振幅');
title('Duffing 幅频响应曲线'); legend;

function dx = duffing(t, x, m, c, k, beta, F0, omega)
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = (F0*cos(omega*t) - c*x(2) - k*x(1) - beta*x(1)^3)/m;
end