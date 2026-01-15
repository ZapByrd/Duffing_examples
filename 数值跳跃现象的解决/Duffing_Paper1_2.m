%% 基于MATLAB的强迫振动达芬方程的非线性幅频响应分析
%% 强迫振动达芬方程幅频响应曲线的数值求解  
clear; clc; close all;
global epsilon gamma1 F1 zeta1 Omega  % 定义全局变量

% 参数设定
epsilon=0.2;      % 小参数
gamma1 = 2;       % 非线性参数
F1 = 0.5;         % 激励幅值
zeta1 = 0.25;     % 阻尼项
np = 400;         % 计算的频率数目

Omega1 = linspace(0.5, 1.5, np);  % 正向扫频时的频率
Omega2 = linspace(1.5, 0.5, np);  % 反向扫频时的频率

yy = [];
Y0 = [0.1 0.1];  % 初始值

% 正向扫频
for i = 1:1:length(Omega1)
    Omega = Omega1(i);
    [T, Y] = ode45(@duffing, [0 400], Y0);
    nn = length(Y(:,1));
    ymax = max(Y(nn-round(nn/2):nn, 1));  % 取稳态幅值
    ymax1 = max(Y(nn-round(nn/2):nn, 2));
    Y0 = [ymax ymax1];  % 下一次计算的初始值
    yy = [yy; ymax];    % 保存计算结果到yy向量
end
figure('Position', [100, 100, 900, 700])
set(gcf, 'Color', 'w')
plot(Omega1, yy, '-ro','DisplayName','正扫描')  % 正向扫频曲线用红色表示
hold on

% 反向扫频
yy1 = [];
ymax1 = max(Y(nn-round(nn/2):nn, 1));  % 反向扫频的初始值
ymax2 = max(Y(nn-round(nn/2):nn, 2));  % 反向扫频的初始值

for j = 1:1:length(Omega2)
    Omega = Omega2(j);
    [T, Y1] = ode45(@duffing, [0 400], [ymax1 ymax2]);
    nn1 = length(Y1(:,1));
    ymax1 = max(Y1(nn1-round(nn1/2):nn1, 1));
    ymax2 = max(Y1(nn1-round(nn1/2):nn1, 2));
    yy1 = [yy1; ymax1];
end

plot(Omega2, yy1, '-bo', 'DisplayName','反扫描')  % 反向扫频曲线用蓝色表示

xlabel('激励频率 ω'); ylabel('稳态振幅');
xlim([0.5, 1.5]); % x轴范围
ylim([0, 1.2]);     % y轴范围
title('Duffing 幅频响应曲线'); legend;

function dy = duffing(t, y)
    global epsilon gamma1 F1 zeta1 Omega
    dy = zeros(2, 1);
    dy(1) = -zeta1*y(1) + F1/2*sin(y(2));
    dy(2) = (Omega-1)/epsilon - 3*gamma1/8*y(1)^2 + F1/(2*y(1))*cos(y(2));
end