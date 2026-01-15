%% 基于MATLAB的强迫振动达芬方程的非线性幅频响应分析————陈赵江
%% 强迫振动达芬方程幅频响应曲线的解析求解
%% 解决了跳跃现象  
clear; clc; close all;
% 参数设定 
epsilon=0.2;      % 小参数
gamma1 = 2;       % 非线性参数
F0 = 0.5;         % 激励幅值
zeta1 = 0.25;     % 阻尼项
a = linspace(0.1, 2.0, 100);  % 计算响应幅值范围

    figure('Position', [100, 100, 900, 700])
    set(gcf, 'Color', 'w')
for ii = 1:1:length(a)
    % 根据(23)式,对应一个响应幅值,可能存在两个激励频率值;
    % Omega1为较小的激励频率值
    Omega1(ii) = (1 + 3*epsilon*gamma1/8*a(ii)^2 - sqrt((epsilon*F0/(2*a(ii)))^2 - (epsilon*zeta1)^2));
    
    % Omega1对应的本征函数值
    lamOmega11 = sqrt(-((Omega1(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega1(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2)) - zeta1;
    lamOmega12 = -sqrt(-((Omega1(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega1(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2)) - zeta1;
    
    % Omega2为较大的激励频率值
    Omega2(ii) = (1 + 3*epsilon*gamma1/8*a(ii)^2 + sqrt((epsilon*F0/(2*a(ii)))^2 - (epsilon*zeta1)^2));
    
    % Omega2对应的本征函数值
    lamOmega21 = sqrt(-((Omega2(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega2(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2)) - zeta1;
    lamOmega22 = -sqrt(-((Omega2(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega2(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2)) - zeta1;

    if lamOmega11 == conj(lamOmega21)
        % 如果两个本征值相等退出计算(即此时计算到频响曲线最高点)
        plot(Omega1(ii), a(ii), '-bo')
        hold on
        break
    else
        FG = ((Omega1(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega1(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2) + zeta1^2; % 对应上文中的(31)式
        if FG < 0
            plot(Omega1(ii), a(ii),'-bo') % 不稳定点
        else
            plot(Omega1(ii), a(ii), '-bo') % 稳定点
            hold on
        end
        
        GF = ((Omega2(ii)-1)/epsilon - 3*gamma1/8*a(ii)^2)*((Omega2(ii)-1)/epsilon - 9*gamma1/8*a(ii)^2) + zeta1^2; % 对应上文中的(31)式
        if GF < 0
            plot(Omega2(ii), a(ii), '-ro') % 不稳定点
        else
            plot(Omega2(ii), a(ii), '-bo') % 稳定点
            hold on
        end
        
xlabel('激励频率 ω'); ylabel('稳态振幅');
        xlim([0.5, 1.5]); % x轴范围
        ylim([0, 1.2]);     % y轴范围
        title('Duffing 幅频响应曲线'); legend;
    end % 满足条件退出循环(计算到幅频响应曲线顶点)
end % 结束对不同a值的计算