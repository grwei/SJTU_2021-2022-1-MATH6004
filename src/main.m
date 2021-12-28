%% 计算方法实验题
% Author: 危国锐(313017602@qq.com)
% Created: 2021-12-28
% Last modified:

%% 配置环境

clc;clear;close all;

%% Q1

x = -1:1e-2:1;
Q1_L_10 = zeros(size(x));
Q1_x_nodes = -1:2/10:1;
for k = 0:10
    Q1_L_10 = Q1_L_10 + 1./(1 + 25 * Q1_x_nodes(k+1).^2) .* Q1_interp_base(x,k,10);
end
% 插值型积分系数
Q1_int_coef = zeros(size(Q1_x_nodes));
for i = 0:length(Q1_int_coef)-1
    Q1_int_coef(i+1) = integral(@(x) Q1_interp_base(x,i,10),-1,1);
end
% 检查代数精度
Q1_I10_f = sum(Q1_int_coef ./ (1 + 25 * Q1_x_nodes.^2)); %
Q1_I_f = atan(5) * 2 / 5;
Q1_I10_x11 = sum(Q1_int_coef .* Q1_x_nodes.^11); % I_x11 = \int_{-1}^{1} x^11 dx = 0
Q1_I10_x12 = sum(Q1_int_coef .* Q1_x_nodes.^12); % I_x12 = \int_{-1}^{1} x^12 dx = 2/13

%% Q2

x = -1:1e-2:1;
Q2_L_10 = zeros(size(x));
Q2_x_nodes = cos(pi * (2 * (0:10) + 1) ./ (2*10 + 2));
for k = 0:10
    Q2_L_10 = Q2_L_10 + 1./(1 + 25 * Q2_x_nodes(k+1).^2) .* Q2_interp_base(x,k,10);
end

%% Q3

x = -1:1e-2:1;
Q3_L_10 = zeros(size(x));
Q3_x_nodes = -1:2/10:1;
for k = 0:10
    Q3_L_10 = Q3_L_10 + 1./(1 + 25 * Q3_x_nodes(k+1).^2) .* Q3_interp_base(x,k,10);
end

%% 绘图

%% Q1
figure('Name','Q1')
h = plot(x,[1./(1 + 25 * x.^2);Q1_L_10]);
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel("\fontname{Times New Roman} \fontsize{10} \it x")
ylabel("\fontname{Times New Roman} \fontsize{10} \it y")
legend( h, '$$f(x) = \frac{1}{1 + 25 x^2}$$','$$L_{10}^{(1)}(x) = \sum_{k = 0}^{10} f(x_k) l_k^{(1)}(x)$$','Location', 'north', 'Interpreter', 'latex','fontsize',10 );
legend('boxoff')
title(sprintf("\\fontname{Times New Roman} \\fontsize{10} \\bf polynomial interpolation-Lagrange"))
exportgraphics(gca,'../doc/fig/Q1.emf','BackgroundColor','none','ContentType','auto','Resolution',800);

%% Q2
figure('Name','Q2')
h = plot(x,[1./(1 + 25 * x.^2);Q2_L_10]);
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel("\fontname{Times New Roman} \fontsize{10} \it x")
ylabel("\fontname{Times New Roman} \fontsize{10} \it y")
legend( h, '$$f(x) = \frac{1}{1 + 25 x^2}$$','$$L_{10}^{(2)}(x) = \sum_{k = 0}^{10} f(x_k) l_k^{(2)}(x)$$','Location', 'south', 'Interpreter', 'latex','fontsize',10 );
legend('boxoff')
title(sprintf("\\fontname{Times New Roman} \\fontsize{10} \\bf polynomial interpolation-Chebyshev"))
exportgraphics(gca,'../doc/fig/Q2.emf','BackgroundColor','none','ContentType','auto','Resolution',800);

%% Q3
figure('Name','Q3')
h = plot(x,[1./(1 + 25 * x.^2);Q3_L_10]);
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel("\fontname{Times New Roman} \fontsize{10} \it x")
ylabel("\fontname{Times New Roman} \fontsize{10} \it y")
legend( h, '$$f(x) = \frac{1}{1 + 25 x^2}$$','$$L_{10}^{(3)}(x) = \sum_{k = 0}^{10} f(x_k) l_k^{(3)}(x)$$','Location', 'best', 'Interpreter', 'latex','fontsize',10 );
legend('boxoff')
title(sprintf("\\fontname{Times New Roman} \\fontsize{10} \\bf piecewise linear interpolation"))
exportgraphics(gca,'../doc/fig/Q3.emf','BackgroundColor','none','ContentType','auto','Resolution',800);

%% 局部函数

%% Q1

function output = Q1_interp_base(x,k,n)
    arguments
        x double
        k uint8 
        n double = 10
    end
    x_nodes = -1:2/n:1;
    output = ones(size(x));
    for i = 0:n
        if i ~= k
            output = output .* (x - x_nodes(i+1)) ./ (x_nodes(k+1) - x_nodes(i+1));
        end
    end
end

%% Q2

function output = Q2_interp_base(x,k,n)
    arguments
        x double
        k uint8 
        n double = 10
    end
    x_nodes = cos(pi * (2 * (0:n) + 1) ./ (2*n + 2));
    output = ones(size(x));
    for i = 0:n
        if i ~= k
            output = output .* (x - x_nodes(i+1)) ./ (x_nodes(k+1) - x_nodes(i+1));
        end
    end
end

%% Q3

function output = Q3_interp_base(x,k,n)
    arguments
        x double
        k uint8 
        n double = 10
    end
    x_nodes = -1:2/n:1;
    output = zeros(size(x));
    if k == 0
        idx = (x >= x_nodes(0+1)) & (x <= x_nodes(1+1));
        output(idx) = (x(idx) - x_nodes(1+1)) ./ (x_nodes(0+1) - x_nodes(1+1));
    elseif k == n
        idx = (x >= x_nodes(n-1+1)) & (x <= x_nodes(n+1));
        output(idx) = (x(idx) - x_nodes(n-1+1)) ./ (x_nodes(n+1) - x_nodes(n-1+1));
    else
        idx_left = (x <= x_nodes(k+1)) & (x >= x_nodes(k-1+1));
        idx_right = (x >= x_nodes(k+1)) & (x <= x_nodes(k+1+1));
        output(idx_left) = (x(idx_left) - x_nodes(k-1+1)) ./ (x_nodes(k+1) - x_nodes(k-1+1));
        output(idx_right) = (x(idx_right) - x_nodes(k+1+1)) ./ (x_nodes(k+1) - x_nodes(k+1+1));
    end
end
