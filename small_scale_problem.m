%% 使用BB方法和最速下降法并比较结果
% 参数说明：
% X：每次迭代时自变量 X 的值
% K：迭代次数
% F: 每次迭代的目标函数值
% dF_norm: 每次迭代的梯度范数

           
clc;clear;close all;
A = [2 3;
     1 -1;
     -2 4];
b = [5; 
     -10; 
     -15];  %x_answer = [3; -1]   (norm(A*[3;-1]-b))^2 = 225


% 函数参数含义：BB_method(start, eps, alpha, alpham, alphaM, M, c1, beta,A,b)
[K, x_min, f_min, X, F, dF_norm] = BB_method([-10; -1], 1e-6, 1, 0.001, 20000, 3, 0.2, 0.3,A,b);
disp(['利用BB方法: 迭代次数 = ', num2str(K), ', 最优值点 x1 = ', num2str(x_min(1)),', x2 = ', num2str(x_min(2)), ', 最小值 = ', num2str(f_min)]);

[K1, x_min1, f_min1, X1, F1, dF_norm1] = DG_method([-10; -1], 1e-6,A,b);
disp(['利用最速下降法：迭代次数 = ', num2str(K1), ', 最优值点 x2 = ', num2str(x_min1(1)),', x2 = ', num2str(x_min(2)), ', 最小值 = ', num2str(f_min1)]);

% 绘制图形
Go_plot( X, F, dF_norm, X1, F1, dF_norm1,A,b);


% Barzilai-Borwein 方法
function [K, x_min, f_min, X, F, dF_norm] = BB_method(start, eps, alpha, alpham, alphaM, M, c1, beta,A,b)
    % 如果没有提供参数，则使用默认值
    if nargin < 9, beta = 0.3; end
    if nargin < 8, c1 = 0.2; end
    if nargin < 7, M = 10; end
    if nargin < 6, alphaM = 2; end
    if nargin < 5, alpham = 0.05; end
    if nargin < 4, alpha = 1; end
    
    % 初始化记录数据
    X = zeros(2, 1000);
    F = zeros(1, 1000);
    dF = zeros(2, 1000);
    dF_norm = zeros(1, 1000);
    
    % 初始值设置
    X(:, 1) = start;
    F(1) = f(X(:, 1),A,b);
    dF(:, 1) = df(X(:, 1),A,b);
    dF_norm(1) = norm(dF(:, 1));
    
    k = 1;
    while dF_norm(k) >= eps  %%%对应算法6.2第2行
        % 对步长进行修正
        while f(X(:, k) - alpha * dF(:, k),A,b) >= max(F(max(1, k-M):k)) - c1 * alpha * dF_norm(k)^2  %%%对应算法6.2第3行
            alpha = beta * alpha;  %%%对应算法6.2第4行
        end  %对应算法6.2第5行
        
        % 更新迭代点
        X(:, k + 1) = X(:, k) - alpha * dF(:, k);  %%%对应算法6.2第6行
        F(k + 1) = f(X(:, k + 1),A,b);
        dF(:, k + 1) = df(X(:, k + 1),A,b);
        dF_norm(k + 1) = norm(dF(:, k + 1));
        
        % 计算新的步长
        s = X(:, k + 1) - X(:, k);
        y = dF(:, k + 1) - dF(:, k);
        alpha = min([max([dot(s, s) / dot(s, y), alpham]), alphaM]);  %%%对应算法6.2第7行
        % alpha = min([max([dot(s, y) / dot(y, y), alpham]), alphaM]);

        k = k + 1;  %对应算法6.2第8行
    
    end  %%%对应算法6.2第9行
    
    K = k;
    x_min = X(:, k);
    f_min = F(k);
    X = X(:,1:K);
end

% 最速下降法
function [K, x_min, f_min, X, F, dF_norm] = DG_method(start, eps,A,b)
    % 初始化记录数据
    X = zeros(2, 1000);
    F = zeros(1, 1000);
    dF = zeros(2, 1000);
    dF_norm = zeros(1, 1000);
    
    % 初始值设置
    X(:, 1) = start;
    F(1) = f(X(:, 1),A,b);
    dF(:, 1) = df(X(:, 1),A,b);
    dF_norm(1) = norm(dF(:, 1));
    
    k = 1;
    while dF_norm(k) >= eps
        % 计算步长
        alpha = dot(dF(:, k), dF(:, k)) / dot( A'*A * dF(:, k), dF(:, k)); %精确线搜索
        
        % 更新迭代点
        X(:, k + 1) = X(:, k) - alpha * dF(:, k);
        F(k + 1) = f(X(:, k + 1),A,b);
        dF(:, k + 1) = df(X(:, k + 1),A,b);
        dF_norm(k + 1) = norm(dF(:, k + 1));
        
        k = k + 1;
    end
    
    K = k;
    x_min = X(:, k);
    f_min = F(k);
    X = X(:,1:K);
end

% 目标函数
function val = f(x,A,b)
    val = (1/2)*(norm(A*x-b,2))^2;
end

% 梯度函数
function grad = df(x,A,b)
    grad = A' * A * x - A' * b;
end

% 绘制图形函数
function Go_plot(X1, F1, dF_norm1,X2, F2, dF_norm2, A, b)
    % 获取精确解
    x_optimal = A \ b;
    f_optimal = 0.5*norm(A*x_optimal - b)^2;
    
    % 创建动态等高线范围
    x1_min = min([X1(1,:), X2(1,:)]) - 1;
    x1_max = max([X1(1,:), X2(1,:)]) + 1;
    x2_min = min([X1(2,:), X2(2,:)]) - 1;
    x2_max = max([X1(2,:), X2(2,:)]) + 1;
    
    % 生成网格
    [XX1, XX2] = meshgrid(linspace(x1_min, x1_max, 100),...
                   linspace(x2_min, x2_max, 100));
    ZZ = arrayfun(@(x,y) 0.5*norm(A*[x;y]-b)^2, XX1, XX2);

    % 配色方案
    color_bb = [0.9290 0.6940 0.1250];  % 金色
    color_sd = [0 0.4470 0.7410];       % 蓝色
    opt_color = [0.6350 0.0780 0.1840]; % 深红色
    
    %% 图1: 优化路径对比
    figure('Position', [100 100 1200 500])
    
    % 子图1: 带等高线的优化路径
    subplot(1,2,1)
    contour(XX1, XX2, ZZ, 15, 'LineWidth', 0.8, 'LineColor', [0.6 0.6 0.6]);
    hold on
    p1 = plot(X1(1,:), X1(2,:), 'o-', 'Color', color_bb,...
        'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', color_bb);
    p2 = plot(X2(1,:), X2(2,:), 's-', 'Color', color_sd,...
        'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', color_sd);
    p3 = plot(x_optimal(1), x_optimal(2), 'pentagram', 'Color', opt_color,...
        'MarkerSize', 15, 'LineWidth', 2, 'MarkerFaceColor', opt_color);
    axis tight equal
    grid on
    set(gca, 'FontSize', 11, 'LineWidth', 1.2)
    xlabel('x_1', 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('x_2', 'FontSize', 12, 'FontWeight', 'bold')
    title('优化路径对比', 'FontSize', 14, 'FontWeight', 'bold')
    legend([p1 p2 p3], {'BB方法','最速下降法','最优解'}, 'Location', 'best')
    colormap(jet)
    
    % 子图2: 梯度范数下降
    subplot(1,2,2)
    semilogy(dF_norm1, 'o-', 'Color', color_bb,...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', color_bb)
    hold on
    semilogy(dF_norm2, 's-', 'Color', color_sd,...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', color_sd)
    grid on
    set(gca, 'FontSize', 11, 'LineWidth', 1.2)
    xlabel('迭代次数', 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('梯度范数', 'FontSize', 12, 'FontWeight', 'bold')
    title('梯度下降过程', 'FontSize', 14, 'FontWeight', 'bold')
    legend({'BB方法','最速下降法'}, 'Location', 'northeast')
    ylim([1e-7, max([dF_norm1, dF_norm2])*2])

    %% 图2: 目标函数收敛过程
    figure('Position', [100 100 600 400])
    F1_diff = F1 - f_optimal;
    F2_diff = F2 - f_optimal;
    
    semilogy(F1_diff, 'o-', 'Color', color_bb,...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', color_bb)
    hold on
    semilogy(F2_diff, 's-', 'Color', color_sd,...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', color_sd)
    grid on
    set(gca, 'FontSize', 11, 'LineWidth', 1.2)
    xlabel('迭代次数', 'FontSize', 12, 'FontWeight', 'bold')
    ylabel('f(x) - f^*', 'FontSize', 12, 'FontWeight', 'bold')
    title('目标函数收敛过程', 'FontSize', 14, 'FontWeight', 'bold')
    legend({'BB方法','最速下降法'}, 'Location', 'best')
    ylim([1e-7, max([F1_diff, F2_diff])*2])
end
