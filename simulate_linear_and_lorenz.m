function [tspan, E_lin_avg, E_lor_avg, X_lin_last, Xx_lor_last] = ...
    simulate_linear_and_lorenz(L, glin, glor, Tend, nRun, nTime, sigma, rho, beta)
%SIMULATE_LINEAR_AND_LORENZ
% 在同一个时间循环里同时仿真：
%   1) 线性扩散系统  xdot_lin = -g * L * x_lin
%   2) 洛伦兹混沌网络，每个节点一个 Lorenz 振子，在 x 分量上用 L 耦合
%
% 输入：
%   L     : N×N 拉普拉斯矩阵（例如 L = D - A）
%   g     : 耦合强度（线性和洛伦兹都用同一个 g，如需区分可自己扩展）
%   Tend  : 仿真结束时间，例如 50 或 100
%   nRun  : 初值重复次数，用于对 E 做平均
%   nTime : 采样时间点个数，例如 5000 或 10000
%   sigma, rho, beta : 洛伦兹参数（可选，缺省为 10,28,8/3）
%
% 输出：
%   tspan      : 1×T 时间向量
%   E_lin_avg  : N×T，线性系统的平均误差 |x_i - mean(x)|
%   E_lor_avg  : N×T，洛伦兹系统 x 分量的平均误差 |x_i - mean(x)|
%   X_lin_last : N×T，最后一次 run 的线性系统轨迹
%   Xx_lor_last: N×T，最后一次 run 的洛伦兹 x_i(t) 轨迹

    if nargin < 9, beta  = 8/3; end
    if nargin < 8, rho   = 28;  end
    if nargin < 7, sigma = 10;  end

    N = size(L, 1);

    % 时间网格
    tspan = linspace(0, Tend, nTime);

    % 误差累积
    E_lin_acc = zeros(N, numel(tspan));
    E_lor_acc = zeros(N, numel(tspan));

    for run = 1:nRun
        % ===== 初值设置 =====
        % 线性系统初值
        x_lin0 = -1 + 2*rand(N,1);

        % 洛伦兹系统初值 (x,y,z)
        x0 = -20 + 40*rand(N,1);
        y0 = -30 + 60*rand(N,1);
        z0 = 0 + 50*rand(N,1);

        % 统一拼成 4N 维状态向量：
        % [ x_lin;
        %   x_lor;
        %   y_lor;
        %   z_lor ]
        X0 = zeros(4*N,1);
        X0(1:N)           = x_lin0;
        X0(N+1:2*N)       = x0;
        X0(2*N+1:3*N)     = y0;
        X0(3*N+1:4*N)     = z0;

        % ===== 一次 ODE45 同时积分两个系统 =====
        odefun = @(t,X) mixed_ode(t, X, L, glin, glor, sigma, rho, beta);

        sol = ode45(odefun, [tspan(1) tspan(end)], X0, ...
            odeset('RelTol',1e-6,'AbsTol',1e-8));

        Xall = deval(sol, tspan);   % 4N × T

        % 拆出线性系统和洛伦兹 x 分量
        X_lin  = Xall(1:N,          :);   % 线性系统 x_i(t)
        Xx_lor = Xall(N+1:2*N,      :);   % 洛伦兹 x_i(t)
        
        % ===== 计算两个系统各自的 |x_i - mean(x)| =====
        s_lin = mean(X_lin, 1);
        E_lin(run,:,:) = abs(X_lin - s_lin);

        s_lor = mean(Xx_lor, 1);
        E_lor(run,:,:) = abs(Xx_lor - s_lor) / 20;

        % 累积
        E_lin_acc = E_lin_acc + squeeze(E_lin(run,:,:));
        E_lor_acc = E_lor_acc + squeeze(E_lor(run,:,:));
    end

    
    % 平均
    E_lin_avg  = E_lin_acc / nRun;
    E_lor_avg  = E_lor_acc / nRun;

    % 最后一次 run 的轨迹也返回，方便后面画图
    X_lin_last  = X_lin;
    Xx_lor_last = Xx_lor;
end


% ================== 本地子函数：联合 ODE ==================
function dXdt = mixed_ode(~, X, L, glin, glor, sigma, rho, beta)
% 一个 4N 维系统，前 N 维是线性扩散，后 3N 维是洛伦兹网络

    N = size(L, 1);

    % 拆分变量
    x_lin = X(1:N);
    x     = X(N+1:2*N);
    y     = X(2*N+1:3*N);
    z     = X(3*N+1:4*N);

    % ---- 1) 线性扩散 ----
    dx_lin = - glin * (L * x_lin);

    % ---- 2) 洛伦兹网络（只在 x 分量上耦合）----
    dx = sigma .* (y - x) - glor * (L * x);
    dy = x .* (rho - z) - y;
    dz = x .* y - beta .* z;

    % 拼回 4N 维
    dXdt              = zeros(4*N,1);
    dXdt(1:N)         = dx_lin;
    dXdt(N+1:2*N)     = dx;
    dXdt(2*N+1:3*N)   = dy;
    dXdt(3*N+1:4*N)   = dz;
end
