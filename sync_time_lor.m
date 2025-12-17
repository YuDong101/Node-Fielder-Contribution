function tau_node = sync_time_lor(E, t, varargin)
% COMPUTE_SYNC_TIME_PER_NODE
% 使用“首达 + 驻留”准则，计算每个节点达到近似同步的特征时间。
%
% 输入：
%   E : N x T 的误差矩阵，E(i,t) = |x_i(t) - mean(x(t))| 或其它误差定义
%   t : 1 x T 的时间向量（等步长）
%
% 可选参数（Name-Value）：
  % 'RelThresh'  : 相对阈值系数 eps_rel，thr_i = eps_rel * E(i,1)
%                  默认 0.05（即衰减到初始误差的 5%）
%   'AbsThresh'  : 绝对阈值 thr，可以是标量或长度为 N 的向量。
%                  如果提供，则忽略 RelThresh。
%   'DwellTime'  : 要求在阈值以下连续停留的时间（同 t 的单位），默认 1
%
% 输出：
%   tau_node : N x 1，每个节点的同步特征时间。
%              若该节点在仿真区间内从未满足准则，则为 NaN。
%
% 用法示例：
%   tau_lor = compute_sync_time_per_node(E_lor, tspan, ...
%               'RelThresh', 0.05, 'DwellTime', 1.0);

    [N, T] = size(E);

    if numel(t) ~= T
        error('E 的列数 (%d) 必须等于 t 的长度 (%d)。', T, numel(t));
    end

    % ---------- 解析可选参数 ----------
    p = inputParser;
    addParameter(p, 'RelThresh', 0.1, @(x)isnumeric(x) && isscalar(x) && x>0); % 0.1
    addParameter(p, 'AbsThresh', [], @(x)isnumeric(x));
    addParameter(p, 'DwellTime', 10.0, @(x)isnumeric(x) && isscalar(x) && x>=0); % 5.0
    parse(p, varargin{:});

    eps_rel   = p.Results.RelThresh;
    abs_thr   = p.Results.AbsThresh;
    dwellTime = p.Results.DwellTime;

    % ---------- 基本量 ----------
    dt = t(2) - t(1);
    if any(abs(diff(t) - dt) > 1e-9*abs(dt))
        warning('t 似乎不是严格等步长，DwellTime 的步数近似处理，请注意。');
    end
    dwellStep = max(1, round(dwellTime / dt));

    tau_node = NaN(N,1);

    % ---------- 主循环：逐节点计算首达时间 ----------
    for i = 1:N
        % 阈值确定
        if ~isempty(abs_thr)
            if isscalar(abs_thr)
                thr_i = abs_thr;
            else
                thr_i = abs_thr(i);
            end
        else
            thr_i = eps_rel; %  * E(i,1)
        end

        below = E(i,:) < thr_i;   % 是否低于阈值
        idx   = find(below, 1, 'first');

        while ~isempty(idx)
            % 检查从 idx 开始是否连续 dwellStep 个点都在阈值下面
            iEnd = min(idx + dwellStep - 1, T);
            if all(below(idx:iEnd))
                tau_node(i) = t(idx);
                break;
            else
                % 找下一段可能满足条件的位置
                idx2 = find(below(idx+1:end), 1, 'first');
                if isempty(idx2)
                    break;
                end
                idx = idx + idx2;
            end
        end
    end
end
