function out = fiedler_directed_weighted(A, varargin)
% Fiedler value & left/right vectors & per-edge contributions
% for directed, weighted, possibly NON-SYMMETRIC adjacency A.
%
% 定义（Jiang et al., PRE 109, 054301, 2024）:
%   L = D_out - A, D_out = diag(sum_j A(i,j))
%   取 L 的“非零特征值中实部最小”的那个 λ2，Fiedler 值 = Re(λ2)
%   左/右 Fiedler 向量 y2, x2 满足 y2'*x2 = 1
%   边贡献: I^{wd}_{i->j} = A(i,j) * Re( y2(j) * ( x2(j) - x2(i) ) )
%
% 用法:
%   out = fiedler_directed_weighted(A)
%   out = fiedler_directed_weighted(A,'convention','row','use_eigs',true)
%
% 参数(name/value):
%   'convention' : 'row'(默认) 行出度; 若你的 A 是列出度，请传 'col' 或先转置
%   'zero_tol'   : 判零阈值(特征值) 默认 1e-10
%   'use_eigs'   : 自动选择(默认 [])；true 强制 eigs；false 强制 eig
%   'k'          : eigs 候选特征数量(默认 20)
%   'sigma'      : eigs 的移位(默认 1e-8)
%   'tol'        : eigs/svds 容差(默认 1e-10)
%   'maxit'      : eigs 最大迭代(默认 2000)
%   'add_selfloop_eps' : 对出度为0的点加自环权 eps (默认 0) —— 稳定数值
%   'teleport'   : [false] 若 true，对 A 做行随机 + teleportation 归一化再回到权图
%   'alpha'      : teleport 参数(默认 0.99)
%   'televec'    : teleport 目标分布，'uniform'或 N×1 向量(默认 'uniform')
%   'verbose'    : true/false 打印自检(默认 true)

% ---- 解析参数 ----
p = inputParser;
addParameter(p,'convention','row');
addParameter(p,'zero_tol',1e-10);
addParameter(p,'use_eigs',[]);
addParameter(p,'k',20);
addParameter(p,'sigma',1e-8);
addParameter(p,'tol',1e-10);
addParameter(p,'maxit',2000);
addParameter(p,'add_selfloop_eps',0);
addParameter(p,'teleport',false);
addParameter(p,'alpha',0.99);
addParameter(p,'televec','uniform');
addParameter(p,'verbose',true);
parse(p,varargin{:});
opt = p.Results;

% ---- 基本清洗 & 方向约定 ----
A = double(A);
A(A<0) = 0;                        % 非负
N = size(A,1);
assert(size(A,2)==N,'A must be square.');

% 若用户的 A 采用“列出度”记法，则转置使其符合行出度约定
if strcmpi(opt.convention,'col')
    A = A.';
end

% 可选：死节点自环 & teleportation (让马氏链遍历/更稳)
if opt.add_selfloop_eps > 0
    drow = sum(A,2);
    dead = (drow==0);
    A(dead,dead) = A(dead,dead) + opt.add_selfloop_eps;
end
if opt.teleport
    % 行随机化 + teleportation，然后再映射回“权图”解释(用于稳定谱)
    s = sum(A,2); s(s==0) = 1;
    P = spdiags(1./s,0,N,N) * sparse(A);
    if ischar(opt.televec) && strcmpi(opt.televec,'uniform')
        v = ones(N,1)/N;
    else
        v = opt.televec(:); v = v / sum(v);
    end
    P = opt.alpha*P + (1-opt.alpha)*(v*ones(1,N));
    % 将 P 转回等“强度”的权：A_tp 使得行和与原 A 相同
    row_sum = sum(A,2);
    A = spdiags(row_sum,0,N,N) * P;          % 仍满足行出度的度和
    A = full(A);
end

% ---- 构造 L = D_out - A ----
dout = sum(A,2);
L = spdiags(dout,0,N,N) - sparse(A);

% ---- 自动选择 eig / eigs ----
if isempty(opt.use_eigs)
    use_eigs = (issparse(L) && N > 600);
else
    use_eigs = logical(opt.use_eigs);
end

% ---- 求右模 (λ, x) —— 非零特征里 Re(λ) 最小者 ----
if ~use_eigs
    % 小规模: 直接 eig(full)
    [V, Dm] = eig(full(L));
    lam = diag(Dm);
else
    % 大规模稀疏: eigs + shift-invert 在 0 附近抓一批
    k = min([opt.k, N-1, 50]); if k < 2, k = min(N-1, 2); end
    sigma = opt.sigma;
    opts = struct('tol',opt.tol,'maxit',opt.maxit);
    try
        [V, Dm] = eigs(L, k, sigma, opts);
        lam = diag(Dm);
    catch
        warning('eigs failed; fallback to eig(full(L)).');
        [V, Dm] = eig(full(L)); lam = diag(Dm);
        use_eigs = false;
    end
end

% 选择 λ2: 在“非零”特征里找 Re(λ) 最小
nz = find(abs(lam) > opt.zero_tol);
assert(~isempty(nz), 'All eigenvalues are numerically zero?');
[~, ix] = min(real(lam(nz)));
idx2 = nz(ix);
lambda2 = lam(idx2);
x2 = V(:, idx2);

% ---- 求左模 y: 解 (L.' - λ2 I) y = 0 的最小奇异向量 ----
M = (L.' - lambda2*speye(N));
y2 = [];

% 优先 svds (最小奇异向量)
try
    [~,~,Vsvd] = svds(M, 1, 'smallest', struct('tol',opt.tol,'maxit',opt.maxit));
    y2 = Vsvd(:, end);
catch
    % 退化: 小规模用 svd，或用 null
    try
        [~,~,Vsvd] = svd(full(M),'econ');
        y2 = Vsvd(:, end);
    catch
        y2 = null(full(M),'r');
        if isempty(y2)
            % 最后手段：在 eig(L.') 中选与 x2 相关性最大的列
            [W,~] = eig(full(L.'));
            [~,jx] = max(abs(W' * x2));
            y2 = W(:,jx);
        else
            y2 = y2(:,1);
        end
    end
end


% ---- 归一化: y2'*x2 = 1 ----
alpha = y2.' * x2;
if abs(alpha) < 1e-12
    x2 = x2 / max(norm(x2), eps);
    alpha = y2.' * x2;
end
y2 = y2 / alpha;   % 现在 y2'*x2 == 1

% ---- 边贡献 I^{wd}_{i->j} ----
[i_idx, j_idx, w] = find(A);
edges   = [j_idx, i_idx];

contrib = w .* ( y2(i_idx) .* ( x2(i_idx) - x2(j_idx) ) );  % 允许是复数
sum_c = sum(contrib);                                       % 复数和
                       % 应≈ real(lambda2)
% ---- 节点贡献 I^{wd}_{i->j} ----
contrib_r = real(contrib);                 % 与 out.contrib 一致
node_fiedler_out = accumarray(j_idx, contrib_r, [N,1], @sum, 0);  % 出边累加
node_fiedler_in  = accumarray(i_idx, contrib_r, [N,1], @sum, 0);  % 入边累加
node_fiedler_all = node_fiedler_out + node_fiedler_in;            % 出入一起

% 可选：一致性检查（应当满足）
check_out = sum(node_fiedler_out);   % ≈ real(lambda2)
check_in  = sum(node_fiedler_in);    % ≈ real(lambda2)
check_all = sum(node_fiedler_all);   % ≈ 2*real(lambda2)

% ---- 打包输出 ----
out = struct();
out.lambda2      = lambda2;
out.alllambda    = lam;
out.lambda2_real = real(lambda2);
out.x2           = x2;
out.y2           = y2;
out.edges        = edges;          % [src, dst]
out.contrib      = real(contrib);        % per-edge contributions
out.check_sum = real(sum_c);      % ≈ Re(lambda2)
out.L            = L;

% 新增：节点层面的 Fiedler 值
out.node_fiedler_out = node_fiedler_out;  % 出边累加
out.node_fiedler_in  = node_fiedler_in;   % 入边累加
out.node_fiedler_all = node_fiedler_all;  % 出+入

% 可选：检查量
out.check_node_sum_out = check_out;       % ≈ out.check_sum
out.check_node_sum_in  = check_in;        % ≈ out.check_sum
out.check_node_sum_all = check_all;       % ≈ 2*out.check_sum
end