function tau_1e = estimate_decay_tau(t, y, varargin)

p = inputParser;
addParameter(p,'method','lsq');
addParameter(p,'robust',true);
addParameter(p,'winFrac',0.2);
addParameter(p,'minAmp',0.05);
addParameter(p,'window',[]);
addParameter(p,'profileC',false);
parse(p,varargin{:});
opt = p.Results;


t = t(:);                     % Nt x 1
if isvector(y)
    y = y(:).';               % 1 x Nt
end
[M, Nt] = size(y);
assert(Nt==numel(t), 'y 的列数必须等于 t 的长度（Nt）。');


out = struct();
out.method    = lower(opt.method);
out.t         = t(:).';                    % 1xNt
out.tau       = nan(M,1);
out.A         = nan(M,1);
out.C         = nan(M,1);
out.ci_tau    = nan(M,2);
out.tau_1e    = nan(M,1);
out.tau_ar1   = nan(M,1);
out.mask_used = false(M,Nt);
out.window    = nan(M,2);
out.yfit_all  = nan(M,Nt);
out.resid_all = nan(M,Nt);
out.notes     = cell(M,1);

if opt.profileC
    out.profileC = struct();
    out.profileC.C_grid   = cell(M,1);
    out.profileC.tau_grid = cell(M,1);
    out.profileC.SSE      = cell(M,1);
    out.profileC.C_star   = nan(M,1);
    out.profileC.tau_star = nan(M,1);
end


for m = 1:M
    yy = y(m,:).';
    s = single_series(t, yy, opt);

    out.tau(m)     = s.tau;
    out.A(m)       = s.A;
    out.C(m)       = s.C;
    out.tau_1e(m)  = s.tau_1e;
    out.tau_ar1(m) = s.tau_ar1;
    out.mask_used(m,:) = s.mask_used(:).';
    out.window(m,:)    = s.window(:).';

    if ~isempty(s.yfit_all)
        out.yfit_all(m,:)  = s.yfit_all(:).';
        tmp_res = yy - s.yfit_all;
        out.resid_all(m,:) = tmp_res(:).';
    end
    if isfield(s,'ci_tau') && ~isempty(s.ci_tau)
        out.ci_tau(m,:) = s.ci_tau(:).';
    end
    out.notes{m} = s.notes;

    if opt.profileC && isfield(s,'profileC')
        out.profileC.C_grid{m}   = s.profileC.C_grid;
        out.profileC.tau_grid{m} = s.profileC.tau_grid;
        out.profileC.SSE{m}      = s.profileC.SSE;
        out.profileC.C_star(m)   = s.profileC.C_star;
        out.profileC.tau_star(m) = s.profileC.tau_star;
    end
end

tau_1e=out.tau_1e;


function out1 = single_series(t, y, opt)
    t = t(:); y = y(:);
    assert(numel(t)==numel(y) && numel(t)>3, 't,y 长度需一致且>3');


    if isempty(opt.window)
        C_tail = median(y(end-max(5,round(opt.winFrac*numel(y)))+1:end));
        tmin = t(1); tmax = t(end);
    else
        tmin = opt.window(1); tmax = opt.window(2);
        maskW = (t>=tmin & t<=tmax);
        C_tail = median(y(maskW & (t>tmin+(tmax-tmin)*0.8)));
    end

    mask_all = (t>=tmin & t<=tmax);
    A0 = y(find(mask_all,1,'first')) - C_tail;
    mask_use = mask_all & (y > C_tail + opt.minAmp*abs(A0));
    if nnz(mask_use) < 3
        mask_use = mask_all;
    end
    t1 = t(mask_use); y1 = y(mask_use);

    out1 = struct(); out1.method = lower(opt.method); out1.mask_used = mask_use;
    out1.notes = '';

    out1.tau_1e = local_tau_one_over_e(t1,y1,C_tail);
    out1.tau_ar1 = local_tau_ar1(t1,y1);

    switch lower(opt.method)
    case 'expfit'
        [tau,A,C,stats] = local_expfit(t1,y1,C_tail);
        out1.tau = tau; out1.A = A; out1.C = C;
        out1.ci_tau = stats.tau_ci95;
        % 全时间拟合
        mdl = @(p,tt) p(3) + p(2)*exp(-tt/p(1));
        out1.yfit_all = mdl([tau,A,C], t);
        out1.window = [tmin tmax];
        out1.notes = sprintf('expfit R2=%.3f', stats.R2);

    case 'lsq'
        [tau,A,C,ci_tau,notes,~,~] = local_lsq(t1,y1,C_tail,opt.robust);
        out1.tau = tau; out1.A = A; out1.C = C; out1.ci_tau = ci_tau;
        out1.window = [tmin tmax];
        out1.notes = notes;

        mdl = @(p,tt) p(3) + p(2)*exp(-tt/p(1));
        out1.yfit_all = mdl([tau,A,C], t);

        if opt.profileC
            span = 0.02*(max(y1)-min(y1));
            prof = local_profileC(t1,y1, C, span);
            out1.profileC = prof;
        end

    case 'one_over_e'
        out1.tau = out1.tau_1e; out1.A = y1(1)-C_tail; out1.C = C_tail;
        out1.ci_tau = []; out1.yfit_all = nan(size(t)); 
        out1.window = [tmin tmax];
        out1.notes = '1/e 粗估';

    case 'ar1'
        out1.tau = out1.tau_ar1; out1.A = NaN; out1.C = mean(y1);
        out1.ci_tau = []; out1.yfit_all = nan(size(t));
        out1.window = [tmin tmax];
        out1.notes = 'AR(1) 相关时间（仅供参考）';

    otherwise
        error('未知 method：%s', opt.method);
    end
end


function tau = local_tau_one_over_e(t,y,C)

    C = local_baseline(t,y,C);
    y0 = y(1);
    target = C + (y0 - C)/exp(1);

    ymin = min(y); ymax = max(y);
    if target < ymin || target > ymax
        tau = NaN; return;
    end

    z = y - C;
    z = max(z, 0);
    z_env = cummin(z);
    y_env = C + z_env;


    k = find(y_env <= target, 1, 'first');

    if isempty(k)
        tau = NaN; return;
    elseif k == 1
        tau = t(1); return;
    else
        t1 = t(k-1); t2 = t(k);
        y1 = y(k-1); y2 = y(k);
        if y2 == y1
            tau = t1;
        else
            tau = t1 + (target - y1) * (t2 - t1) / (y2 - y1);
        end
    end
end


function tau = local_tau_ar1(t,y)
    dt = median(diff(t));
    x  = y - mean(y);
    if numel(x) < 3, tau = NaN; return; end
    phi = corr(x(1:end-1), x(2:end), 'rows','complete');
    if ~isfinite(phi) || phi<=0 || phi>=1
        tau = NaN;
    else
        tau = -dt / log(phi);
    end
end

function C = local_baseline(~,y,C0)
    if nargin<3 || isempty(C0)
        m = max(5, round(0.2*numel(y)));
        C = median(y(end-m+1:end));
    else
        C = C0;
    end
end

function [tau,A,C,stats] = local_expfit(t,y,C0)
    C = local_baseline(t,y,C0);
    z = y - C; mask = z>0;
    t1 = t(mask); z1 = z(mask);
    X = [ones(numel(t1),1) t1];
    b = X \ log(zk_safe(z1));
    A   = exp(b(1));
    tau = -1/b(2);
    res = log(zk_safe(z1)) - X*b;
    s2  = sum(res.^2)/(numel(t1)-2);
    covb= s2*inv(X.'*X);
    se_slope = sqrt(covb(2,2));
    stats.tau_ci95 = tau + [-1,1] * (se_slope/abs(b(2))^2) * 1.96;
    stats.R2 = 1 - sum(res.^2)/sum((log(zk_safe(z1))-mean(log(zk_safe(z1)))).^2);
end

function [tau,A,C,ci_tau,notes,yfit,resid] = local_lsq(t,y,C0,robustFlag)

    C = local_baseline(t,y,C0);
    t = t(:); y = y(:);
    dt = median(diff(t)); if ~isfinite(dt) || dt<=0, dt = max(t)-min(t); end
    if dt<=0, dt = 1; end
    t0 = t(1);
    ts = (t - t0) / dt;
    
    tau1e = local_tau_one_over_e(t,y,C);
    if ~isfinite(tau1e) || tau1e<=0, tau1e = max(t)-min(t); end
    tau0 = max(tau1e, eps);
    A0   = y(1) - C;
    p0s  = [max(tau0/dt, eps), A0, C];
    
    mdl_s = @(p,tt) p(3) + p(2)*exp(-tt/p(1));
    notes = '';
    
    try
        if exist('nlinfit','file')
            opts = statset('nlinfit');
            if robustFlag, opts.RobustWgtFun = 'bisquare'; end
            opts.MaxIter = 5000;
            opts.TolX   = 1e-10;
            opts.TolFun = 1e-12;
            [p_hat, R, J] = nlinfit(ts, y, mdl_s, p0s, opts);
            ci = nlparci(p_hat, R, 'jacobian', J);
            tau  = p_hat(1)*dt; A = p_hat(2); C = p_hat(3);
            ci_tau = ci(1,:)*dt;
            yfit = mdl_s(p_hat, (t - t0)/dt); resid = y - yfit;
            notes = 'nlinfit(稳健, 归一化时间)';
        else
            error('no_nlinfit');
        end
    catch
        try
            if exist('lsqcurvefit','file')
                opts = optimoptions('lsqcurvefit','Display','off','MaxIterations',5000,...
                                    'FunctionTolerance',1e-12,'StepTolerance',1e-10);
                lb = [0, -Inf, -Inf]; ub = [Inf, Inf, Inf];
                p_hat = lsqcurvefit(mdl_s, p0s, ts, y, lb, ub, opts);
                tau  = p_hat(1)*dt; A = p_hat(2); C = p_hat(3);
                ci_tau = []; yfit = mdl_s(p_hat, (t - t0)/dt); resid = y - yfit;
                notes = 'lsqcurvefit(归一化时间)';
            else
                error('no_lsqcurvefit');
            end
        catch
            fSSE = @(q) sum((mdl_s([exp(q(1)), q(2), q(3)], ts) - y).^2);
            q0 = [log(max(p0s(1),eps)), p0s(2), p0s(3)];
            qhat = fminsearch(fSSE, q0, optimset('Display','off'));
            p_hat = [exp(qhat(1)), qhat(2), qhat(3)];
            tau  = p_hat(1)*dt; A = p_hat(2); C = p_hat(3);
            ci_tau = []; yfit = mdl_s(p_hat, (t - t0)/dt); resid = y - yfit;
            notes = 'fminsearch fallback(归一化时间)';
        end
    end
end


function prof = local_profileC(t,y, C_center, span)
    C_grid = linspace(C_center - span, C_center + span, 41);
    SSE = zeros(size(C_grid));
    Tau = zeros(size(C_grid));
    for k = 1:numel(C_grid)
        Ck = C_grid(k);
        z  = y - Ck; m = z>0;
        tk = t(m); zk = z(m);
        if numel(tk) < 3, SSE(k)=NaN; Tau(k)=NaN; continue; end
        X = [ones(numel(tk),1) tk];
        b = X \ log(zk_safe(zk));
        Tau(k) = max(eps, -1/b(2));
        yhat = Ck + exp(b(1))*exp(-tk/Tau(k));
        SSE(k) = sum((yhat - y(m)).^2);
    end
    [~,kmin] = min(SSE);
    prof = struct('C_grid',C_grid,'tau_grid',Tau,'SSE',SSE, ...
                  'C_star',C_grid(kmin),'tau_star',Tau(kmin));
end

function [yfit,resid] = local_yfit_resid(t,y,p)
    mdl = @(pp,tt) pp(3) + pp(2)*exp(-tt/pp(1));
    yfit = mdl(p,t); resid = y - yfit;
end

function z1 = zk_safe(z)
    z1 = z;
    epsz = 1e-12;
    z1(z1<=epsz) = epsz;
end
end
