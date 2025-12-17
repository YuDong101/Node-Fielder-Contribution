clear; clc;
% load C.ele_chemical_adjacency.mat;
rng(1)

% Directed cycle
% A = [0 1 0 0 0 0 0;
%      0 0 1 0 0 0 0;
%      0 0 0 1 0 0 0;
%      0 0 0 0 1 0 0;
%      0 0 0 0 0 1 0;
%      0 0 0 0 0 0 1
%      1 0 0 0 0 0 0];

% Random cycle
% A = [0 1 0 0 0 0 0;
%      0 0 1 0 0 0 0;
%      1 0 0 1 0 0 0;
%      0 0 0 0 1 1 0;
%      0 1 0 0 0 1 0;
%      0 1 0 0 0 0 1
%      1 0 1 0 0 0 0];

% Hub motif
A = [0 1 0 0 0 0 0;
     0 0 1 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 1 0 1;
     0 0 0 0 0 1 1;
     1 0 0 0 0 0 1
     1 1 1 0 0 0 0];

N = size(A,1); glin=0.1; glor=10;

out = fiedler_directed_weighted(A);
[~, ord] = sort(out.node_fiedler_in, 'descend');
max_Nodes = ord';
max_Nodes_Contrib = out.node_fiedler_in(ord);


D = diag(sum(A, 2));
L = D - A;

Tend  = 100;  nTime=10000; nRun  = 1000;
[tspan, E_lin, E_lor, X_lin_last, Xx_lor_last] = ...
    simulate_linear_and_lorenz(L, glin, glor, Tend, nRun, nTime);

tau_lin = estimate_decay_tau(tspan, E_lin(:,:));
tau_lor = sync_time_lor(E_lor, tspan);
%%
figure(101);
subplot(2,2,1); semilogy(tspan, E_lin); xlim([0 100]);
title('Linear diffusion  |x_i - mean(x)|');

subplot(2,2,2); plot(tspan, E_lor); xlim([0 10]);
title('Lorenz network  |x_i - mean(x)|');

subplot(2,2,3); bar(tau_lin); xlabel('node'); ylabel('\tau_{lin}');
subplot(2,2,4); bar(tau_lor); xlabel('node'); ylabel('\tau_{Lorenz}');

%%
Node_feilder = ...
out.node_fiedler_in;

[~,idxx]=sort(Node_feilder); rank_feilder = zeros(size(Node_feilder)); rank_feilder(idxx) = 1:numel(Node_feilder);
[~,idxx]=sort(tau_lin); rank_tao = zeros(size(tau_lin)); rank_tao(idxx) = 1:numel(tau_lin);

figure(2)
subplot(311),bar(((Node_feilder)));
subplot(312),bar(tau_lin);
subplot(313),bar(tau_lor);

figure(32)
subplot(211),plot(log10(abs(Node_feilder)),tau_lin,'ro');
subplot(212),plot(log10(abs(Node_feilder)),tau_lor,'ro');

[rlin,plin]=corr(log10(abs(Node_feilder)),tau_lin,"Type","Pearson");
[rlor,plor]=corr(log10(abs(Node_feilder)),tau_lor,"Type","Pearson");
fprintf('lin r %4f p %d mean %f.\n', rlin, plin, mean(tau_lin));
fprintf('lor r %4f p %d mean %f.\n', rlor, plor, mean(tau_lor));
%%

figure(1);
subplot(211),plot(tspan, E_lin, 'LineWidth', 0.9);
axis([0 30,-inf inf]);
grid on; box on;
xlabel('t'); ylabel('e_i(t) = |\xi_i - mean(\xi)|');
figure(1);
subplot(212),plot(tspan, X_lin_last, 'LineWidth', 0.9);
grid on; box on;
xlabel('t'); ylabel('e_i(t) = |\xi_i - mean(\xi)|');
figure(222);
subplot(211),plot(tspan, E_lor, 'LineWidth', 0.9);
axis([0 30,-inf inf]);
grid on; box on;
xlabel('t'); ylabel('e_i(t) = |\xi_i - mean(\xi)|');
subplot(212),plot(tspan, Xx_lor_last, 'LineWidth', 0.9);
axis([0 30,-inf inf]);
grid on; box on;
xlabel('t'); ylabel('e_i(t) = |\xi_i - mean(\xi)|');
%
figure(4001)
h = heatmap(X_lin_last);
colormap(jet);
colorbar;
h.GridVisible = 'off';
h.Title = 'Linear System State at Last Time Step';
h.XLabel = 'Node Index';
h.YLabel = 'Node Index';
figure(4002)
h = heatmap(Xx_lor_last);
colormap(jet);
colorbar; 
h.GridVisible = 'off';
h.Title = 'Linear System State at Last Time Step';
h.XLabel = 'Node Index';
h.YLabel = 'Node Index';