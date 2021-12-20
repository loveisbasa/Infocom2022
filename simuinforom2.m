clear all 
clc
T = 1e5;
sigma = 1.5;
mu = 1;
exp_num = 100;
AoI_zw_dataset = zeros(exp_num, T);
AoI_opt_dataset = zeros(exp_num, T);
AoI_online_dataset = zeros(exp_num, T);
AoI_online_l_dataset = zeros(exp_num, T);
C = 200;

for exp_idx = 1:exp_num

%% generate delay distribution
pd = makedist('Lognormal','mu',log(mu),'sigma',sigma)
delay = random(pd,T,1);
exp_idx

%% zero wait
R = zeros(1, T);
FrameLength = zeros(1, T);
Tau = zeros(1, T);
k = 0;
t = 0;
prev_delay = 0;
AoI_zw = zeros(1, T);
Frame_zw = zeros(1, T);
tot_aoi_prev = 0;
tot_frame_len = 0;
while k < T
    k = k + 1;
    current_delay = delay(k);
    R(k) = (2 * prev_delay + current_delay) * current_delay / 2 + C;
    FrameLength(k) = current_delay;
    prev_delay = current_delay;
    tot_aoi_prev = tot_aoi_prev + R(k);
    tot_frame_len = tot_frame_len + current_delay;
    AoI_zw(k) = tot_aoi_prev;
    Frame_zw(k) = tot_frame_len;
end
avg_AoI_zw = sum(R) / FrameLength(T);

%% Opt by Yin

% step 1. search for tau
l = 0;
r = 1e8;
mu_bar = log(mu);
logn_dens = @(x) 1 / (sqrt(2 * pi) .* x * sigma) * exp(-1/(2*sigma*sigma) * (log(x)-mu_bar).^2);
logn_x = @(x) x ./ (sqrt(2 * pi) .* x * sigma) .* exp(-1/(2*sigma*sigma) .* (log(x)-mu_bar).^2);
logn_x_2 = @(x) x.^2 ./ (sqrt(2 * pi) .* x * sigma) .* exp(-1/(2*sigma*sigma) .* (log(x)-mu_bar).^2);
while (r - l) > 0.0001
    b = (l + r) / 2;
    term1 = exp(mu_bar + sigma^2 / 2) + b * logncdf(b, mu_bar, sigma) - integral(logn_x, 0, b);
    term2 = (exp(mu_bar + sigma^2 / 2))^2 + (exp(sigma^2) - 1) * exp( 2 * mu_bar + sigma^2 ) + b^2 * logncdf(b, mu_bar, sigma) - integral(logn_x_2, 0, b) + 2 * C;
    if term1 - term2 / 2 / b > 0
        r = b;
    else
        l = b;
    end
end

%step 2, sumu
R = zeros(1, T);
Tau = zeros(1, T);
k = 0;
t = 0;
prev_delay = 0;
AoI_opt = zeros(1, T);
Frame_opt = zeros(1, T);
prev_dpw = 0;
tot_aoi_prev = 0;
tot_frame_len = 0;
while k < T
    k = k + 1;
    current_delay = delay(k);
    w = max(b - current_delay, 0);
    R(k) = (2 * prev_dpw + current_delay) * current_delay / 2 + (2 * current_delay + w) * w / 2 + C;
    prev_dpw = current_delay + w;
    tot_aoi_prev = tot_aoi_prev + R(k);
    tot_frame_len = tot_frame_len + prev_dpw;
    AoI_opt(k) = tot_aoi_prev;
    Frame_opt(k) = tot_frame_len;
end
avg_AoI_opt = sum(R) / t;
b_precise = b;

%% Online Strategy

%step 1, maintains b
b = 0;

%step 2, sumu
R = zeros(1, T);
Tau = zeros(1, T);
k = 0;
t = 0;
prev_delay = 0;
AoI_online = zeros(1, T);
Frame_online = zeros(1, T);
res = 0;
t_cum = 0;
prev_aoi = 0;
current_aoi = 0;
t_prev = 0;
prev_dpw = 0;
eta = 1.;
tot_aoi_prev = 0;
tot_frame_len = 0;
while k < T
    k = k + 1;
    current_delay = delay(k);
    w = max(b - current_delay, 0);
    R(k) = (2 * prev_dpw + current_delay) * current_delay / 2 + (2 * current_delay + w) * w / 2 + C;
    prev_dpw = current_delay + w;
    grad = (current_delay + w) ^ 2 / 2 - b * (current_delay + w) + C;
    b = max(b + eta * grad, 0);
    eta = 1 / k;
    tot_aoi_prev = tot_aoi_prev + R(k);
    tot_frame_len = tot_frame_len + prev_dpw;
    AoI_online(k) = tot_aoi_prev;
    Frame_online(k) = tot_frame_len;
end
avg_AoI_online = sum(R) / t;

t_rng = 1:T;
AoI_zw_dataset(exp_idx, 1:T) = AoI_zw(1:T) ./ Frame_zw(1:T);
AoI_opt_dataset(exp_idx, 1:T) = AoI_opt(1:T) ./ Frame_opt(1:T);
AoI_online_dataset(exp_idx, 1:T) = AoI_online(1:T);
AoI_online_l_dataset(exp_idx, 1:T) = Frame_online(1:T);
end
% figure;

AoI_zw_mean = mean(AoI_zw_dataset, 1);
AoI_opt_mean = mean(AoI_opt_dataset, 1);
AoI_online_tot = mean(AoI_online_dataset, 1);
AoI_online_l = mean(AoI_online_l_dataset, 1);
AoI_online_mean = AoI_online_tot(1:T) ./ AoI_online_l(1:T);
plot(1:T, mean(AoI_zw_mean) * ones(1, T), ':', 'linewidth', 2, 'Color', 'b');
hold on
% semilogx(1:K,lower_bound*ones(K,1));
% semilogx(1:K,Reward);
plot(1:T, mean(AoI_opt_mean) * ones(1, T), 'linewidth', 2, 'Color', 'g');
plot(1:T, AoI_online_mean, '-.', 'linewidth', 2, 'Color', 'r');
t_rng = 1:100:T;
AoI_online_ub = t_rng;
AoI_online_lb = t_rng;
% patch([t_rng fliplr(t_rng)], [AoI_online_mean(t_rng)-AoI_online_std(t_rng) fliplr(AoI_online_mean(t_rng)+AoI_online_std(t_rng))], 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
h = legend('$\pi_{\mathsf{zw}}$', '$\pi^*$', '$\pi_{\mathsf{online}}$', 'Interpreter', 'Latex');
set(h,'Fontsize',14);
xlabel('Cycle, $k$', 'Interpreter', 'Latex');
ylabel('Average Cost $h_{\pi}$', 'Interpreter', 'Latex');
ylim([0 150])