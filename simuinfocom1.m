clear all 
clc
T = 1e6;
sigma = 1;
mu = 1;
exp_num = 20;
AoI_zw_dataset = zeros(exp_num, T);
AoI_opt_dataset = zeros(exp_num, T);
AoI_online_dataset = zeros(exp_num, T);

for exp_idx = 1:exp_num

%% generate delay distribution
pd = makedist('Lognormal','mu',log(mu),'sigma',sigma)
delay = random(pd,T,1);
exp_idx

%% zero wait
R = zeros(1, T);
Tau = zeros(1, T);
k = 0;
t = 0;
prev_delay = 0;
AoI_zw = zeros(1, T);
res = 0;
t_cum = 0;
prev_aoi = 0;
current_aoi = 0;
t_prev = 0;
while t < T - 1
    k = k + 1;
    current_delay = delay(k);
    t_prev = t;
    t = t + current_delay;
    if t_cum > t - 1
        res = res + (2 * prev_aoi + current_delay) * current_delay / 2;
        prev_aoi = current_delay;
    else
        t_cum = t_cum + 1;
        duration = t_cum - t_prev;
        AoI_zw(t_cum) = res + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = prev_aoi + duration;
        while t_cum < t - 1 && t_cum < T
            t_cum = t_cum + 1;
            AoI_zw(t_cum) = prev_aoi + 0.5;
            prev_aoi = prev_aoi + 1;
        end
        duration = t - t_cum;
        res = (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = current_delay;
    end
    R(k) = (2 * prev_delay + current_delay) * current_delay / 2;
    prev_delay = current_delay;
end
avg_AoI_zw = sum(R) / t;
cum_AoI_zw = sum(AoI_zw) / t;

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
    term2 = (exp(mu_bar + sigma^2 / 2))^2 + (exp(sigma^2) - 1) * exp( 2 * mu_bar + sigma^2 ) + b^2 * logncdf(b, mu_bar, sigma) - integral(logn_x_2, 0, b);
    if term1 - term2 / 2 / b > 0
        r = b;
    else
        l = b;
    end
end
R = zeros(1, T);
Tau = zeros(1, T);
k = 0;
t = 0;
prev_delay = 0;
AoI_opt = zeros(1, T);
res = 0;
t_cum = 0;
prev_aoi = 0;
current_aoi = 0;
t_prev = 0;
prev_dpw = 0;
while t < T
    k = k + 1;
    current_delay = delay(k);
    w = max(b - current_delay, 0);
    R(k) = (2 * prev_dpw + current_delay) * current_delay / 2 + (2 * current_delay + w) * w / 2;
    % R(k)
    prev_dpw = current_delay + w;
    t_prev = t;
    
    % step 2.1, add delay
    cycle_rwd = 0;
    t = t + current_delay;
    if t_cum > t - 1
        res = res + (2 * prev_aoi + current_delay) * current_delay / 2;
        cycle_rwd = cycle_rwd + (2 * prev_aoi + current_delay) * current_delay / 2;
        prev_aoi = current_delay;
        t_prev = t;
    else
        t_cum = t_cum + 1;
        duration = t_cum - t_prev;
        AoI_opt(t_cum) = res + (2 * prev_aoi + duration) * duration / 2;
        cycle_rwd = cycle_rwd + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = prev_aoi + duration;
        while t_cum < t - 1 && t_cum < T
            t_cum = t_cum + 1;
            AoI_opt(t_cum) = prev_aoi + 0.5;
            cycle_rwd = cycle_rwd + (2 * prev_aoi + duration) * duration / 2;
            prev_aoi = prev_aoi + 1;
            
        end
        duration = t - t_cum;
        res = (2 * prev_aoi + duration) * duration / 2;
        cycle_rwd = cycle_rwd + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = current_delay;
        t_prev = t;
    end
    
    % step 2.2 add waiting time
    t = t + w;
    if t_cum > t - 1
        res = res + (2 * current_delay + w) * w / 2;
        cycle_rwd = cycle_rwd + (2 * current_delay + w) * w / 2;
        prev_aoi = current_delay + w;
    else
        t_cum = t_cum + 1;
        duration = t_cum - t_prev;
        AoI_opt(t_cum) = res + (2 * prev_aoi + duration) * duration / 2;
        cycle_rwd = cycle_rwd + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = prev_aoi + duration;
        while t_cum < t - 1 && t_cum < T
            t_cum = t_cum + 1;
            AoI_opt(t_cum) = prev_aoi + 0.5;
            cycle_rwd = cycle_rwd + prev_aoi + 0.5;
            prev_aoi = prev_aoi + 1;
            
        end
        duration = t - t_cum;
        res = (2 * prev_aoi + duration) * duration / 2;
        cycle_rwd = cycle_rwd + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = current_delay + w;
    end
    % cycle_rwd
end
avg_AoI_opt = sum(R) / t;
cum_AoI_opt = sum(AoI_opt) / t;
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
res = 0;
t_cum = 0;
prev_aoi = 0;
current_aoi = 0;
t_prev = 0;
prev_dpw = 0;
eta = 1.;
while t < T
    k = k + 1;
    current_delay = delay(k);
    w = max(b - current_delay, 0);
    R(k) = (2 * prev_dpw + current_delay) * current_delay / 2 + (2 * current_delay + w) * w / 2;
    prev_dpw = current_delay + w;
    t_prev = t;
    
    % step 2.1, add delay
    t = t + current_delay;
    if t_cum > t - 1
        res = res + (2 * prev_aoi + current_delay) * current_delay / 2;
        prev_aoi = current_delay;
        t_prev = t;
    else
        t_cum = t_cum + 1;
        duration = t_cum - t_prev;
        AoI_online(t_cum) = res + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = prev_aoi + duration;
        while t_cum < t - 1 && t_cum < T
            t_cum = t_cum + 1;
            AoI_online(t_cum) = prev_aoi + 0.5;
            prev_aoi = prev_aoi + 1;
        end
        duration = t - t_cum;
        res = (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = current_delay;
        t_prev = t;
    end
    
    % step 2.2 add waiting time
    t = t + w;
    if t_cum > t - 1
        res = res + (2 * current_delay + w) * w / 2;
        prev_aoi = current_delay + w;
    else
        t_cum = t_cum + 1;
        duration = t_cum - t_prev;
        AoI_online(t_cum) = res + (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = prev_aoi + duration;
        while t_cum < t - 1 && t_cum < T
            t_cum = t_cum + 1;
            AoI_online(t_cum) = prev_aoi + 0.5;
            prev_aoi = prev_aoi + 1;
        end
        duration = t - t_cum;
        res = (2 * prev_aoi + duration) * duration / 2;
        prev_aoi = current_delay + w;
    end
    grad = (current_delay + w) ^ 2 - b * (current_delay + w) * 2;
    b = max(b + eta * grad, 0);
    eta = 1 / k;
end
avg_AoI_online = sum(R) / t;
cum_AoI_online = sum(AoI_online) / t;

t_rng = 1:T;
for t=2:T
    AoI_zw(t) = AoI_zw(t-1) + AoI_zw(t);
    AoI_opt(t) = AoI_opt(t-1) + AoI_opt(t);
    AoI_online(t) = AoI_online(t-1) + AoI_online(t);
end
AoI_zw_dataset(exp_idx, 1:T) = AoI_zw(1:T) ./ t_rng;
AoI_opt_dataset(exp_idx, 1:T) = AoI_opt(1:T) ./ t_rng;
AoI_online_dataset(exp_idx, 1:T) = AoI_online(1:T) ./ t_rng;
end
% figure;
subplot(1, 3, 1)

AoI_zw_mean = mean(AoI_zw_dataset, 1);
AoI_opt_mean = mean(AoI_opt_dataset, 1);
AoI_online_mean = mean(AoI_online_dataset, 1);
AoI_online_std = std(AoI_online_dataset, 1);
plot(1:T, mean(AoI_zw_mean) * ones(1, T),':',  'linewidth', 2, 'Color', 'b');
hold on
% semilogx(1:K,lower_bound*ones(K,1));
% semilogx(1:K,Reward);
plot(1:T, mean(AoI_opt_mean) * ones(1, T), 'linewidth', 2, 'Color', 'g');
plot(1:T, AoI_online_mean, '-.', 'linewidth', 2, 'Color', 'r');
t_rng = 1:100:T;
AoI_online_ub = t_rng;
AoI_online_lb = t_rng;
patch([t_rng fliplr(t_rng)], [AoI_online_mean(t_rng)-AoI_online_std(t_rng) fliplr(AoI_online_mean(t_rng)+AoI_online_std(t_rng))], 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
h = legend('$\pi_{\mathsf{zw}}$', '$\pi^*$', '$\pi_{\mathsf{online}}$', 'Interpreter', 'Latex');
set(h,'Fontsize',14);
xlabel('Time, $t$', 'Interpreter', 'Latex');
ylabel('Average AoI, $\overline{A}_{\pi,t}$', 'Interpreter', 'Latex');
ylim([0 30])
