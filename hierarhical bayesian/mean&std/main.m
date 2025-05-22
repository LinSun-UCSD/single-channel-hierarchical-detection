clc;
close all;
clear all;
%%
sigma = 1;
temp_x1 = "";
temp_x2 = "";
temp_x1 = temp_x1.monthlyValuesMean{1,5}';
temp_x2 = temp_x2.monthlyValuesSTD{1,5}';
data = [temp_x1 temp_x2];
% Number of samples to draw
n_samples = 1000;
cov_matrix = cov(temp_x1, temp_x2);
% Generate samples from the joint Gaussian distribution
samples = mvnrnd([mean(temp_x1); mean(temp_x2)], cov_matrix, n_samples);
%%

figure()
scatter(samples(:, 1), samples(:,2),"Marker", ".","MarkerEdgeColor", "blue");
grid on;
box on;
xlabel("Normalized Mean");
ylabel("Normalized STD");
% title("Samples from Joint Guassian Distribution");
xlim([5 15]);
ylim([0 10]);

set(gca, "FontSize", 10, "FontName", "Times New Roman");
hold on;
scatter(temp_x1, temp_x2,"Marker", ".","MarkerEdgeColor", "k");
grid minor;
box on;
xlabel("Normalized Mean");
ylabel("Normalized STD");
% title("Samples from Events")
% legend({"Joint Gaussian Distribution", "Samples from Population Probability"})
set(gca, "FontSize", 12, "FontName", "Times New Roman");
xlim([5 15]);
ylim([0 10]);
print_plot("1.png", 4,2.5,800);
%% get the mapping
close all;
temp_x = -10: 0.01: 10;
k = 100;
[temp_z, nominator, denominator, theta1] = getMapping(sigma, data, k, temp_x);
figure()
plot(temp_x, denominator, "LineWidth", 2);
hold on;
plot(temp_x, nominator, "LineWidth", 2)

hold on;
k = 100;
[temp_z, nominator, denominator, theta2] = getMapping(sigma, samples, k, temp_x);

hold on;
plot(temp_x, nominator, "LineWidth", 2, "LineStyle",'--');
grid on;
box on;
xlabel("x");
ylabel("Probability Density")
set(gca, "FontName", "Times New Roman", "FontSize", 12)
legend({"H_{0}", "H_{1} from samples of events", "H_{1} from sample of Gaussian Dist."}, "FontSize",  8)
print_plot("1.png", 4, 2.5, 800)
%% plot theta
figure()
[f,xi] = ksdensity(theta1);
plot(xi, f, "Color", "k", "LineWidth", 1.2)
hold on;
[f,xi] = ksdensity(theta2);
plot(xi, f, "Color", "b", "LineWidth", 1.2, "LineStyle","--" )
xlabel("\theta");
ylabel("Probability Density");
box on;
grid on;
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png", 4, 2.5, 800)

%% plot the mapping
figure()
plot(temp_x, temp_z, "LineWidth", 2, "Color",'k');
grid on;
xlabel("x");
ylabel("z");
set(gca, "FontName", "Times New Roman", "FontSize", 10);
% legend({"\theta based on ", "\theta based Events"})
box on;
print_plot("1.png", 4, 2.5, 800)

%% samples from H1 hypothesis
type = "numerical";
log_target_pdf = [temp_x' log(nominator)];
n_samples = 100000;
initial_state = [10];
proposal_std = 10;
n_dim = 1;
log_target_pdf = log_target_pdf';
[samples, acceptance_rate] = metropolis(log_target_pdf, n_samples, initial_state, proposal_std, n_dim, type);
%%
burn_in = 2000;
samples = samples(burn_in:end);
figure()
histogram(samples, "NumBins", 30);
xlabel("x;H_{1}");
ylabel("Count");
title("Samples of x under H_{1} from Metropolis Algorithm");
grid on;
box on;
set(gca, "FontSize", 10, "FontName", "Times New Roman");
print_plot("1.png", 4, 2.5, 800)
%%
pfa = 0.000001;
N = 10000000;
x = randn(N,1) *sigma;
z = zeros(N,1);
for i = 1:N
    z(i) = interp1(temp_x, temp_z, x(i), "Linear", "extrap");
end
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);
z_value = z_sorted(index);
figure()
plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "r");
hold on;
plot(z_sorted(index), temp_pfas(index), 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r','HandleVisibility', 'off');
N = length(samples);


z = zeros(N,1);
for i = 1:N
    z(i) = interp1(temp_x, temp_z, samples(i), "Linear", "extrap");
end
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "b");
[z_sorted_unique, idx] = unique(z_sorted, 'stable'); % Keep the original order
temp_pfas_unique = temp_pfas(idx); % Select corresponding elements
y = interp1(z_sorted_unique, temp_pfas_unique, z_value, "linear");
% y = interp1(z_sorted, temp_pfas, z_value, "linear");
plot(z_value, y, 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','HandleVisibility', 'off');
xlim([-10 50]);
xlabel("x");
ylabel("Probability");
legend({"p_{FA}", "p_{D}"});
grid on;
box on;
set(gca, "FontSize", 10, "FontName", "Times New Roman");
print_plot("1.png", 4,2.5, 800)
%%
close all;
path = 'D:\UCSD Post-doc\pythonCode\SooLocks\signal processing\00272281_00272575\'
data = readNPY([path 'ai211.npy']);
fs = 10000;
range = 20*fs : 150 * fs;
% range = 1:length(data);
x = data(:, 4);
x = x - mean(x);
x = x';

pfas = [0.000001];
figure()
len = 1*fs;
N = 100000;
for j = 1 : length(pfas)
    subplot(4, 1, j)
    pfa = pfas(j);
    indices = [];
    t = (1:296*fs) / fs;
    y = [];
    indices = [];
    T = [];
    pfa = 0.000001;
    x_noise = randn(N,1) *sigma;
    z = zeros(N,1);
    for i = 1:N
        z(i) = interp1(temp_x, temp_z, x_noise(i), "Linear", "extrap");
    end
    z_sorted = sort(z);
    cdf_values = (1:N) / N;
    temp_pfas = 1-cdf_values;
    index = round((1-pfa) * N);
    z_value = z_sorted(index);
    for i = 1:296
        window = x((i-1)*len+1: (i) * len);
        window = window - mean(window);
        
        window = window./std(window);
        % y = [y window];
        window_z = zeros(length(window), 1);

        for ii = 1:length(window)
            window_z(ii) = interp1(temp_x, temp_z, window(ii), "Linear", "extrap");
        end
        clear index;
        y = [y; window_z];
        index = find(window_z> z_value);
        shift = (i-1) *len;
        index = index + shift;
        indices = [indices; index];
    end
end
%%
close all;
figure()
subplot(2,1,1)
plot(t, y, 'color', 'b');
yline(z_value, "LineWidth", 2, "color", "m")
grid on;
hold on;
xlabel("Times [sec]");
ylabel("z");
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
plot(t(indices), y(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');

subplot(2,1,2)
plot(t, x, 'k')
hold on;
plot(t(indices), x(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
% plot(t, T, 'r');
xlim([t(1) t(end)])
xlabel("Time [sec]");
ylabel("Accel. [g]");
grid on;
set(gca, 'FontName', "Times New Roman", "FontSize", 10);
% xlim([0 296])
print_plot("1.png",5,3,800)

% temp_x = randn(N,1) * sqrt((sigma^2 + sigma_A^2 + sigma0^2)) + mu0;
% 
% z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
% z_sorted = sort(z);
% cdf_values = (1:N) / N;
% temp_pfas = 1-cdf_values;
% 
% plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "b");
% xlim([-10 100])
% legend({"p_{FA}", "p_{D}"})
% xlabel("z");
% ylabel("Probability")
% grid on;
% set(gca, "FontSize", 10, "FontName", "Times New Roman");
% print_plot("1.png",5,2,800);
%%

%%
function [temp_z, nominator, denominator, theta] = getMapping(sigma,samples, k, temp_x)
n = length(samples);
theta = zeros(k, n);
for i = 1:n
    muA = samples(i, 1);
    sigmaA = samples(i, 2);
    theta(:, i) = randn(k,1) * sigmaA + muA;

end
theta = reshape(theta,[],1);
% figure()
% histogram(theta)
nominator = zeros(length(temp_x), 1);
denominator = zeros(length(temp_x), 1);
temp_z = zeros(length(temp_x), 1);
for i = 1:length(temp_x)
    nominator(i) = mean((1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((temp_x(i) - theta) / sigma).^2));
    denominator(i) = mean((1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((temp_x(i)) / sigma).^2));
    temp_z(i) = log(nominator(i)) - log(denominator(i));
end

end