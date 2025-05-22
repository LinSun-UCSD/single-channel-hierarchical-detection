clear all;
clc;
close all;

% load data
data = "";
fs = 10000;
range = 34*fs : 37* fs;
% range = 1:length(data);
x = data(:, 4);
x = x - mean(x);
x = x';
t = (1:length(x)) / fs;
N = length(x);

%% get mu0, sigma0
temp = load("D:\UCSD Post-doc\paper\paper 1\figures\code\statistics\mean.mat");
temp = temp.monthlyValuesMean{1,5};
sigma0 = std(temp);
mu0 = mean(temp);
%% get the sigmA
indices = [];
t = (1:280*fs) / fs;
y = [];
indices = [];
len = fs;
T = [];
for i = 1:280
    clear index
    window = x((i-1)*len+1: (i) * len);

    window = window - mean(window);
    window = window ./std(window);
    shift = (i-1) *len;
    index = find(window > 7);

    indices = [indices index + shift];
    y = [y window];
end
sigma_A = std(y(indices));
%% plot H0, H1 distribution
close all;
figure;
mu = 0;
sigma = 1;
x = linspace(-4, 30, 1000);
pdf_values = normpdf(x, mu, sigma);
plot(x, pdf_values, 'LineWidth', 2);

hold on;

x = linspace(-4, 30, 1000);
pdf_values = normpdf(x, mu0, sqrt(sigma_A^2 + sigma0 ^ 2 + sigma^1));
plot(x, pdf_values, 'LineWidth', 2);
% xline(0, "color","r", "LineWidth", 2)
% xline(mu0, "color","r", "LineWidth", 2)
xlabel('x');
ylabel('Probability Density');
% title('PDF of Binary Hypothesis');
grid on;
set(gca, 'FontName', "Times New Roman", "FontSize", 12);
legend({"H_{0}", "H_{1}"})
print_plot("1.png", 5,3,800);
%% plot the mapping between x and z
close all;
temp_x = -10:1:10;
z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
figure()
plot(temp_x, z, "LineWidth", 2, "Color", "k");
grid on;
xlabel("x");
ylabel("z");
grid minor;
set(gca, "FontSize", 12, "FontName", "Times New Roman");
print_plot("1.png",5,2,800);
%% pfa

N = 10000000;
pfa = 0.000001;
temp_x = randn(N,1) *sigma;

z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
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
temp_x = randn(N,1) * sqrt((sigma^2 + sigma_A^2 + sigma0^2)) + mu0;

z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;

plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "b");
y = interp1(z_sorted, temp_pfas, z_value, "linear");
plot(z_value, y, 'pentagram', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b','HandleVisibility', 'off');
xlim([-10 100])
legend({"p_{FA}",  "p_{D}"}, "location", "SouthEast")
xlabel("z");
ylabel("Probability")
grid on;
set(gca, "FontSize", 12, "FontName", "Times New Roman");
print_plot("1.png",5,2,800);
%%

close all;

data = "";
fs = 10000;
range = 20*fs : 150 * fs;
% range = 1:length(data);
x = data(:, 4);
x = x - mean(x);
x = x';

pfas = [0.000001];
figure()
len = 1*fs;
N = 100000000;
for j = 1 : length(pfas)
    subplot(4, 1, j)
    pfa = pfas(j);
    indices = [];
    t = (1:296*fs) / fs;
    y = [];
    indices = [];
    T = [];
    temp_x = randn(N,1) *sigma;

    z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
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

        window_z = window.^2 / (sigma^2) - (window - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
        % window_z = window.^2 / sigma^2 - (window-mu).^2/(sigma^2 + sigma_A^2);
        clear index;
        y = [y window_z];
        % figure()
        % plot(window_z)
        index = find(window_z> z_value);
        shift = (i-1) *len;
        index = index + shift;
        indices = [indices index];
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
set(gca, 'FontName', "Times New Roman", "FontSize", 12);
% xlim([0 296])
print_plot("1.png",5,3,800)
%%
figure();
