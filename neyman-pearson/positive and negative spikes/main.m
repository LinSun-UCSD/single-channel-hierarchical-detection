clear all;
clc;
close all;

% load data

data = "";
fs = 10000;

% range = 156*fs : 250 * fs;
range = 1:length(data);
figure()
for i = 1:4
    subplot(4,1,i)
    plot(data(:,i))
end
x = data(:, 4);
x = x - mean(x);
%%

data = "";
fs = 10000;
range = 20*fs : 150 * fs;
% range = 1:length(data);
x = data(:, 4);
x = x - mean(x);
x = x';
A = 7;
pfas = [0.000001];
figure()
len = 1*fs;
N = 1000000;
transformed_z = [];
for j = 1 : length(pfas)
    subplot(2, 1, j)
    pfa = pfas(j);
    indices = [];
    t = (1:296*fs) / fs;
    y = [];
    indices = [];
    T = [];
    sigma = 1;
    temp_x = randn(N,1) *sigma;
    temp_x = abs(temp_x);
    z = -2 + log(exp(-(temp_x - A).^2 / 2 / sigma^2) + exp(-(temp_x + A).^2 / 2 / sigma^2)) ...
        + temp_x.^2 / 2 / sigma^2;
    z_sorted = sort(z);
    cdf_values = (1:N) / N;
    temp_pfas = 1-cdf_values;
    index = round((1-pfa) * N);
    z_value = z_sorted(index);
    for i = 1:296
        window = x((i-1)*len+1: (i) * len);
        window = window/std(window);
        sigma = sqrt(var(window));
        window = window - mean(window);

        temp_window = abs(window);
        window_z = -2 + log(exp(-(temp_window - A).^2 / 2 / sigma^2) + exp(-(temp_window + A).^2 / 2 / sigma^2)) ...
            + temp_window.^2 / 2 / sigma^2;
        clear index;
        transformed_z = [transformed_z window_z];
        % figure()
        % plot(window_z)
        index = find(window_z > z_value);
        shift = (i-1) *len;
        index = index + shift;
        indices = [indices index];
        y = [y x((i-1)*len+1: (i) * len)];

    end



end
%%
close all;
figure()
plot(t, transformed_z, 'b');
hold on;
yline(z_value, "Color", 'm', "LineWidth", 2)
plot(t(indices), transformed_z(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
xlabel("Time [sec]");
ylabel('z');
grid on;
set(gca, "FontSize", 12, "FontName", "Times New Roman");
print_plot("1.png", 5, 1.5, 800)
%%
close all;

figure()
temp_x = randn(N,1) *sigma;
temp_x = abs(temp_x);
z = -2 + log(exp(-(temp_x - A).^2 / 2 / sigma^2) + exp(-(temp_x + A).^2 / 2 / sigma^2)) ...
    + temp_x.^2 / 2 / sigma^2;
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);
z_value = z_sorted(index);
plot(z_sorted, temp_pfas, "LineWidth", 2, "Color", "k");
xlabel("z");
ylabel("Probability");
grid on;

hold on;

for A = 1:3:9
temp_x = randn(N,1) *sigma + A;
temp_x = abs(temp_x);
z = -2 + log(exp(-(temp_x - A).^2 / 2 / sigma^2) + exp(-(temp_x + A).^2 / 2 / sigma^2)) ...
    + temp_x.^2 / 2 / sigma^2;
z_sorted = sort(z);
cdf_values = (1:N) / N;
temp_pfas = 1-cdf_values;
index = round((1-pfa) * N);
z_value = z_sorted(index);
plot(z_sorted, temp_pfas, "LineWidth", 2);
end
set(gca, "FontSize", 12, "FontName", "Times New Roman");
lgd = legend({"p_{FA}", "p_{D} where \theta = 1", "p_{D} where \theta = 4", "p_{D} where \theta = 7"})

% Set legend background transparency (50% transparent)
lgd.FontSize = 9

print_plot("1.png", 5,3,800)
%%


% print_plot("1.png",12,6,800)
%%
figure();
plot(t, y, 'k')
hold on;
plot(t(indices), y(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
% plot(t, T, 'r');
xlim([t(1) t(end)])
xlabel("Time [sec]");
ylabel("Accel. [g]");
% title(['p_{fa} = ' num2str(pfa)])
grid on;
set(gca, 'FontName', "Times New Roman", "FontSize", 12);
print_plot("1.png", 5,1.5,800)