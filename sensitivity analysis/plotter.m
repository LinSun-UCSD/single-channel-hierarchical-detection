close all;
clear all;
clc;
figure()

subplot(3,1,1)
mu = -10:0.01:10;
hold on;
res = load("mean_sensitivity.mat").res1;
plot(mu, res(:,1)/(2*0.01), "LineWidth", 1.2, "Color", "r");
hold on;
plot(mu, res(:,2)/(2*0.01), "LineWidth", 1.2, "Color", "b");
legend({"\mu_{0}","\mu_{1}"})
ylim([-4000 4000]);
grid on;
xlabel("Normalized Mean");
ylabel("LSA");
box on;
set(gca, "FontSize", 10, "FontName", "Times New Roman");
ax = gca; % Get the current axes
ax.YAxis.Exponent = 3; % Force y-axis labels to show as x10^3
ax.YAxis.TickLabelFormat = '%.1f';

subplot(3,1,2)
std= 0.1:0.01:2;
res = load("std_sensitivity.mat").res2;
plot(std, res(:,1)/(2*0.001), "LineWidth", 1.2, "Color", "r");
hold on;
plot(std, res(:,2)/(2*0.001), "LineWidth", 1.2, "Color", "b");
grid on;
xlim([std(1) std(end)]);
ylim([-4000 4000]);
legend({"\sigma_{0}","\sigma_{1}"})
xlabel("Normalized STD");
ylabel("LSA");
set(gca, "FontSize", 10, "FontName", "Times New Roman");
ax = gca; % Get the current axes
ax.YAxis.Exponent = 3; % Force y-axis labels to show as x10^3
ax.YAxis.TickLabelFormat = '%.1f';

subplot(3,1,3)
gamma = -1:0.01:1;;
res = load("gamma_sensitivity.mat").res1;
plot(gamma, res/(2*0.001), "LineWidth", 1.2, "Color", "k");
ylim([-4000 4000]);
grid on;
xlabel("\rho");
ylabel("LSA");
set(gca, "FontSize", 10, "FontName", "Times New Roman");
ax = gca; % Get the current axes
ax.YAxis.Exponent = 3; % Force y-axis labels to show as x10^3
ax.YAxis.TickLabelFormat = '%.1f';

print_plot("1.png",6,6,800)