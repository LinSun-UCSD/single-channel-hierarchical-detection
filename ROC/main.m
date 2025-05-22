%% plot ROC
close all;
clear all;
clc;
figure()

sigma = 1;
pd = [];
pfa = 0.00:0.001:1;
A = 0.5;


pd = zeros(length(pfa), 1);
for i = 1:length(pfa)

    pd(i) =  0.5 * erfc(erfcinv(2*pfa(i)) - A / sqrt(2) / sigma);

end
plot(pfa, pd, 'b', 'LineWidth', 1.5);
auc1 = trapz(pfa, pd)
hold on;


% plot traditional Bayesian
N = 100000;
mu = A;
sigma_A = 1;
temp_x = randn(N,1) *sigma;
z = temp_x.^2 / sigma^2 - (temp_x-mu).^2/(sigma^2 + sigma_A^2);
z1 = sort(z);
cdf_values = (1:N) / N;
pfas = 1-cdf_values;


temp_x = randn(N,1) * sqrt(sigma_A^2 + sigma^2) + mu;
z = temp_x.^2 / sigma^2 - (temp_x-mu).^2/(sigma^2 + sigma_A^2);
z2 = sort(z);
cdf_values = (1:N) / N;
pds = 1-cdf_values;
pd = zeros(length(pfa), 1);
hold on;
for i = 1:length(pfa)-1
    index = round((1-pfa(i)) * N);
    z_value = z1(index);
    pd(i) = interp1(z2, pds, z_value, "linear", "extrap");

end
plot(pfa(1:length(pfa ) -1 ), pd(1:length(pfa ) -1 ), 'k', "LineWidth", 1.5)
auc2 = trapz(pfa(1:length(pfa ) -1 ), pd(1:length(pfa ) -1 ))
% plot the hierarhical Bayesian model

mu0 = A;
sigma_A = 1;
sigma0 = 1;
temp_x = randn(N,1) *sigma;
z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
z1 = sort(z);
cdf_values = (1:N) / N;
pfas = 1-cdf_values;


temp_x = randn(N,1) * sqrt(sigma^2 + sigma_A^2 + sigma0^2) + mu;
z = temp_x.^2 / (sigma^2) - (temp_x - mu0).^2 / (sigma^2 + sigma_A^2 + sigma0^2);
z2 = sort(z);
cdf_values = (1:N) / N;
pds = 1-cdf_values;
pd = zeros(length(pfa), 1);
hold on;
for i = 1:length(pfa)-1
    index = round((1-pfa(i)) * N);
    z_value = z1(index);
    pd(i) = interp1(z2, pds, z_value, "linear", "extrap");

end
plot(pfa(1:length(pfa ) -1 ), pd(1:length(pfa ) -1 ),'r', "LineWidth", 1.5);
plot([0 1], [0, 1], "LineStyle", "--", "Color", 'k', "LineWidth", 1.5)
legend({"Neyman-Pearson", "Traditional Bayesian", "Hierarchical Bayesian"}, 'Location', 'southeast');
xlim([0, 1]);
ylim([0, 1]);
grid on;
grid minor
xlabel("p_{FA}")
ylabel("p_{D}")
set(gca, "FontSize", 12, "FontName", "Times New Roman");
auc3 = trapz(pfa(1:length(pfa ) -1 ), pd(1:length(pfa ) -1 ))
print_plot("1.png", 4, 3, 800)