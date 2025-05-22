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
mu0 = mean(temp_x1);
mu1 = mean(temp_x2);
sigma0 = sqrt(cov_matrix(1,1));
sigma1 = sqrt(cov_matrix(2,2));
pfa = 0.000001;
gamma = cov_matrix(2,1) / sigma0 / sigma1;
N = 1000000;
k = 100;
temp_mu0 = -10:0.01:10;
temp_mu1 = -10:0.01:10;

res1 = zeros(length(temp_mu0),2);
temp_x = -10: 0.01: 10;
%% local sensitivy analysis of mu0
close all
delta = [0.01 -0.01];
values = zeros(2, 1);


for j = 1:length(temp_mu0)
    for ii = 1:2
        samples = mvnrnd([temp_mu0(j)+delta(ii); mu1;], cov_matrix, n_samples);
        [temp_z, nominator, denominator, theta2] = getMapping(sigma, samples, k, temp_x);
        x = randn(N,1) *sigma;
        z = zeros(N,1);
        for i = 1:N
            z(i) = interp1(temp_x, temp_z, x(i), "Linear", "extrap");
        end
        z_sorted = sort(z);
        cdf_values = (1:N) / N;
        temp_pfas = 1-cdf_values;
        index = round((1-pfa) * N);
        values(ii) = z_sorted(index);
    end
    res1(j,1) = values(1) - values(2);
end
%% local sensitivy analysis of mu1
close all
delta = [0.01 -0.01];
values = zeros(2, 1);

for j = 1:length(temp_mu1)
    for ii = 1:2
        samples = mvnrnd([mu0; temp_mu1(j)+delta(ii);], cov_matrix, n_samples);
        [temp_z, nominator, denominator, theta2] = getMapping(sigma, samples, k, temp_x);
        x = randn(N,1) *sigma;
        z = zeros(N,1);
        for i = 1:N
            z(i) = interp1(temp_x, temp_z, x(i), "Linear", "extrap");
        end
        z_sorted = sort(z);
        cdf_values = (1:N) / N;
        temp_pfas = 1-cdf_values;
        index = round((1-pfa) * N);
        values(ii) = z_sorted(index);
    end
    res1(j,2) = values(1) - values(2);
end
save("mean_sensitivity.mat", "res1")

%% utility function
function [temp_z, nominator, denominator, theta] = getMapping(sigma,samples, k, temp_x)
n = length(samples);
theta = zeros(k, n);
for i = 1:n
    muA = samples(i, 1);
    sigmaA = samples(i, 2);
    theta(:, i) = randn(k,1) * sigmaA + muA;

end
theta = reshape(theta,[],1);
nominator = zeros(length(temp_x), 1);
denominator = zeros(length(temp_x), 1);
temp_z = zeros(length(temp_x), 1);
for i = 1:length(temp_x)
    nominator(i) = mean((1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((temp_x(i) - theta) / sigma).^2));
    denominator(i) = mean((1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((temp_x(i)) / sigma).^2));
    temp_z(i) = log(nominator(i)) - log(denominator(i));
end

end