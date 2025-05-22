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
%%
close all;

fs = 10000;
figure()
pfas = [0.00001 0.000001]
len = fs;
for j = 1 : length(pfas)
    subplot(2, 1, j)
    pfa = pfas(j);
    indices = [];
    t = (1:296*fs) / fs;
    y = [];
    indices = [];
    T = [];
    for i = 1:296
        window = x((i-1)*len+1: (i) * len);
        y = [y; window];
        window = window ;
        sigma = sqrt(var(window));
        threshold = sqrt(2) * erfcinv(2 * pfa) * sigma ;
        index = find(window> threshold);
        shift = (i-1) *len;
        index = index + shift;
        indices = [indices; index];
        
        T = [T; threshold * ones(length(window), 1)];
    end

    plot(t, y, 'k')
    hold on;
    plot(t(indices), y(indices),'.', 'MarkerSize',10, 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r');
    plot(t, T, 'm');
    xlim([t(1) t(end)])
    xlabel("Time [sec]");
    ylabel("Accel. [g]");
    % title(['Pfa = ' num2str(pfa)])
    grid on;
    set(gca, 'FontName', "Times New Roman", "FontSize", 12);
    
end

print_plot("1.png", 6, 4, 800)

%% plot ROC
close all;
clear all;
clc;
figure()
A = 2;
sigma = 1;
pd = [];
for A = 1:1:4
    pfas = 0.00:0.001:1;
    pd = zeros(length(pfas), 1);
    for i = 1:length(pfas)
        pfa = pfas(i);
        pd(i) =  0.5 * erfc(erfcinv(2*pfa) - A / sqrt(2) / sigma);

    end
    plot(pfas, pd, 'b', 'LineWidth', 1.2);
    hold on;
end
grid on;
xlabel("p_{fa}");
ylabel("p_{d}");
title("ROC Curve");
set(gca, 'FontName', "Times New Roman", "FontSize", 12);
print_plot("1.png", 6,4, 800)
%% ROC
close all;
clear all;
clc;
figure()
sigma = 1;
pd = [];
A = 0:0.01:20;
pd = zeros(length(A),1);
pfa = 0.02;
for pfa = 0.000001:0.000005:0.001
for i = 1:length(A)

    pd(i) =  0.5 * erfc(erfcinv(2*pfa) - (A(i) / sqrt(2) / sigma));

end
plot(A, pd, 'b')
hold on;
end
grid on;
xlabel("\theta/\sigma");
ylabel("P_{d}");
title("ROC Curve");
set(gca, 'FontName', "Times New Roman", "FontSize", 12);
print_plot("1.png", 6,4, 800)
%%
% Define the range for A and PFA
A = linspace(1, 10, 20); % Example range for A
pfa = logspace(-6, -2, 20); % Example range for PFA (logarithmic scale)

% Preallocate PD matrix
pd = zeros(length(pfa), length(A));

% Define sigma
sigma = 1;

% Compute PD values for each combination of A and PFA
for i = 1:length(pfa)
    for j = 1:length(A)
        pd(i, j) = 0.5 * erfc(erfcinv(2 * pfa(i)) - (A(j) / sqrt(2) / sigma));
    end
end

% Create a meshgrid for 3D plotting
[X, Y] = meshgrid(A, pfa);

% Plot the 3D surface
figure;
surf(X, Y, pd, 'FaceColor', 'interp', 'EdgeColor', 'k'); % 'interp' for interpolated color, 'k' for black grid lines
xlabel('\theta/\sigma'); % Label for X-axis
ylabel('p_{FA}'); % Label for Y-axis
zlabel('p_{D}'); % Label for Z-axis
% title('3D Plot of PD as a function of A and PFA');
colorbar; % Optional: Add a colorbar for better visualization
set(gca, 'YScale', 'log'); % Use logarithmic scale for PFA

% Adjust appearance
grid on; % Ensure grid lines are enabled
shading faceted; % Add grid lines on the surface
view([8 31]);
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png", 5, 4, 800)