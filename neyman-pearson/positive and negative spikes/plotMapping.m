clear all;
close all;
%%
close all
% Define the range for A and x
A = linspace(2, 10, 20); % Example range for A
x = linspace(-5, 5, 20);  % Example range for x

% Define sigma
sigma = 1;

% Create a meshgrid for A and x
[X, A_grid] = meshgrid(x, A);

% Compute Z values using the provided formula
Z = -2 + log(exp(-(X - A_grid).^2 / (2 * sigma^2)) + exp(-(X + A_grid).^2 / (2 * sigma^2))) ...
    + X.^2 / (2 * sigma^2);

% Plot the surface
figure;
surf(A_grid, X, Z, 'EdgeColor', 'none'); % Create the surface plot
xlabel('\theta'); % Label for A-axis
ylabel('x'); % Label for x-axis
zlabel('z'); % Label for z-axis
% title('Surface Plot of z with respect to \theta and x');
colorbar; % Add a colorbar for visualization
set(gca, "FontSize",12, "FontName", "Times New Roman");
grid on; % Ensure grid lines are enabled
box on;
shading faceted; % Add grid lines on the surface
print_plot("1.png", 4, 3, 600)