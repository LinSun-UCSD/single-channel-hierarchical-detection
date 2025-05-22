clc;
close all;
clear all;
%%
figure;
% Parameters for the folded normal distribution
mus = [0, 2]; % Mean of the normal distribution

for i = 1:length(mus)
    mu = mus(i);
    sigma = 1; % Standard deviation
    
    % Define the range of x
    x = linspace(0, 10, 1000);
    
    % Calculate the PDF of the folded normal distribution
    pdf = (1 / (sqrt(2 * pi) * sigma)) .* ...
        (exp(-((x - mu).^2) / (2 * sigma^2)) + exp(-((x + mu).^2) / (2 * sigma^2)));
    
    % Plot the PDF

        plot(x, pdf, 'LineWidth', 2);
   
    % title('Folded Normal Distribution');
    xlabel('x');
    ylabel('PDF');
    grid on;
    grid minor;
    
    hold on
end
% legend({['A = ' num2str(0) ', \sigma = ' num2str(sigma)], ['\theta = ' num2str(2) ', \sigma = ' num2str(sigma)], ['A = ' num2str(-2) ', \sigma = ' num2str(sigma)]});
legend({"H_{0}", "H_{1}"})
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png", 4, 3, 600)

%%
