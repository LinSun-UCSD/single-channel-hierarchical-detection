clear all;
close all;
clc;


data = "";
Fs = 10000;
t = (1 : length(data)) / Fs;
%%
close all;
fs = Fs;
len = Fs;
x = data(:, 4);
t = (1:288*fs) / fs;
y = [];
figure()
count = 0;
for i = 1:60:300
    
    window = x((i-1)*len+1: (i) * len);
    window = window - mean(window);
    [f,xi] = ksdensity(window, -0.02:0.0001:0.02);
    plot(xi, f, "LineWidth", 1.2, "LineStyle", "-")
    % xlim([-0.1 0.2])
    count = count + 1;
    hold on;
end
legend({"1^{st} window", "60^{th} window", "120^{th} window", "180^{th} window", "240^{th} window"})
grid on;
xlabel("Accel. [g]");
ylabel("PDF");
set(gca, "FontSize", 12, "FontName", "Times New Roman"); 
print_plot("1.png", 8, 2.5, 800)
%%
figure;
imagesc(1:280, x, pdfs);  % Rows: x, Columns: PDFs
axis xy; % Correct axis orientation
colorbar;
xlabel('PDF Index');
ylabel('X-axis');
title('Heatmap of 280 PDFs');
