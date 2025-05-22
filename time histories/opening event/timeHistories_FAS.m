clear all;
close all;
clc;

data = ""
Fs = 10000;
t = (1 : length(data)) / Fs;
%%
titleName = {"Motor", "Speed Reducer Bearing Block",...
    "Jack Shaft Bearing Block", "Sector Gear Bearing Block"}
figure()
for i = 1:4
    subplot(4 ,1,i)
    plot(t, data(:,i),'Color', 'k');
    grid on;
    xlim([t(1) t(end)]);
    % xlabel("Time [sec]");
    % ylabel("Accel. [g]");
    title(titleName{i});
    set(gca, "FontName", "Times New Roman", "FontSize", 12);
end
print_plot("1.png", 6, 6.5, 800)
%%
figure()
inclinometer = readNPY([path 'ai217.npy']);
inclinometer = inclinometer(:,7);
plot(t, inclinometer,'Color', 'b', "LineWidth", 1.5);
grid on;
xlim([t(1) t(end)]);
% xlabel("Time [sec]");
% ylabel("Accel. [g]");
title("Sector Gear Inclinometer Time History");
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png", 6, 1.4, 800)
%%
x = data(:,4);
x = x - mean(x);
win_len = 256; % Window length
overlap = 128; % Overlap length
window = hann(win_len);

% Local mean removal for each window
x_demeaned = buffer(x, win_len, overlap, 'nodelay');
x_demeaned = x_demeaned - mean(x_demeaned, 1); % Subtract mean for each frame

% Reconstruct the signal after mean removal
x_clean = x_demeaned(:); % Flatten to 1D signal
[s, f, t_stft] = stft(x, Fs, 'Window', hann(256), 'OverlapLength', 128, 'FFTLength', 512*4);

% Normalize the magnitude spectrum
s_mag = abs(s);                 % Take the magnitude
s_norm = s_mag / max(s_mag(:)); % Normalize to range [0, 1]

% Plot normalized STFT magnitude spectrum
figure;
imagesc(t_stft, f, s_norm);     % Time vs frequency with normalized magnitude
axis xy; % Correct axis orientation
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 5000])
% title('Normalized STFT Magnitude Spectrum of Sector Gear Bearing Block');
colorbar;
caxis([0 0.1]);                 % Set color range from 0 to 0.5
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png", 6, 4, 800)
%%
close all
figure()
subplot(2,1,1)
plot(t, data(:, 4), 'k')
title("Sector Gear Bearing Block")
ylim([-0.3 0.3])
set(gca, "FontName", "Times New Roman", "FontSize", 12)
xlim([30 62]);
grid minor;

% ylabel("Accel. [g]")
subplot(2,1,2)
plot(t, inclinometer,'b', "LineWidth", 1);
xlim([30 62]);
title("Sector Gear Inclinometer Time History")
set(gca, "FontName", "Times New Roman", "FontSize", 12);
grid minor ;
% ylabel("Degree [^{o}]")
% ylabel("Accel. [g]")
print_plot("1.png",6,4,800);
%%
figure()
subplot(2,1,1)
plot(t, data(:, 4), 'k')
title("Sector Gear Bearing Block")
ylim([-0.3 0.3])
set(gca, "FontName", "Times New Roman", "FontSize", 12)
xlim([35 37]);
grid minor;

% ylabel("Accel. [g]")
subplot(2,1,2)
plot(t, inclinometer,'b', "LineWidth", 1);
xlim([35 37]);
title("Sector Gear Inclinometer Time History")
set(gca, "FontName", "Times New Roman", "FontSize", 12);
grid minor ;
% ylabel("Degree [^{o}]")
% ylabel("Accel. [g]")
print_plot("1.png",4,4,800);
