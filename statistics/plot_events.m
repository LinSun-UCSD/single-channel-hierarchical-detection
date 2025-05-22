close all;
clear all;
clc;


%% read the directory folders


% Get the list of all items in the directory
allItems = dir(path);

% Filter only directories (ignoring files)
isFolder = [allItems.isdir];  % Logical array indicating folders
folderNames = {allItems(isFolder).name};  % Get names of folders

% Remove '.' and '..' (current and parent directories)
folderNames = folderNames(~ismember(folderNames, {'.', '..'}));

% Display the folder names
disp('Folders:');
disp(folderNames);
%%
colors = {'k', 'b', 'r','c','green'};
legendNames = {"May", "June", "July", "Sept.", "Nov."};
monthlyValuesMean = cell(1, 12); % 12 months
monthlyValuesSTD = cell(1, 12); % 12 months
figure()
for i = 1:5
    norm_window1_mean = [];
    norm_window1_std = [];
    norm_window2_mean = [];
    norm_window2_std = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    
    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;
        datetimeString = dataToStore.time_stamp;

        % Convert to datetime object
        dt = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd''T''HH_mm_ss');
        
        % Get the month name
        monthNumber = month(dt);

        % Display the result
        disp(['Month Number: ', num2str(monthNumber)]);
        
       
        if strcmp(dataToStore.event_label, 'full-open')
            for l = 1:length(dataToStore.spikeTime)
                clear temp;
                if iscell(dataToStore.data1)
                    if length(dataToStore.spikeTime) > 1
                        temp1 = dataToStore.data1(1,l);
                        temp1 = temp1{1};
                        temp2 = dataToStore.data2(1,l);
                        temp2 = temp2{1};
                    else
                        temp1 = dataToStore.data1;
                        temp2 = dataToStore.data2;
                    end
                    positiveValues = temp1(temp1 > 0);
                    if length(positiveValues) == 0
                        continue;
                    end
                    % norm_window1_mean = [norm_window1_mean; mean(positiveValues)];
                    % norm_window1_std = [norm_window1_std; std(positiveValues)];
                    monthlyValuesSTD{monthNumber} = [monthlyValuesSTD{monthNumber}, std(positiveValues)];
                    monthlyValuesMean{monthNumber} = [monthlyValuesMean{monthNumber}, mean(positiveValues)];
                    % positiveValues = temp2(temp2 > 0);
                    % norm_window2_mean = [norm_window2_mean; mean(positiveValues)];
                  
                end
            end
        end
    end


end
%%
figure()
count = 1;
close all
for i = 1:12
    
    if (length(monthlyValuesMean{1,i}) > 1)
        subplot(1,5,count);
        scatter(monthlyValuesMean{1,i}, monthlyValuesSTD{1,i}, '.', 'MarkerEdgeColor', colors{count})
        count = count + 1;
        box on;
        grid on;
        xlabel("Normalized Mean");
        ylabel("Normalized STD");
        set(gca, "FontSize", 12, "FontName", "Times New Roman");
        xlim([5 15]);
        ylim([0 10])
    end
end
print_plot("1.png",12,2,800)
%%
figure()
count = 1;
close all
legendNames = {"May", "June", "July", "Sept.", "Nov."};
for i = 1:12
    
    if (length(monthlyValuesMean{1,i}) > 1)
        subplot(1,5,count);
        h = histogram(monthlyValuesMean{1,i}, 'NumBins', 20, "Normalization", "PDF");
        

        % Change the bin color
        h.FaceColor = colors{count}; % Change to red ('r')
        h.EdgeColor = 'k'; % Optional: Change the edge color to black ('k')
        count = count + 1;
        box on;
        grid on;
        xlabel("Normalized Mean");
        ylabel("p");
        
        set(gca, "FontSize", 12, "FontName", "Times New Roman");
        
    end
end
print_plot("1.png",12,2,800)
%%

colors = {'k', 'b', 'r','c','green'};
legendNames = {"May", "June", "July", "Sept.", "Nov."};
monthlyValuesOpen = zeros(12,1); 
monthlyValuesTotal = zeros(12,1);
figure()
for i = 1:5
    norm_window1_mean = [];
    norm_window1_std = [];
    norm_window2_mean = [];
    norm_window2_std = [];
    folderPath = [path '\' folderNames{i}];
    % Get a list of all .mat files in the folder
    matFiles = dir(fullfile(folderPath, '*.mat'));
    
    % Loop through each file and load the data
    for k = 1:length(matFiles)
        % Get the full path of the .mat file
        matFilePath = fullfile(folderPath, matFiles(k).name);

        % Load the .mat file
        data = load(matFilePath);
        dataToStore = data.dataToStore;
        datetimeString = dataToStore.time_stamp;

        % Convert to datetime object
        dt = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd''T''HH_mm_ss');
        
        % Get the month name
        monthNumber = month(dt);
        if strcmp(dataToStore.event_label, 'full-open') || strcmp(dataToStore.event_label, 'open') 
            monthlyValuesOpen(monthNumber) = monthlyValuesOpen(monthNumber)  + 1;

        end
         monthlyValuesTotal(monthNumber) = monthlyValuesTotal(monthNumber)  + 1;
    end


end 
save("mean.mat", "monthlyValuesMean");
save("std.mat", "monthlyValuesSTD")
%%
clc;
close all;
figure()
index = monthlyValuesTotal > 0;
plot(monthlyValuesTotal(index), "LineWidth",2, "Color", "k", "Marker", "s");
hold on;
plot(monthlyValuesOpen(index), "LineWidth",2, "Color", "r", "Marker", "h");
plot(monthlyValuesTotal(index) - monthlyValuesOpen(index), "LineWidth",2, "Color", "b", "Marker", "+");
legend({"Total", "Opening", "Closing"})
xticks(1:5); % Set x-tick positions
xticklabels(["May", "June", "July", "September", "November"]); % Set corresponding labels
grid on;
box on;
xlabel("Month");
ylabel("Count");
title("The Number of Events Over Time")
set(gca, "FontName", "Times New Roman", "FontSize", 12);
print_plot("1.png",5,3,800)
