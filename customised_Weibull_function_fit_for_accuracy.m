close all; clear
filePath = matlab.desktop.editor.getActiveFilename;
for i = length(filePath):-1:1
    if filePath(i) == filesep
        slash(i,1) = 1;
    end
end
filePath = filePath(1:length(slash) ); cd (filePath);   % change to code's path to load the dataset

load([filePath filesep 'dataset and figures' filesep 'dataset.mat']) 
%% defining variables
myfittype = fittype('1 - 0.5*exp(-(x/a)^b)',...         % Defining a customised Weibull function (1 - 0.5*exp(-(x/a)^b)')
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b'});
case_number = {'i', 'ii', 'iii', 'iv'};
num_of_cases = length(case_number);
Mo_Strength=[0:0.01:0.03 0.032 0.04:0.01:0.12 0.128 0.13:0.01:0.25 0.256 0.26:0.01:1]; % the the motion strength
x=100*Mo_Strength;
fit_lower_band = [0 0];

figure1 = figure('units','normalized','outerposition',[0 0 1 1]);
for accuracy_var_num = 1 : num_of_cases
    ac_temp = eval(['ac' num2str(accuracy_var_num)]);               % reading accuracy variables for different cases
    myfit{accuracy_var_num} = fit(x',ac_temp',myfittype, 'Lower', fit_lower_band);           % fitting accuracy data to the custom Weibull function
    % plotting the data and fitted curve
    fit_subplot(accuracy_var_num) = subplot(num_of_cases, 1 , accuracy_var_num);
    plot(myfit{accuracy_var_num},x,ac_temp)
    hold on
    Spacing_lines = 3;
    h = plot(nan(size(x,1),Spacing_lines));
    hold off
    set(h,{'Color'},{'k'}); % clear the dummy lines
    title(['Case ' case_number{accuracy_var_num}])
    xlabel('Motion strength (% coh)')
    ylabel('Accuracy')
    % place the legend:
    hl = legend([{'data','fitted curve'}]);
    % add your text:
    annotation('textbox',hl.Position - [0.08 0.01  0 0],'String',{['Scale (a) = ' num2str(myfit{accuracy_var_num}.a)];['Shape (b) = ' num2str(myfit{accuracy_var_num}.b)]},...
        'VerticalAlignment','top','Edgecolor','none', 'FitBoxToText', 'on');
%     ac_temp = [];
end
sgtitle(figure1,'Fitted psychometric functions to a customised Weilbull function (1 - 0.5*exp(-(c/a)^b))')

%% saving figures
figure_save_directory_name = 'Figures';
data_save_path = [filePath filesep 'dataset and figures'];
figure_save_path = [data_save_path filesep figure_save_directory_name];

% Clear unnecessary plot variables
clear *plot* *legend* han *lh ht

% creating folders for saving data
if ~isfolder(data_save_path)
    mkdir(data_save_path);
end

if ~isfolder(figure_save_path)
    mkdir(figure_save_path);
end

figure_name = 'customised_Weibull_fit';
exportgraphics(figure1,[fullfile(figure_save_path, figure_name) '.png'],'Resolution',300);
savefig(figure1,fullfile(figure_save_path, figure_name));