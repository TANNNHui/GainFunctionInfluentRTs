close all; clear all;

scratch_main_dir = pwd; % Specify the root directory to save data and figures in (Current Folder)
tic;
%Parameter setting
z1= 1; % Upper decision threshold (choice 1)
z2= -1; % Lower decision threshol (choice 2)
dt=0.1; % Time step 5ms
name_of_cases = {'time-variant gain on both drift rate and noise term with noise depending on c', ...
    'time-variant gain on both drift rate and noise term with noise independent of c', ...
    'time-variant gain only on noise term', 'time-variant gain only on drift rate'};
trials=1000000; % Trial number can be reduced to 100

error_trial_number_threshold = floor(trials/1000) + 1; % mean reaction time of error trials gets equaled to NaN if the number of error trials is below threshold in each motion strength

interpolation_method = 'linear';    % interpolation method for subtituting zero values in mean reaction times of error trials up until the last non-zero element if needed

k=9.66*10^(-3); % the proportionality factor to form the mean of the drift rate
Mo_Strength=[0:0.01:0.03 0.032 0.04:0.01:0.12 0.128 0.13:0.01:0.25 0.256 0.26:0.01:1]; % the the motion strength
ti=0:dt:2000; % Time < 2s
Sy=7.34; % when G both: 7.34, sigle 0.734
Sx=0.0028;
d=584; % by ms
G=Sy.*exp(Sx*(ti-d))./(1+exp(Sx*(ti-d)))+(1+(1-Sy)*exp(-Sx*d))/(1+exp(-Sx*d));
sigma0=0.0188; % Size of the noise 

T_res=261;  % Delay between decision time and reaction time based on t_residual of Model 4 in (Ditterich, 2006)

AveDTc_initial=zeros(1,length(Mo_Strength));
AveDTe_initial=zeros(1,length(Mo_Strength));
ac1=zeros(1,length(Mo_Strength));

AveDTc_both=zeros(1,length(Mo_Strength));
AveDTe_both=zeros(1,length(Mo_Strength));
ac2=zeros(1,length(Mo_Strength));

AveDTc_onlyNoise=zeros(1,length(Mo_Strength));
AveDTe_onlyNoise=zeros(1,length(Mo_Strength));
ac3=zeros(1,length(Mo_Strength));

AveDTc_onlyDrift=zeros(1,length(Mo_Strength));
AveDTe_onlyDrift=zeros(1,length(Mo_Strength));
ac4=zeros(1,length(Mo_Strength));

AveMu=zeros(1,length(Mo_Strength));
AveNoi=zeros(1,length(Mo_Strength));
%% time-variant gain on both drift rate and noise term with noise depending on c
disp(['Analysing time-variant gain on both drift rate and noise term with noise depending on c at ' char(datetime)])

j=1;

for c=Mo_Strength
    DTc1=zeros(1,trials); % Correct reaction time
    DTe1=zeros(1,trials); % Error reaction time 
    Mu0=k*c;% Drift rate
    sigma=sigma0*sqrt(1+c); 
    parfor i=1:trials % Trial number
        rng(i); % Different random seed for different Trials but the seeds are the same all four scenarios
        x1=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x1(t+1)=x1(t) + dt*Mu0*G(t) + sqrt(dt)*sigma*G(t)*randn;% Update x;
            if x1(t) >= 1 % Record the correct reaction time
                DTc1(i)=ti(t)+T_res;
                break;
            end
            if x1(t) <= -1 % Record the error reaction time
                DTe1(i)=ti(t)+T_res;
                break;
            end
        end
%         xtemp(j, i,:) = x1;
    end
    DTc1(DTc1 < T_res)=[]; % Removing trials with reaction time below T_res
    DTe1(DTe1 < T_res)=[]; % Removing trials with reaction time below T_res     
    if(c==0.032)
        DTc1_032=DTc1;
        DTe1_032=DTe1;
        DTc1_032(DTc1_032==T_res)=[];
        DTe1_032(DTe1_032==T_res)=[];
    end
    if(c==0.128)
        DTc1_128=DTc1;
        DTe1_128=DTe1;
        DTc1_128(DTc1_128==0)=[];
        DTe1_128(DTe1_128==0)=[];      
    end
    if(c==0.256)
        DTc1_256=DTc1;
        DTe1_256=DTe1;
        DTc1_256(DTc1_256==0)=[];
        DTe1_256(DTe1_256==0)=[];      
    end
    AveDTc_initial(j)=mean(DTc1); %Calculate the average reaction time
    correct_trial_number(1,j) = length(find(DTc1 > 0)); % storing the number of correct trials in each case
    error_trial_number(1,j) = length(find(DTe1 > 0)); % storing the number of error trials in each case
    
    if error_trial_number(1,j) < error_trial_number_threshold                    
        AveDTe_initial(j) = 0;      % Assigning zeto to mean reaction time if there is no error trials
    else
        AveDTe_initial(j)=mean(DTe1);
    end


    ac1(j)=length(find(DTc1 > 0))/(length(find(DTc1 > 0)) + length(find(DTe1 > 0))); % Changed this to just include non-zero reaction times
    
    j=j+1;
    
end
% AveDTe_initial_temp = AveDTe_initial;
% AveDTe_initial = interpolate_vector(AveDTe_initial,interpolation_method);       % Interpolating zero values up until the last non-zero element
AveDTe_initial(AveDTe_initial==0) = NaN;
%% time-variant gain on both drift rate and noise term with noise independent of c
disp(['Analysing time-variant gain on both drift rate and noise term with noise independent of c at ' char(datetime)])

j=1;
for c=Mo_Strength
    DTc2=zeros(1,trials); % Correct reaction time
    DTe2=zeros(1,trials); % Error reaction time 
    temp1=0;
    temp2=0;
    Mu0=k*c;% Drift rate
    parfor i=1:trials % Trial number
        rng(i);
        x2=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            
            x2(t+1)=x2(t) + dt*Mu0*G(t) + sqrt(dt)*sigma0*G(t)*randn;% Update x;
            if x2(t) >= 1 % Record the correct reaction time
                DTc2(i)=ti(t)+T_res;
                temp1=trapz(Mu0*G(1:t))/t+temp1;
                temp2=trapz(sigma0*G(1:t))/t+temp2;
                break;
            end
            if x2(t) <= -1 % Record the error reaction time
                DTe2(i)=ti(t)+T_res;
                temp1=trapz(Mu0*G(1:t))/t+temp1;
                temp2=trapz(sigma0*G(1:t))/t+temp2;
                break;
            end
        end
    end
    DTc2(DTc2 < T_res)=[]; % Removing trials with reaction time below T_res
    DTe2(DTe2 < T_res)=[]; % Removing trials with reaction time below T_res
    if(c==0.032)
        DTc2_032=DTc2;
        DTe2_032=DTe2;
        DTc2_032(DTc2_032==0)=[];
        DTe2_032(DTe2_032==0)=[]; 
    end
    if(c==0.128)
        DTc2_128=DTc2;
        DTe2_128=DTe2;
        DTc2_128(DTc2_128==0)=[];
        DTe2_128(DTe2_128==0)=[]; 
    end
    if(c==0.256)
        DTc2_256=DTc2;
        DTe2_256=DTe2;
        DTc2_256(DTc2_256==0)=[];
        DTe2_256(DTe2_256==0)=[];      
    end
    AveDTc_both(j)=mean(DTc2); %Calculate the average reaction time

    correct_trial_number(2,j) = length(find(DTc2 > 0)); % storing the number of correct trials in each case
    error_trial_number(2,j) = length(find(DTe2 > 0)); % storing the number of error trials in each case

    if error_trial_number(2,j) < error_trial_number_threshold                   
        AveDTe_both(j) = 0;      % Assigning zeto to mean reaction time if there is no error trials
    else
        AveDTe_both(j)=mean(DTe2);
    end

    ac2(j)= length(find(DTc2 > 0))/(length(find(DTc2 > 0)) + length(find(DTe2 > 0))); % Changed this to just include non-zero reaction times
    AveMu(j)=temp1/trials;
    AveNoi(j)=temp2/trials;
    j=j+1;
end

% AveDTe_both_temp = AveDTe_both;
% AveDTe_both = interpolate_vector(AveDTe_both,interpolation_method); % Interpolating zero values up until the last non-zero element
AveDTe_both(AveDTe_both==0) = NaN;
%% time-variant gain only on noise term
disp(['Analysing time-variant gain only on noise term at ' char(datetime)])

j=1;
for c=Mo_Strength
    DTc3=zeros(1,trials); % Correct reaction time
    DTe3=zeros(1,trials); % Error reaction time 
    Mu0=AveMu(j);% Drift rate
    parfor i=1:trials % Trial number
        rng(i);
        x3=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x3(t+1)=x3(t) + dt*Mu0 + sqrt(dt)*G(t)*sigma0*randn;% Update x;
            if x3(t) >= 1 % Record the correct reaction time
                DTc3(i)=ti(t)+T_res;
                break;
            end
            if x3(t) <= -1 % Record the error reaction time
                DTe3(i)=ti(t)+T_res;
                break;
            end
        end
    end
    DTc3(DTc3 < T_res)=[]; % Removing trials with reaction time below T_res
    DTe3(DTe3 < T_res)=[]; % Removing trials with reaction time below T_res
    if(c==0.032)
        DTc3_032=DTc3;
        DTe3_032=DTe3;
        DTc3_032(DTc3_032==0)=[];
        DTe3_032(DTe3_032==0)=[]; 
    end
    if(c==0.128)
        DTc3_128=DTc3;
        DTe3_128=DTe3;
        DTc3_128(DTc3_128==0)=[];
        DTe3_128(DTe3_128==0)=[]; 
    end
    if(c==0.256)
        DTc3_256=DTc3;
        DTe3_256=DTe3;
        DTc3_256(DTc3_256==0)=[];
        DTe3_256(DTe3_256==0)=[];      
    end
    AveDTc_onlyNoise(j)=mean(DTc3); %Calculate the average reaction time

    correct_trial_number(3,j) = length(find(DTc3 > 0)); % storing the number of correct trials in each case
    error_trial_number(3,j) = length(find(DTe3 > 0)); % storing the number of error trials in each case

    if error_trial_number(3,j) < error_trial_number_threshold                  
        AveDTe_onlyNoise(j) = 0;      % Assigning zeto to mean reaction time if there is no error trials
    else
        AveDTe_onlyNoise(j)=mean(DTe3);
    end

%     AveDTe_onlyNoise(j)=dt.*mean(DTe3);
    ac3(j)=length(find(DTc3 > 0))/(length(find(DTc3 > 0)) + length(find(DTe3 > 0))); % Changed this to just include non-zero reaction times
    j=j+1;
end

% AveDTe_onlyNoise_temp = AveDTe_onlyNoise;
% AveDTe_onlyNoise = interpolate_vector(AveDTe_onlyNoise,interpolation_method);   % Interpolating zero values up until the last non-zero element
AveDTe_onlyNoise(AveDTe_onlyNoise==0) = NaN;
%% time-variant gain only on drift rate
disp(['Analysing time-variant gain only on drift rate at ' char(datetime)])

j=1;
for c=Mo_Strength
    DTc4=zeros(1,trials); % Correct reaction time
    DTe4=zeros(1,trials); % Error reaction time 
    Mu0=k*c;% Drift rate
    sigma0=AveNoi(j);
    parfor i=1:trials % Trial number
        rng(i);
        x4=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x4(t+1)=x4(t) + dt*Mu0*G(t) + sqrt(dt)*sigma0*randn;% Update x;
            if x4(t) >= 1 % Record the correct reaction time
                DTc4(i)=ti(t)+T_res;
                break;
            end
            if x4(t) <= -1 % Record the error reaction time
                DTe4(i)=ti(t)+T_res;
                break;
            end
        end
    end
    DTc4(DTc4 < T_res)=[]; % Removing trials with reaction time below T_res
    DTe4(DTe4 < T_res)=[]; % Removing trials with reaction time below T_res
    if(c==0.032)
        DTc4_032=DTc4;
        DTe4_032=DTe4;
        DTc4_032(DTc4_032==0)=[];
        DTe4_032(DTe4_032==0)=[];
    end
    if(c==0.128)
        DTc4_128=DTc4;
        DTe4_128=DTe4;
        DTc4_128(DTc4_128==0)=[];
        DTe4_128(DTe4_128==0)=[];
    end
    if(c==0.256)
        DTc4_256=DTc4;
        DTe4_256=DTe4;
        DTc4_256(DTc4_256==0)=[];
        DTe4_256(DTe4_256==0)=[];      
    end
    AveDTc_onlyDrift(j)=mean(DTc4); %Calculate the average reaction time

    correct_trial_number(4,j) = length(find(DTc4 > 0)); % storing the number of correct trials in each case
    error_trial_number(4,j) = length(find(DTe4 > 0)); % storing the number of error trials in each case

    if error_trial_number(4,j) < error_trial_number_threshold                  
        AveDTe_onlyDrift(j) = 0;      % Assigning zeto to mean reaction time if there is no error trials
    else
        AveDTe_onlyDrift(j)=mean(DTe4);
    end

%     AveDTe_onlyDrift(j)=dt.*mean(DTe4);
    ac4(j)= length(find(DTc4 > 0))/(length(find(DTc4 > 0)) + length(find(DTe4 > 0))); % Changed this to just include non-zero reaction times
    j=j+1;
end

% AveDTe_onlyDrift_temp = AveDTe_onlyDrift;
% AveDTe_onlyDrift = interpolate_vector(AveDTe_onlyDrift,interpolation_method);   % Interpolating zero values up until the last non-zero element 
AveDTe_onlyDrift(AveDTe_onlyDrift==0) = NaN;
%% Plot the 4x2 MotionStrength vs DT and Accuracy figure
disp(['Plotting figures at ' char(datetime)])

X1{1}=100*Mo_Strength;
X1{2}=log10(X1{1});

for accuracy_var_num = 1 : 4
    ac_temp = eval(['ac' num2str(accuracy_var_num)]);
    YMatrix{accuracy_var_num, 2} = ac_temp;
    ac_temp = [];
end
YMatrix{1, 1}=[AveDTc_initial; AveDTe_initial];
YMatrix{2, 1}=[AveDTc_both; AveDTe_both];
YMatrix{3, 1}=[AveDTc_onlyNoise; AveDTe_onlyNoise];
YMatrix{4, 1}=[AveDTc_onlyDrift; AveDTe_onlyDrift];
case_number = {'i', 'ii', 'iii', 'iv'};
common_xlabel = {'Motion strength (% coh)', 'Motion strength (Log10)'};
common_ylabel = {'Mean reaction time (ms)','% correct' };

num_of_subplots = numel(YMatrix);
subplot_names = cellstr(('a':char('a'+(num_of_subplots-1)))');

subplot_num = 1;

figure1 = figure('units','normalized','outerposition',[0 0 1 1]);
for row_num = 1 : num_of_subplots/2
    for column_num = 1 : length(X1)
% create subplot
        subplot1(subplot_num) = subplot(num_of_subplots/2,length(X1),subplot_num);
        hold(subplot1(subplot_num),'on');
        switch column_num
            case 1
                plot1{subplot_num} = plot(X1{column_num},YMatrix{row_num, column_num},'LineWidth',5);
                set(plot1{subplot_num}(1),'DisplayName','Correct','Color',[0 0 1]);
                set(plot1{subplot_num}(2),'DisplayName','Error','Color',[1 0 0]);
                title(['Case (' case_number{row_num} ')']);
                % Set the rest of the axes properties
                set(subplot1(subplot_num),'FontSize',16,'LineWidth',2,'TickDir','out');
                % create legend
                if subplot_num == 1
                    legend1 = legend(subplot1(subplot_num),'show');
                    set(legend1,...
                        'Position',[0.364390839087796 0.865014403540914 0.0896107984202805 0.0765746516829394]);
                end
                if row_num == num_of_subplots/2
                    xlabel(common_xlabel{column_num});
                end
            case 2
                plot(X1{column_num},YMatrix{row_num, column_num},'LineWidth',3.5,'Color',[0.717647058823529 0.274509803921569 1]);
                % Set the rest of the axes properties
                set(subplot1(subplot_num),'FontSize',16,'LineWidth',2,'TickDir','out','YTick',...
                    [0.5 0.75 1],'YTickLabel',{'50','75','100'});
                if row_num == num_of_subplots/2
                    xlabel(common_xlabel{column_num});
                end
        end

            pos(row_num, column_num) = {get(subplot1(subplot_num), 'position')};
            dim = cellfun(@(x) x.*[1 1 1 1], pos(row_num, column_num), 'uni',0);
            annotation(figure1, 'textbox', dim{1} + [-0.05 0.03 0 0], 'String', ['(' subplot_names{subplot_num} ')'], 'vert', 'top',...
                'FontWeight','bold',...
                'FontSize', 16,...
                'FitBoxToText','off',...
                 'EdgeColor','none');
            subplot_num = subplot_num + 1;
    end
end

mean_pos = mean(cell2mat(pos),1);
common_ylabel_pos = mat2cell(mean_pos, 1, 4 * ones(1, length(X1)));
ylabel_pos_bias = {[0.01 -0.05 0 0], [0.01 0 0 0]};
for ylabel_pos_num = 1 : length(common_ylabel_pos)
                ht(ylabel_pos_num) = annotation(figure1, 'textbox', common_ylabel_pos{ylabel_pos_num} + ylabel_pos_bias{ylabel_pos_num}, 'String', common_ylabel{ylabel_pos_num},...
                'FontSize', 18,...
                'FitBoxToText','off',...
                 'EdgeColor','none');
                set(ht(ylabel_pos_num),'Rotation',90)
                annotation_pos{ylabel_pos_num} = get(ht(ylabel_pos_num),'Position');

end

%% plot the RT distribution
% creating variables
responses = {'c','e'};
plotting_motion_strengths = {'032','128'};
figure2_column_titles = {'correct 3.2%', 'error 3.2%';
                        'correct 12.8%', 'error 12.8%'};
num_of_cases = length(name_of_cases);
color_of_plots = {[0 0 1], [1 0 0]};
var_num = 1;
for case_num = 1 : num_of_cases
    for motion_num = 1 : length(plotting_motion_strengths)
        for res_num = 1 : length(responses)
            variable_name{var_num} = ['DT' responses{res_num} num2str(case_num) '_' plotting_motion_strengths{motion_num}];
            var_num = var_num + 1;
        end
    end
end

num_of_subplots = length(variable_name);
subplot_names = cellstr(('a':char('a'+(num_of_subplots-1)))');
subplots_with_title = 1:num_of_cases
mo_strength032_correct_ylim = trials/20;
mo_strength032_error_ylim = mo_strength032_correct_ylim;
mo_strength128_correct_ylim = trials/10;
mo_strength128_error_ylim = trials/200;

mo_strength_ylim_temp = [mo_strength032_correct_ylim mo_strength032_error_ylim mo_strength128_correct_ylim mo_strength128_error_ylim];
mo_strength_ylim = repmat(mo_strength_ylim_temp, 1, length(plotting_motion_strengths) * length(responses));

% create figure
figure2 = figure('units','normalized','outerposition',[0 0 1 1]);

subplot_num = 1;
for case_num = 1 : num_of_cases
    for motion_num = 1 : length(plotting_motion_strengths)
        for res_num = 1 : length(responses)

            % create subplot
            subplots(subplot_num) = subplot(num_of_cases, length(plotting_motion_strengths) * length(responses),subplot_num);
            hold(subplots(subplot_num),'on');
            % create histogram
            clear eval
            histogram(eval(variable_name{subplot_num}),'DisplayName',subplot_names{subplot_num} ,'FaceColor',color_of_plots{res_num},'NumBins',50);
            % Set the rest of the axes properties
            set(subplots(subplot_num),'FontSize',14,'LineWidth',2,'TickDir','out','YTick',[0 mo_strength_ylim(subplot_num)/2 mo_strength_ylim(subplot_num)]);
            xlim(subplots(subplot_num),[0 2000]);
            ylim(subplots(subplot_num),[0 mo_strength_ylim(subplot_num)]);
            if ismember(subplot_num,subplots_with_title)
                title(figure2_column_titles{motion_num, res_num});
            end
            % create textbox 
            pos = {get(subplots(subplot_num), 'position')};
            dim = cellfun(@(x) x.*[1 1 1 1], pos, 'uni',0);
            annotation(figure2, 'textbox', dim{1} + [-0.04 0.04 0 0], 'String', ['(' subplot_names{subplot_num} ')'], 'vert', 'top',...
                'FontWeight','bold',...
                'FontSize', 16,...
                'FitBoxToText','off',...
                 'EdgeColor','none');
            subplot_num = subplot_num + 1;
        end
    end
end



% creating common xlabel and ylabel

han=axes(figure2,'visible','off'); 
% han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylh = ylabel(han,'Frequency',  'FontSize', 24, 'FontWeight','bold');
xlh = xlabel(han,'Mean reaction time (ms)', 'FontSize', 24, 'FontWeight','bold');
% title(han,'yourTitle');
xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.2);
ylh.Position(1) = ylh.Position(1) - abs(ylh.Position(1) * 3);

ElapsedTime = toc

disp(['Analysis finished at ' char(datetime)])


%% saving figures
figure_save_directory_name = 'Figures';
data_save_path = [scratch_main_dir filesep 'Hui paper'];
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

figure_name = 'Figure_1';
savefig(figure1,fullfile(figure_save_path, figure_name));
figure_name = 'Figure_2';
savefig(figure2,fullfile(figure_save_path, figure_name));

clear figure*

save(fullfile(data_save_path,'dataset.mat'));
% save('dataset.mat');
