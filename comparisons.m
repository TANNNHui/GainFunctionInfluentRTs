close all; clear
filePath = matlab.desktop.editor.getActiveFilename;
for i = length(filePath):-1:1
    if filePath(i) == filesep
        slash(i,1) = 1;
    end
end
filePath = filePath(1:length(slash) ); cd (filePath);   % change to code's path to load the dataset
load([filePath filesep 'dataset and figures' filesep 'dataset.mat']) 
%% Assigning variables
coherence_level = {'128'; '032'};
correct_error = {'c'; 'e'};
case_numbers = 4;
group_number = 1;
trials = [];
group_name = [];
for case_number = 1 : case_numbers
    for coherence = 1 : length(coherence_level)
        for correctness = 1 : length(correct_error)
            case_length = length(eval(['DT' correct_error{correctness} num2str(case_number) '_' coherence_level{coherence}]));
            if case_length ~= 0
                trials_tmp(:,1) = group_number * ones(case_length, 1);
                group_name_tmp = (['DT' correct_error{correctness} num2str(case_number) '_' coherence_level{coherence}]);
                group_name = [group_name; group_name_tmp];
                trials_tmp(:,2) = (eval(['DT' correct_error{correctness} num2str(case_number) '_' coherence_level{coherence}]))';
            end
            trials = [trials;trials_tmp];

            trials_tmp = [];
            group_number = group_number + 1;
        end
    end
end

%% Performing Kruskal-Wallis test
final_group_name = cellstr(group_name);
[p,tbl,stats] = kruskalwallis(trials(:,2),trials(:,1), 'on');
figure
c = multcompare(stats, "CriticalValueType","tukey-kramer");     % pairwise comparisons
tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
non_significant_groups = find(c(:,6) > 0.001);
non_significant_groups_pair(:,1) = final_group_name(c(non_significant_groups,1));
non_significant_groups_pair(:,2) = final_group_name(c(non_significant_groups,2));
disp('Non-significant pairs:')
non_significant_groups_pair

