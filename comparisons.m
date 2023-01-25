clearvars; close all;
current_path = pwd;
load([current_path filesep 'Hui paper' filesep 'dataset.mat']) % change path to load the dataset
coherence_level = {'256'; '032'};
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
%                 trials_tmp(:,2) = case_number * ones(case_length, 1);
%                 trials_tmp(:,3) = coherence * ones(case_length, 1);
%                 trials_tmp(:,4) = correctness * ones(case_length, 1);
                trials_tmp(:,2) = (eval(['DT' correct_error{correctness} num2str(case_number) '_' coherence_level{coherence}]))';
            end
            trials = [trials;trials_tmp];

            trials_tmp = [];
            group_number = group_number + 1;
        end
    end
end

% [p, observeddifference, effectsize] = permutationTest(DTc4_512, DTc4_032, 20000, 'plotresult', 1, 'showprogress', 250);
% group_name = cellstr(num2str(trials(:,1)));
final_group_name = cellstr(group_name);
[p,tbl,stats] = kruskalwallis(trials(:,2),trials(:,1), 'on');
figure
c = multcompare(stats, "CriticalValueType","tukey-kramer");
non_significant_groups = find(c(:,6) > 0.001);
non_significant_groups_pair(:,1) = final_group_name(c(non_significant_groups,1));
non_significant_groups_pair(:,2) = final_group_name(c(non_significant_groups,2))
