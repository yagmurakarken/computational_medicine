% sensitivity_tornado_plots_nodal.m
% Nodal-specific sensitivity tornado plots sorted by magnitude

clc; clear; close all;

%% Load Nodal Biomarkers and Parameters
load('results/full_pom_analysis_nodal.mat');

nodal_biomarkers = [APD90_nodal, Ca_amp_nodal, Ca_decay50_nodal];
nodal_bio_names = {'Nodal APD90', 'Nodal Ca Amplitude', 'Nodal Ca Decay50'};

num_params_nodal = length(param_names_nodal);
num_biomarkers_nodal = length(nodal_bio_names);

nodal_sensitivities = zeros(num_params_nodal, num_biomarkers_nodal);
nodal_p_values = zeros(num_params_nodal, num_biomarkers_nodal);

%% Calculate Spearman correlations
for j = 1:num_params_nodal
    for k = 1:num_biomarkers_nodal
        [rho, pval] = corr(param_sets_nodal(:,j), nodal_biomarkers(:,k), ...
                          'Type', 'Spearman', 'Rows', 'complete');
        nodal_sensitivities(j,k) = rho;
        nodal_p_values(j,k) = pval;
    end
end

%% Generate Explicitly Sorted Tornado Plots
for k = 1:num_biomarkers_nodal
    figure;

    [~, idx_sort] = sort(abs(nodal_sensitivities(:,k)), 'ascend');
    sorted_params = param_names_nodal(idx_sort);
    barh(nodal_sensitivities(idx_sort,k), 'FaceColor', 'cyan');

    yticks(1:num_params_nodal);
    yticklabels(sorted_params);
    ytickangle(45);

    xlabel('Spearman Correlation');
    title(['Tornado Plot: Sensitivity of ', nodal_bio_names{k}]);
    grid on;

    saveas(gcf, ['results/', nodal_bio_names{k}, '_Tornado_nodal.png']);
end
