% sensitivity_tornado_plots_atrial.m
% Atrial-specific sensitivity tornado plots sorted by magnitude

clc; clear; close all;

%% Load Atrial Biomarkers and Parameters
load('results/full_pom_analysis_atrial.mat');

atrial_biomarkers = [APD90_atrial, Ca_amp_atrial, Ca_decay50_atrial];
atrial_bio_names = {'Atrial APD90', 'Atrial Ca Amplitude', 'Atrial Ca Decay50'};

num_params_atrial = length(param_names_atrial);
num_biomarkers_atrial = length(atrial_bio_names);

atrial_sensitivities = zeros(num_params_atrial, num_biomarkers_atrial);
atrial_p_values = zeros(num_params_atrial, num_biomarkers_atrial);

%% Calculate Spearman correlations for atrial biomarkers
for j = 1:num_params_atrial
    for k = 1:num_biomarkers_atrial
        [rho, pval] = corr(param_sets_atrial(:,j), atrial_biomarkers(:,k), 'Type','Spearman');
        atrial_sensitivities(j,k) = rho;
        atrial_p_values(j,k) = pval;
    end
end

%% Generate Explicitly Sorted Atrial Tornado Plots
for k = 1:num_biomarkers_atrial
    figure;

    [~, idx_sort] = sort(abs(atrial_sensitivities(:,k)), 'ascend');
    sorted_params = param_names_atrial(idx_sort);
    barh(atrial_sensitivities(idx_sort,k), 'FaceColor', 'cyan');

    yticks(1:num_params_atrial);
    yticklabels(sorted_params);
    ytickangle(45);

    xlabel('Spearman Correlation');
    title(['Tornado Plot: Sensitivity of ', atrial_bio_names{k}]);
    grid on;

    saveas(gcf, ['results/', atrial_bio_names{k}, '_Tornado_atrial.png']);
end