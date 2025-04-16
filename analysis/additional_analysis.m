% sensitivity_tornado_plots.m
% Adjusted sensitivity tornado plots sorted by magnitude, with SERCA/RyR analysis

clc; clear; close all;

%% Load filtered biomarker results
load('results/biomarkers_ventricular_filtered.mat');

%% Plot Biomarker Distributions
figure;
subplot(2,2,1); histogram(APD90_valid); xlabel('APD90 (ms)'); title('APD90 Distribution'); grid on;
subplot(2,2,2); histogram(DMP_valid); xlabel('DMP (mV)'); title('DMP Distribution'); grid on;
subplot(2,2,3); histogram(peakVm_valid); xlabel('Peak Vm (mV)'); title('Peak Vm Distribution'); grid on;
subplot(2,2,4); histogram(spontInterval_valid); xlabel('Spont Interval (s)'); title('Spontaneous Interval Distribution'); grid on;
saveas(gcf,'results/Biomarker_distributions_filtered.png');

%% SERCA/RyR Analysis and Validation
figure;
scatter(Iup_peaks_valid, Ca_amp_valid,'filled');
xlabel('I_{up} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('SERCA vs. Ca Transient Amplitude'); grid on;
saveas(gcf,'results/SERCA_Ca_amp_filtered.png');

figure;
scatter(Irel_peaks_valid, Ca_amp_valid,'filled');
xlabel('I_{rel} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('RyR vs. Ca Transient Amplitude'); grid on;
saveas(gcf,'results/RyR_Ca_amp_filtered.png');

%% Sensitivity Analysis (Refined)
biomarkers = [APD90_valid, Ca_amp_valid, Ca_decay50_valid];
bio_names = {'APD90','Ca Amplitude','Ca Decay50'};

num_params = size(param_sets_valid,2);
num_biomarkers = length(bio_names);

sensitivities = zeros(num_params, num_biomarkers);
p_values = zeros(num_params, num_biomarkers);

%% Calculate Spearman correlations
for j = 1:num_params
    for k = 1:num_biomarkers
        [rho, pval] = corr(param_sets_valid(:,j), biomarkers(:,k), 'Type','Spearman');
        sensitivities(j,k) = rho;
        p_values(j,k) = pval;
    end
end

%% Generate Sorted Tornado plots
for k = 1:num_biomarkers
    figure;

    [~, idx_sort] = sort(abs(sensitivities(:,k)), 'ascend');
    sorted_params = param_names(idx_sort);
    barh(sensitivities(idx_sort,k), 'FaceColor', 'cyan');

    yticks(1:num_params);
    yticklabels(sorted_params);
    ytickangle(45);

    xlabel('Spearman Correlation');
    title(['Tornado Plot: Sensitivity of ', bio_names{k}]);
    grid on;

    saveas(gcf, ['results/', bio_names{k}, '_Tornado_filtered.png']);
end

%% Sensitivity Scatter Plots
figure;
for j = 1:num_params
    for k = 1:num_biomarkers
        subplot(num_params, num_biomarkers, (j-1)*num_biomarkers+k);
        scatter(param_sets_valid(:,j), biomarkers(:,k),'filled');
        xlabel(param_names{j}); ylabel(bio_names{k});
        title(sprintf('\\rho=%.2f, p=%.3f',sensitivities(j,k),p_values(j,k)));
        grid on;
    end
end
sgtitle('Refined Sensitivity Analysis (Spearman)');
saveas(gcf,'results/Sensitivity_Analysis_filtered.png');
