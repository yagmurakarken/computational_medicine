% analyze_pom_results_atrial.m
% Integrated Analysis of Atrial Population-of-Models (PoM), SERCA/RyR, and Sensitivity Analysis

clc; clear; close all;

%% Load individual atrial PoM results
atrial_model_files = dir('results/atrial_models/model_*.mat');
numAtrialModels = length(atrial_model_files);

%% Load time from the first atrial model
first_atrial_model = load(fullfile(atrial_model_files(1).folder, atrial_model_files(1).name));
t_fixed_atrial = first_atrial_model.t_fixed;
analysis_window_atrial = (t_fixed_atrial >= 780 & t_fixed_atrial <= 800);
t_analysis_atrial = t_fixed_atrial(analysis_window_atrial);

%% Preallocate atrial biomarkers and parameters
APD90_atrial = zeros(numAtrialModels,1); 
DMP_atrial = zeros(numAtrialModels,1);
peakVm_atrial = zeros(numAtrialModels,1); 
spontInterval_atrial = zeros(numAtrialModels,1);
Ca_peak_atrial = zeros(numAtrialModels,1); 
Ca_base_atrial = zeros(numAtrialModels,1);
Ca_amp_atrial = zeros(numAtrialModels,1); 
Ca_decay50_atrial = zeros(numAtrialModels,1);

param_sets_atrial = zeros(numAtrialModels, length(fieldnames(first_atrial_model.params)));
param_names_atrial = fieldnames(first_atrial_model.params);

Iup_peaks_atrial = zeros(numAtrialModels,1);
Irel_peaks_atrial = zeros(numAtrialModels,1);

%% Biomarkers calculation for atrial phenotype
parfor i = 1:numAtrialModels
    data = load(fullfile(atrial_model_files(i).folder, atrial_model_files(i).name));
    param_values = cellfun(@(f) data.params.(f), param_names_atrial);
    param_sets_atrial(i,:) = param_values;

    Vm_analysis = data.Vm(analysis_window_atrial)*1e3;
    Cai_analysis = data.Cai(analysis_window_atrial)*1e6;

    [pks, locs] = findpeaks(Vm_analysis);

    if isempty(pks)
        spontInterval_atrial(i) = NaN;
        peakVm_atrial(i) = NaN;
        DMP_atrial(i) = mean(Vm_analysis);
        APD90_atrial(i) = NaN;
    else
        spontInterval_atrial(i) = mean(diff(t_analysis_atrial(locs)));
        peakVm_atrial(i) = mean(pks);
        DMP_atrial(i) = mean(Vm_analysis(Vm_analysis < 0));

        APD90_list = [];
        for j = 1:length(locs)-1
            segment = Vm_analysis(locs(j):locs(j+1));
            repol90 = pks(j) - 0.9*(pks(j)-min(segment));
            crossIdx = find(segment <= repol90, 1);
            if ~isempty(crossIdx)
                APD90_list(end+1) = (t_analysis_atrial(locs(j)+crossIdx-1)-t_analysis_atrial(locs(j)))*1000;
            end
        end
        APD90_atrial(i) = mean(APD90_list, 'omitnan');
    end

    Ca_peak_atrial(i) = max(Cai_analysis);
    Ca_base_atrial(i) = min(Cai_analysis);
    Ca_amp_atrial(i) = Ca_peak_atrial(i)-Ca_base_atrial(i);
    peakIdx = find(Cai_analysis == Ca_peak_atrial(i),1);
    decay50Idx = find(Cai_analysis(peakIdx:end) <= (Ca_peak_atrial(i)-0.5*Ca_amp_atrial(i)),1) + peakIdx - 1;
    Ca_decay50_atrial(i) = t_analysis_atrial(decay50Idx) - t_analysis_atrial(peakIdx);

    Iup_peaks_atrial(i) = max(data.Iup(analysis_window_atrial));
    Irel_peaks_atrial(i) = max(data.Irel(analysis_window_atrial));
end

%% Plot Atrial Biomarkers explicitly labeled
figure;
subplot(2,2,1); histogram(APD90_atrial); xlabel('APD90 (ms)'); title('Atrial APD90 Distribution'); grid on;
subplot(2,2,2); histogram(DMP_atrial); xlabel('DMP (mV)'); title('Atrial DMP Distribution'); grid on;
subplot(2,2,3); histogram(peakVm_atrial); xlabel('Peak Vm (mV)'); title('Atrial Peak Vm Distribution'); grid on;
subplot(2,2,4); histogram(spontInterval_atrial); xlabel('Interval (s)'); title('Atrial Spontaneous Interval'); grid on;
saveas(gcf,'results/Biomarker_distributions_atrial.png');

%% Atrial SERCA/RyR Analysis
figure;
scatter(Iup_peaks_atrial, Ca_amp_atrial,'filled');
xlabel('I_{up} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('Atrial SERCA vs. Ca Amplitude'); grid on;
saveas(gcf,'results/SERCA_Ca_amp_atrial.png');

figure;
scatter(Irel_peaks_atrial, Ca_amp_atrial,'filled');
xlabel('I_{rel} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('Atrial RyR vs. Ca Amplitude'); grid on;
saveas(gcf,'results/RyR_Ca_amp_atrial.png');

%% Atrial Sensitivity Analysis
atrial_biomarkers = [APD90_atrial, Ca_amp_atrial, Ca_decay50_atrial];
atrial_bio_names = {'APD90','Ca Amplitude','Ca Decay50'};

sensitivities_atrial = zeros(length(param_names_atrial), length(atrial_bio_names));
p_values_atrial = zeros(length(param_names_atrial), length(atrial_bio_names));

for j = 1:length(param_names_atrial)
    for k = 1:length(atrial_bio_names)
        [rho, pval] = corr(param_sets_atrial(:,j), atrial_biomarkers(:,k), 'Type','Spearman');
        sensitivities_atrial(j,k) = rho;
        p_values_atrial(j,k) = pval;
    end
end

save('results/full_pom_analysis_atrial.mat','APD90_atrial','DMP_atrial','peakVm_atrial',...
    'spontInterval_atrial','Ca_peak_atrial','Ca_base_atrial','Ca_amp_atrial','Ca_decay50_atrial',...
    'sensitivities_atrial','p_values_atrial','param_sets_atrial','param_names_atrial');