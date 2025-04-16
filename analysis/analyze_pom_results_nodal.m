% analyze_pom_results_nodal.m
% Integrated Analysis of Nodal Population-of-Models (PoM), SERCA/RyR, and Sensitivity Analysis

clc; clear; close all;

%% Load individual nodal PoM results
nodal_model_files = dir('results/nodal_models/model_*.mat');
numNodalModels = length(nodal_model_files);

%% Load time from the first nodal model
first_nodal_model = load(fullfile(nodal_model_files(1).folder, nodal_model_files(1).name));
t_fixed_nodal = first_nodal_model.t_fixed;
analysis_window_nodal = (t_fixed_nodal >= 780 & t_fixed_nodal <= 800);
t_analysis_nodal = t_fixed_nodal(analysis_window_nodal);

%% Preallocate biomarker and parameter storage
APD90_nodal = NaN(numNodalModels,1); 
DMP_nodal = NaN(numNodalModels,1);
peakVm_nodal = NaN(numNodalModels,1); 
spontInterval_nodal = NaN(numNodalModels,1);
Ca_peak_nodal = NaN(numNodalModels,1); 
Ca_base_nodal = NaN(numNodalModels,1);
Ca_amp_nodal = NaN(numNodalModels,1); 
Ca_decay50_nodal = NaN(numNodalModels,1);

param_names_nodal = fieldnames(first_nodal_model.params);
param_sets_nodal = NaN(numNodalModels, length(param_names_nodal));

Iup_peaks_nodal = NaN(numNodalModels,1);
Irel_peaks_nodal = NaN(numNodalModels,1);

%% Biomarkers and parameter extraction loop
parfor i = 1:numNodalModels
    try
        data = load(fullfile(nodal_model_files(i).folder, nodal_model_files(i).name));

        % Safely extract and assign parameters
        param_values = zeros(1, length(param_names_nodal));
        for p = 1:length(param_names_nodal)
            value = data.params.(param_names_nodal{p});
            if isnumeric(value) && isscalar(value)
                param_values(p) = value;
            else
                error('Non-numeric or non-scalar parameter in model %d: %s', i, param_names_nodal{p});
            end
        end
        param_sets_nodal(i,:) = param_values;

        % Analyze traces
        Vm_analysis = data.Vm(analysis_window_nodal)*1e3;
        Cai_analysis = data.Cai(analysis_window_nodal)*1e6;

        [pks, locs] = findpeaks(Vm_analysis);
        if isempty(pks)
            DMP_nodal(i) = mean(Vm_analysis);
        else
            spontInterval_nodal(i) = mean(diff(t_analysis_nodal(locs)));
            peakVm_nodal(i) = mean(pks);
            DMP_nodal(i) = mean(Vm_analysis(Vm_analysis < 0));

            APD90_list = [];
            for j = 1:length(locs)-1
                segment = Vm_analysis(locs(j):locs(j+1));
                repol90 = pks(j) - 0.9*(pks(j)-min(segment));
                crossIdx = find(segment <= repol90, 1);
                if ~isempty(crossIdx)
                    APD90_list(end+1) = ...
                        (t_analysis_nodal(locs(j)+crossIdx-1)-t_analysis_nodal(locs(j)))*1000;
                end
            end
            APD90_nodal(i) = mean(APD90_list, 'omitnan');
        end

        % Calcium biomarkers
        Ca_peak_nodal(i) = max(Cai_analysis);
        Ca_base_nodal(i) = min(Cai_analysis);
        Ca_amp_nodal(i) = Ca_peak_nodal(i) - Ca_base_nodal(i);

        peakIdx = find(Cai_analysis == Ca_peak_nodal(i), 1);
        decay50Idx = find(Cai_analysis(peakIdx:end) <= ...
                        (Ca_peak_nodal(i) - 0.5*Ca_amp_nodal(i)), 1) + peakIdx - 1;

        if ~isempty(decay50Idx)
            Ca_decay50_nodal(i) = t_analysis_nodal(decay50Idx) - t_analysis_nodal(peakIdx);
        end

        Iup_peaks_nodal(i) = max(data.Iup(analysis_window_nodal));
        Irel_peaks_nodal(i) = max(data.Irel(analysis_window_nodal));

    catch ME
        warning('Skipping model %d: %s', i, ME.message);
        continue;
    end
end

%% Plot distributions of biomarkers
figure;
subplot(2,2,1); histogram(APD90_nodal); xlabel('APD90 (ms)'); title('Nodal APD90 Distribution'); grid on;
subplot(2,2,2); histogram(DMP_nodal); xlabel('DMP (mV)'); title('Nodal DMP Distribution'); grid on;
subplot(2,2,3); histogram(peakVm_nodal); xlabel('Peak Vm (mV)'); title('Nodal Peak Vm Distribution'); grid on;
subplot(2,2,4); histogram(spontInterval_nodal); xlabel('Interval (s)'); title('Nodal Spontaneous Interval'); grid on;
saveas(gcf,'results/Biomarker_distributions_nodal.png');

%% SERCA vs Ca Amp
figure;
scatter(Iup_peaks_nodal, Ca_amp_nodal,'filled');
xlabel('I_{up} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('Nodal SERCA vs. Ca Amplitude'); grid on;
saveas(gcf,'results/SERCA_Ca_amp_nodal.png');

%% RyR vs Ca Amp
figure;
scatter(Irel_peaks_nodal, Ca_amp_nodal,'filled');
xlabel('I_{rel} Peak (A/F)'); ylabel('Ca Amplitude (µM)');
title('Nodal RyR vs. Ca Amplitude'); grid on;
saveas(gcf,'results/RyR_Ca_amp_nodal.png');

%% Sensitivity analysis (Spearman correlation)
nodal_biomarkers = [APD90_nodal, Ca_amp_nodal, Ca_decay50_nodal];
nodal_bio_names = {'APD90','Ca Amplitude','Ca Decay50'};

sensitivities_nodal = zeros(length(param_names_nodal), length(nodal_bio_names));
p_values_nodal = zeros(length(param_names_nodal), length(nodal_bio_names));

for j = 1:length(param_names_nodal)
    for k = 1:length(nodal_bio_names)
        [rho, pval] = corr(param_sets_nodal(:,j), nodal_biomarkers(:,k), 'Type','Spearman', 'Rows','complete');
        sensitivities_nodal(j,k) = rho;
        p_values_nodal(j,k) = pval;
    end
end

%% Save results
save('results/full_pom_analysis_nodal.mat',...
     'APD90_nodal','DMP_nodal','peakVm_nodal','spontInterval_nodal',...
     'Ca_peak_nodal','Ca_base_nodal','Ca_amp_nodal','Ca_decay50_nodal',...
     'Iup_peaks_nodal','Irel_peaks_nodal',...
     'sensitivities_nodal','p_values_nodal',...
     'param_sets_nodal','param_names_nodal');
