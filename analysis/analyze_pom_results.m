% analyze_pom_results.m
% Robust Analysis of Population-of-Models (PoM) with Handling of Non-Physiological Cases

clc; clear; close all;

%% Load individual PoM results
model_files = dir('results/individual_models/model_*.mat');
numModels = length(model_files);

%% Load simulation time once from the first model
first_model = load(fullfile(model_files(1).folder, model_files(1).name));
t_fixed = first_model.t_fixed;
t = t_fixed;
analysis_window = (t >= 780 & t <= 800);
t_analysis = t(analysis_window);

%% Preallocate Biomarkers and parameters
APD90 = zeros(numModels,1); 
DMP = zeros(numModels,1);
peakVm = zeros(numModels,1); 
spontInterval = zeros(numModels,1);
Ca_peak = zeros(numModels,1); 
Ca_base = zeros(numModels,1);
Ca_amp = zeros(numModels,1); 
Ca_decay50 = zeros(numModels,1);
param_sets = zeros(numModels, length(fieldnames(first_model.params)));
param_names = fieldnames(first_model.params);

Iup_peaks = zeros(numModels,1);
Irel_peaks = zeros(numModels,1);

%% Biomarkers Calculation with checks
for i = 1:numModels
    data = load(fullfile(model_files(i).folder, model_files(i).name));
    param_values = cellfun(@(f) data.params.(f), param_names);
    param_sets(i,:) = param_values;

    Vm_analysis = data.Vm(analysis_window)*1e3; % in mV
    Cai_analysis = data.Cai(analysis_window)*1e6; % in µM

    [pks, locs] = findpeaks(Vm_analysis);

    % check for sufficient number of peaks
    if length(pks) < 2
        spontInterval(i) = NaN;
        peakVm(i) = NaN;
        DMP(i) = mean(Vm_analysis);
        APD90(i) = NaN;
    else
        intervals = diff(t_analysis(locs));
        spontInterval(i) = mean(intervals);
        peakVm(i) = mean(pks);
        DMP(i) = mean(Vm_analysis(Vm_analysis < 0));

        % APD90 calculation
        APD90_list = [];
        for j = 1:length(locs)-1
            segment = Vm_analysis(locs(j):locs(j+1));
            repol90 = pks(j) - 0.9*(pks(j)-min(segment));
            crossIdx = find(segment <= repol90, 1, 'first');
            if ~isempty(crossIdx)
                APD90_temp = (t_analysis(locs(j)+crossIdx-1)-t_analysis(locs(j)))*1000;
                APD90_list(end+1) = APD90_temp; %#ok<AGROW>
            end
        end
        APD90(i) = mean(APD90_list,'omitnan');
    end

    % Calcium biomarkers checked
    [Ca_peak_val, peakIdx] = max(Cai_analysis);
    [Ca_base_val, ~] = min(Cai_analysis);

    Ca_peak(i) = Ca_peak_val;
    Ca_base(i) = Ca_base_val;
    Ca_amp(i) = Ca_peak_val - Ca_base_val;

    decay_threshold = Ca_peak_val - 0.5*Ca_amp(i);
    decay50_relative_idx = find(Cai_analysis(peakIdx:end) <= decay_threshold, 1, 'first');

    if ~isempty(decay50_relative_idx)
        decay50Idx = peakIdx + decay50_relative_idx - 1;
        Ca_decay50(i) = t_analysis(decay50Idx) - t_analysis(peakIdx);
    else
        Ca_decay50(i) = NaN;
    end

    % Current peaks handled
    Iup_peaks(i) = max(data.Iup(analysis_window));
    Irel_peaks(i) = max(data.Irel(analysis_window));
end

%% filtering of non-physiological models
valid_idx = find(peakVm > 0 & APD90 < 1000 & spontInterval < 5 & ~isnan(APD90));

APD90_valid = APD90(valid_idx);
DMP_valid = DMP(valid_idx);
peakVm_valid = peakVm(valid_idx);
spontInterval_valid = spontInterval(valid_idx);
Ca_amp_valid = Ca_amp(valid_idx);
Ca_decay50_valid = Ca_decay50(valid_idx);
param_sets_valid = param_sets(valid_idx,:);
Iup_peaks_valid = Iup_peaks(valid_idx); 
Irel_peaks_valid = Irel_peaks(valid_idx); 

%% Save filtered results including currents
save('results/biomarkers_ventricular_filtered.mat', ...
    'APD90_valid', 'DMP_valid', 'peakVm_valid', 'spontInterval_valid', ...
    'Ca_amp_valid', 'Ca_decay50_valid', 'param_sets_valid', 'valid_idx', ...
    'Iup_peaks_valid', 'Irel_peaks_valid', 'param_names');

fprintf('biomarker analysis and filtering completed successfully.\n');
