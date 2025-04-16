% identify_risk_parameters_nodal.m
% Identify high-risk and low-risk nodal models based on biomarker thresholds

clc; clear; close all;

%% Load existing PoM analysis
load('results/full_pom_analysis_nodal.mat'); % contains APD90_nodal, Ca_amp_nodal, Ca_decay50_nodal, param_sets_nodal, param_names_nodal

%% Combine Biomarkers explicitly
biomarkers = [APD90_nodal, Ca_amp_nodal, Ca_decay50_nodal];

%% Define thresholds explicitly for each biomarker
high_thresholds = [prctile(APD90_nodal, 75), prctile(Ca_amp_nodal, 75), prctile(Ca_decay50_nodal, 75)];
low_thresholds  = [prctile(APD90_nodal, 25), prctile(Ca_amp_nodal, 50), prctile(Ca_decay50_nodal, 25)];

%% Identify High-risk models (at least 2 criteria)
high_risk_criteria = (biomarkers(:,1) >= high_thresholds(1)) + ...
                     (biomarkers(:,2) >= high_thresholds(2)) + ...
                     (biomarkers(:,3) >= high_thresholds(3));
high_risk_idx = find(high_risk_criteria >= 2);

%% Identify Low-risk models (at least 2 criteria)
low_risk_criteria = (biomarkers(:,1) <= low_thresholds(1)) + ...
                    (abs(biomarkers(:,2) - low_thresholds(2)) <= 0.05 * low_thresholds(2)) + ...
                    (biomarkers(:,3) <= low_thresholds(3));
low_risk_idx = find(low_risk_criteria >= 2);

%% Handle no matches by combined scoring (High-risk)
if isempty(high_risk_idx)
    fprintf('No high-risk set matches 2 criteria, using combined scoring...\n');
    risk_score_high = zscore(biomarkers(:,1)) + zscore(biomarkers(:,2)) + zscore(biomarkers(:,3));
    [~, high_risk_idx] = max(risk_score_high);
else
    high_risk_idx = high_risk_idx(1);
end

%% Handle no matches by combined scoring (Low-risk)
if isempty(low_risk_idx)
    fprintf('No low-risk set matches 2 criteria, using combined scoring...\n');
    risk_score_low = -zscore(biomarkers(:,1)) - abs(zscore(biomarkers(:,2))) - zscore(biomarkers(:,3));
    [~, low_risk_idx] = max(risk_score_low);
else
    low_risk_idx = low_risk_idx(1);
end

%% Display selected High-risk results
fprintf('\n--- Selected High-risk Nodal Parameter Set ---\n');
fprintf('Model ID: %d\n', high_risk_idx);
fprintf('APD90: %.2f ms, Ca Amplitude: %.2f µM, Ca Decay50: %.2f ms\n', ...
    APD90_nodal(high_risk_idx), Ca_amp_nodal(high_risk_idx), Ca_decay50_nodal(high_risk_idx));
fprintf('Parameter values:\n');
for p = 1:length(param_names_nodal)
    fprintf('\t%s: %.6f\n', param_names_nodal{p}, param_sets_nodal(high_risk_idx,p));
end

%% Display selected Low-risk results
fprintf('\n--- Selected Low-risk Nodal Parameter Set ---\n');
fprintf('Model ID: %d\n', low_risk_idx);
fprintf('APD90: %.2f ms, Ca Amplitude: %.2f µM, Ca Decay50: %.2f ms\n', ...
    APD90_nodal(low_risk_idx), Ca_amp_nodal(low_risk_idx), Ca_decay50_nodal(low_risk_idx));
fprintf('Parameter values:\n');
for p = 1:length(param_names_nodal)
    fprintf('\t%s: %.6f\n', param_names_nodal{p}, param_sets_nodal(low_risk_idx,p));
end

%% Save selected parameters explicitly
save('results/selected_parameters_nodal.mat', 'high_risk_idx', 'low_risk_idx', ...
     'param_sets_nodal', 'param_names_nodal', 'biomarkers');
