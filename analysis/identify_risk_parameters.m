clc; clear; close all;

%% Load filtered biomarker data
load('results/biomarkers_ventricular_filtered.mat'); % contains: valid_idx, param_sets_valid, param_names, APD90_valid, Ca_amp_valid, Ca_decay50_valid

%% Combine Biomarkers
biomarkers = [APD90_valid, Ca_amp_valid, Ca_decay50_valid];

%% Define thresholds explicitly (percentiles)
high_thresholds = prctile(biomarkers, [75 75 75]);
low_thresholds = prctile(biomarkers, [25 50 25]);

%% Identify High-risk models (at least 2 criteria)
high_risk_idx = find(...
    (biomarkers(:,1) >= high_thresholds(1)) + ...
    (biomarkers(:,2) >= high_thresholds(2)) + ...
    (biomarkers(:,3) >= high_thresholds(3)) >= 2);

%% Identify Low-risk models (at least 2 criteria)
low_risk_idx = find(...
    (biomarkers(:,1) <= low_thresholds(1)) + ...
    (abs(biomarkers(:,2)-low_thresholds(2)) <= 0.05*low_thresholds(2)) + ...
    (biomarkers(:,3) <= low_thresholds(3)) >= 2);

%% Fallback scoring if no matches
if isempty(high_risk_idx)
    fprintf('No high-risk set matches 2 criteria, using scoring...\n');
    risk_score_high = zscore(biomarkers(:,1)) + zscore(biomarkers(:,2)) + zscore(biomarkers(:,3));
    [~, high_risk_idx] = max(risk_score_high);
else
    high_risk_idx = high_risk_idx(1);
end

if isempty(low_risk_idx)
    fprintf('No low-risk set matches 2 criteria, using scoring...\n');
    risk_score_low = -zscore(biomarkers(:,1)) - abs(zscore(biomarkers(:,2))) - zscore(biomarkers(:,3));
    [~, low_risk_idx] = max(risk_score_low);
else
    low_risk_idx = low_risk_idx(1);
end

%% Display High-risk info
fprintf('\n--- Selected High-risk parameter set ---\n');
fprintf('Model ID: %d\n', valid_idx(high_risk_idx));
fprintf('APD90: %.2f ms, Ca Amplitude: %.2f µM, Ca Decay50: %.2f ms\n', ...
    biomarkers(high_risk_idx,1), biomarkers(high_risk_idx,2), biomarkers(high_risk_idx,3));
for p = 1:length(param_names)
    fprintf('\t%s: %.6f\n', param_names{p}, param_sets_valid(high_risk_idx,p));
end

%% Display Low-risk info
fprintf('\n--- Selected Low-risk parameter set ---\n');
fprintf('Model ID: %d\n', valid_idx(low_risk_idx));
fprintf('APD90: %.2f ms, Ca Amplitude: %.2f µM, Ca Decay50: %.2f ms\n', ...
    biomarkers(low_risk_idx,1), biomarkers(low_risk_idx,2), biomarkers(low_risk_idx,3));
for p = 1:length(param_names)
    fprintf('\t%s: %.6f\n', param_names{p}, param_sets_valid(low_risk_idx,p));
end

%% Save IDs and info
save('results/selected_parameters_ventricular.mat', 'high_risk_idx', 'low_risk_idx', ...
    'param_sets_valid', 'param_names', 'biomarkers', 'valid_idx');

%% Load and plot high-risk model
hr_model = load(sprintf('results/ventricular_models/model_098.mat', valid_idx(high_risk_idx)));
t = hr_model.t_fixed;
Vm = hr_model.Vm * 1e3;
Cai = hr_model.Cai * 1e6;

figure;
yyaxis left
plot(t, Vm, 'LineWidth', 1.2); ylabel('Vm (mV)');
yyaxis right
plot(t, Cai, 'LineWidth', 1.2); ylabel('Ca_i (\muM)');
xlabel('Time (s)');
title('High-risk AP/Ca - Ventricular');
xlim([790 800]); grid on;
saveas(gcf, 'results/High_risk_AP_ventricular.png');

%% Load and plot low-risk model
lr_model = load(sprintf('results/ventricular_models/model_008.mat', valid_idx(low_risk_idx)));
t = lr_model.t_fixed;
Vm = lr_model.Vm * 1e3;
Cai = lr_model.Cai * 1e6;

figure;
yyaxis left
plot(t, Vm, 'LineWidth', 1.2); ylabel('Vm (mV)');
yyaxis right
plot(t, Cai, 'LineWidth', 1.2); ylabel('Ca_i (\muM)');
xlabel('Time (s)');
title('Low-risk AP/Ca - Ventricular');
xlim([790 800]); grid on;
saveas(gcf, 'results/Low_risk_AP_ventricular.png');

fprintf('\n Figure 3 components (ventricular) successfully generated and saved.\n');
