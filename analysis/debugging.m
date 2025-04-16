% finalize_selected_parameters.m
% saves corrected high-risk and low-risk parameter sets from PoM analysis

clc; clear; close all;

%% Load existing results
load('results/full_pom_analysis.mat');

%% defined Model IDs (from previous analysis)
high_risk_id = 11;
low_risk_id = 1;

%% Retrieve parameters and biomarkers
selected_params = struct();

selected_params.high_risk.model_id = high_risk_id;
selected_params.high_risk.APD90 = APD90(high_risk_id);
selected_params.high_risk.Ca_amp = Ca_amp(high_risk_id);
selected_params.high_risk.Ca_decay50 = Ca_decay50(high_risk_id);
selected_params.high_risk.g_CaL = param_sets(high_risk_id,1);
selected_params.high_risk.kNaCa = param_sets(high_risk_id,2);
selected_params.high_risk.VmaxUp = param_sets(high_risk_id,3);
selected_params.high_risk.g_irel_max = param_sets(high_risk_id,4);

selected_params.low_risk.model_id = low_risk_id;
selected_params.low_risk.APD90 = APD90(low_risk_id);
selected_params.low_risk.Ca_amp = Ca_amp(low_risk_id);
selected_params.low_risk.Ca_decay50 = Ca_decay50(low_risk_id);
selected_params.low_risk.g_CaL = param_sets(low_risk_id,1);
selected_params.low_risk.kNaCa = param_sets(low_risk_id,2);
selected_params.low_risk.VmaxUp = param_sets(low_risk_id,3);
selected_params.low_risk.g_irel_max = param_sets(low_risk_id,4);

%% Display the selected parameters clearly
fprintf('\n--- Final High-risk parameter set ---\n');
disp(selected_params.high_risk);

fprintf('\n--- Final Low-risk parameter set ---\n');
disp(selected_params.low_risk);

%% save results clearly
save('results/final_selected_params.mat','selected_params');

fprintf('Selected parameters saved to results/final_selected_params.mat\n');
