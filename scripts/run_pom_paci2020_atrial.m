% run_pom_paci2020_atrial.m
% Population-of-Models (PoM) simulation using Paci2020 model for Atrial Phenotype

clc; clear; close all;

%% Add path explicitly
addpath('cellular_models/');

%% Simulation settings
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);
stimFlag = 0;
tDrugApplication = 10000;
phenotype = 'atrial'; % Changed to atrial phenotype
t_fixed = 0:0.001:800;

%% Initial conditions (same as previous, standard initial conditions)
load('Y0_atrial.mat');  % Contains variable Y0_atrial
Y0 = Y0_atrial;

%% Number of models
numModels = 200;

%% Parameter variation setup (±30%)
param_names = {'g_CaL', 'kNaCa', 'VmaxUp', 'g_irel_max'};
baseline_values = [8.635702e-5, 6514.47574, 0.82205, 55.808061];
variation_range = 0.3; % ±30%

%% Generate LHS
lhs_matrix = lhsdesign(numModels, length(param_names));
param_sets = zeros(numModels, length(param_names));

for i = 1:length(param_names)
    lower_bound = baseline_values(i)*(1 - variation_range);
    upper_bound = baseline_values(i)*(1 + variation_range);
    param_sets(:,i) = lower_bound + lhs_matrix(:,i)*(upper_bound - lower_bound);
end

%% Output directory creation (specific for atrial models)
output_dir = 'results/atrial_models/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

parfor model_idx = 1:numModels
    filename = sprintf('%s/model_%03d.mat', output_dir, model_idx);

    if exist(filename, 'file')
        fprintf('Model %03d already exists. Skipping...\n', model_idx);
        continue;
    end

    fprintf('Running atrial model %d/%d...\n', model_idx, numModels);

    % Temporary parameter structure
    temp_params = struct('g_CaL', param_sets(model_idx,1), ...
                         'kNaCa', param_sets(model_idx,2), ...
                         'VmaxUp', param_sets(model_idx,3), ...
                         'g_irel_max', param_sets(model_idx,4));

    try
        [~, Yc] = ode15s(@Paci2020, t_fixed, Y0, options, stimFlag, tDrugApplication, ...
                         phenotype, 1, 1, 1, 1, 1, 1, 1, temp_params);

        Vm = Yc(:,1);
        Cai = Yc(:,3);

        ICaL = zeros(size(t_fixed));
        INaCa = zeros(size(t_fixed));
        If = zeros(size(t_fixed));
        IKr = zeros(size(t_fixed));
        IKs = zeros(size(t_fixed));
        IK1 = zeros(size(t_fixed));
        Iup = zeros(size(t_fixed));
        Irel = zeros(size(t_fixed));

        for i = 1:length(t_fixed)
            [~, dati] = Paci2020(t_fixed(i), Yc(i,:), stimFlag, tDrugApplication, phenotype, ...
                                 1, 1, 1, 1, 1, 1, 1, temp_params);

            ICaL(i)  = dati(3);
            If(i)    = dati(2);
            IKr(i)   = dati(6);
            IKs(i)   = dati(5);
            IK1(i)   = dati(7);
            INaCa(i) = dati(8);
            Irel(i)  = dati(13);
            Iup(i)   = dati(14);
        end

        % Save results
        parsave_model(filename, t_fixed, temp_params, Vm, Cai, ICaL, INaCa, If, IKr, IKs, IK1, Iup, Irel);

    catch ME
        warning('Model %d failed: %s', model_idx, ME.message);
    end
end


fprintf('Atrial PoM generation completed successfully.\n');

function parsave_model(filename, t_fixed, params, Vm, Cai, ICaL, INaCa, If, IKr, IKs, IK1, Iup, Irel)
    save(filename, 't_fixed', 'params', 'Vm', 'Cai', 'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel', '-v7.3');
end

%% Visualization from individual files (atrial-specific)
Vm_all = zeros(length(t_fixed), numModels);

for model_idx = 1:numModels
    data = load(sprintf('%s/model_%03d.mat', output_dir, model_idx));
    Vm_all(:,model_idx) = data.Vm;
end

figure;
plot(t_fixed, Vm_all * 1e3, 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('Vm (mV)');
title(sprintf('Population of %d Atrial Models (LHS ±30%%)', numModels));
grid on;
xlim([790 800]);
results_dir = 'results/';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
saveas(gcf, fullfile(results_dir, 'Population_of_Atrial_Models.png'));

fprintf('Visualization completed successfully.\n');
