clc; clear; close all;

%% Add path explicitly
addpath('cellular_models/');

%% Simulation settings explicitly defined
options = odeset('MaxStep',1e-3,'InitialStep',2e-5);
stimFlag = 0;
tDrugApplication = 10000;
phenotype = 'ventricular';
t_fixed = 0:0.001:800;

%% Initial conditions
Y0 = [-0.070, 0.32, 0.0002, 0, 0, 1, 1, 1, 0, 1, 0, ...
      0.75, 0.75, 0, 0.1, 1, 0, 9.2, 0, 0.75, 0.3, 0.9, 0.1];

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

%% Output directory creation
output_dir = 'results/individual_models/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Run simulations with parfor and save individually
parfor model_idx = 1:numModels
    fprintf('Running model %d/%d...\n', model_idx, numModels);

    % Temporary parameter structure
    temp_params = struct('g_CaL', param_sets(model_idx,1), ...
                         'kNaCa', param_sets(model_idx,2), ...
                         'VmaxUp', param_sets(model_idx,3), ...
                         'g_irel_max', param_sets(model_idx,4));

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

    % Save results individually
    filename = sprintf('%s/model_%03d.mat', output_dir, model_idx);
    parsave_model(filename, t_fixed, temp_params, Vm, Cai, ICaL, INaCa, If, IKr, IKs, IK1, Iup, Irel);
end

fprintf('PoM generation completed successfully.\n');

function parsave_model(filename, t_fixed, params, Vm, Cai, ICaL, INaCa, If, IKr, IKs, IK1, Iup, Irel)
    save(filename, 't_fixed', 'params', 'Vm', 'Cai', 'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel', '-v7.3');
end

%% Visualization from individual files explicitly
Vm_all = zeros(length(t_fixed), numModels);

for model_idx = 1:numModels
    data = load(sprintf('%s/model_%03d.mat', output_dir, model_idx));
    Vm_all(:,model_idx) = data.Vm;
end

figure;
plot(t_fixed, Vm_all * 1e3, 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('Vm (mV)');
title(sprintf('Population of %d Ventricular Models (LHS ±30%%)', numModels));
grid on;
xlim([790 800]);
saveas(gcf, 'results/Population_of_Ventricular_Models.png');

fprintf('Visualization completed successfully.\n');
