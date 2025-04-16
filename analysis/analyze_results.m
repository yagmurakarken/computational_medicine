%% analyze_ionic_currents.m 
clc; clear; close all;

%% Load saved data 
ventricular = load('results/ventricular_results.mat');
atrial = load('results/atrial_results.mat');
nodal = load('results/nodal_results.mat');
mixed = load('results/mixed_results.mat');

%% Define analysis window (last 5 s)
time_window = [795, 800];

%% Ionic and Calcium currents of interest explicitly
current_names = {'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel'};
phenotypes = {'Ventricular', 'Atrial', 'Nodal', 'Mixed'};
data = {ventricular, atrial, nodal, mixed};

%% Extract currents for each phenotype within analysis window
currents = struct();

for p = 1:length(phenotypes)
    sim_data = data{p};
    t = sim_data.t;

    idx_analysis = (t >= time_window(1)) & (t <= time_window(2));
    currents.(phenotypes{p}).time = t(idx_analysis);

    if strcmp(phenotypes{p}, 'Mixed')
        phenotype_currents = zeros(sum(idx_analysis), length(current_names));
        for c = 1:length(current_names)
            phenotype_currents(:, c) = sim_data.mixed_currents.(current_names{c})(idx_analysis);
        end
    else
        phenotype_currents = zeros(sum(idx_analysis), length(current_names));
        for c = 1:length(current_names)
            phenotype_currents(:, c) = sim_data.(current_names{c})(idx_analysis);
        end
    end

    currents.(phenotypes{p}).currents = phenotype_currents;
end

%% plot Ionic and Calcium Currents
figure('Name', 'Detailed Ionic and Calcium Currents Comparison', 'Position', [100, 100, 1400, 900]);

for c = 1:length(current_names)
    subplot(4, 2, c);
    hold on;

    % plotting for each phenotype
    for p = 1:length(phenotypes)
        plot(currents.(phenotypes{p}).time, currents.(phenotypes{p}).currents(:, c), ...
            'LineWidth', 1.2, 'DisplayName', phenotypes{p});
    end

    title(sprintf('Current: %s', current_names{c}));
    xlabel('Time (s)');
    ylabel('Current (A/F)');
    xlim(time_window);
    legend('Location', 'best');
    grid on;
end

sgtitle('Ionic and Calcium Current Analysis (Last 5 s of Simulation)');

%% save figure
saveas(gcf, 'results/ionic_calcium_current_comparison.png');