%% quantitative_ionic_analysis.m
clc; clear; close all;

phenotypes = {'ventricular', 'atrial', 'nodal'};
current_names = {'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel'};
current_indices = [3, 8, 2, 6, 5, 7, 14, 13];
time_window = [795, 800];

fraction_weights = struct('ventricular', 0.8, 'atrial', 0.1, 'nodal', 0.1);

results = struct();

%% Load phenotype data
phenotype_currents = struct();

for p = 1:length(phenotypes)
    sim_file = ['results/' phenotypes{p} '_results.mat'];
    sim_data = load(sim_file);

    Yc = sim_data.Yc;
    t = sim_data.t;
    idx_analysis = (t >= time_window(1)) & (t <= time_window(2));

    analysis_time = t(idx_analysis);
    currents = zeros(sum(idx_analysis), length(current_names));

    for i = find(idx_analysis)'
        [~, dati] = Paci2020(t(i), Yc(i,:), 0, 10000, phenotypes{p}, 1,1,1,1,1,1,1);
        currents(i - find(idx_analysis,1) + 1, :) = dati(current_indices);
    end

    phenotype_currents.(phenotypes{p}) = currents;

    % Individual phenotype quantitative analysis
    for c = 1:length(current_names)
        curr = currents(:, c);
        avg_amplitude = mean(abs(curr));
        auc = trapz(analysis_time, abs(curr));
        [~, peak_idx] = max(abs(curr));
        peak_time = analysis_time(peak_idx);

        results.(phenotypes{p}).(current_names{c}) = struct(...
            'AverageAmplitude', avg_amplitude, ...
            'AreaUnderCurve', auc, ...
            'PeakTiming', peak_time ...
        );
    end
end

%% mixed currents analysis
currents_mixed = fraction_weights.ventricular * phenotype_currents.ventricular + ...
                 fraction_weights.atrial * phenotype_currents.atrial + ...
                 fraction_weights.nodal * phenotype_currents.nodal;

for c = 1:length(current_names)
    curr = currents_mixed(:, c);
    avg_amplitude = mean(abs(curr));
    auc = trapz(analysis_time, abs(curr));
    [~, peak_idx] = max(abs(curr));
    peak_time = analysis_time(peak_idx);

    results.mixed.(current_names{c}) = struct(...
        'AverageAmplitude', avg_amplitude, ...
        'AreaUnderCurve', auc, ...
        'PeakTiming', peak_time ...
    );
end

%% Save computed results
save('results/ionic_calcium_currents_quantitative.mat', 'results');

%% Results Display
fprintf('Quantitative Ionic and Calcium Current Analysis:\n');
for phenotype = [{'ventricular'}, {'atrial'}, {'nodal'}, {'mixed'}]
    pheno_name = phenotype{1};
    fprintf('\nPhenotype: %s\n', pheno_name);
    for c = 1:length(current_names)
        r = results.(pheno_name).(current_names{c});
        fprintf('%s - Avg: %.4f A/F, AUC: %.4f, Peak at %.3fs\n', ...
                current_names{c}, r.AverageAmplitude, r.AreaUnderCurve, r.PeakTiming);
    end
end