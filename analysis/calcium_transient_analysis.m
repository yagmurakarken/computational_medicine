%% calcium_transient_analysis.m
clc; clear; close all;

phenotypes = {'ventricular', 'atrial', 'nodal', 'mixed'};
results = struct();
time_window = [795, 800];

figure('Position', [100, 100, 1200, 800]);

for p = 1:length(phenotypes)
    data = load(['results/', phenotypes{p}, '_results.mat']);
    t = data.t;

    if isfield(data, 'Yc')
        Yc = data.Yc;
    else
        Yc = data.Yc_mixed;
    end

    Cai = Yc(:,3) * 1e6; % Convert to µM
    idx_analysis = (t >= time_window(1)) & (t <= time_window(2));
    t_analysis = t(idx_analysis);
    Cai_analysis = Cai(idx_analysis);

    % transient metrics calculation
    Cai_peak = max(Cai_analysis);
    Cai_baseline = min(Cai_analysis);
    transient_amplitude = Cai_peak - Cai_baseline;

    [~, peak_idx] = max(Cai_analysis);
    decay_target = Cai_peak - 0.5 * transient_amplitude; % 50% decay
    decay_idx = find(Cai_analysis(peak_idx:end) <= decay_target, 1) + peak_idx - 1;
    decay_time = t_analysis(decay_idx) - t_analysis(peak_idx);

    results.(phenotypes{p}) = struct('PeakCai', Cai_peak, ...
                                     'BaselineCai', Cai_baseline, ...
                                     'Amplitude', transient_amplitude, ...
                                     'DecayTime50', decay_time);

    % plotting
    subplot(2, 2, p);
    plot(t_analysis, Cai_analysis, 'LineWidth', 1.5);
    title(sprintf('%s Phenotype Ca^{2+} Transient', capitalize(phenotypes{p})));
    xlabel('Time (s)');
    ylabel('[Ca^{2+}]_i (\muM)');
    grid on;
end

%% Save results
save('results/calcium_transients.mat', 'results');

%% Display calculated metrics
fprintf('\nCalcium Transient Metrics:\n');
for p = 1:length(phenotypes)
    res = results.(phenotypes{p});
    fprintf('%s: Peak=%.2f µM, Baseline=%.2f µM, Amplitude=%.2f µM, 50%% Decay=%.3f s\n',...
        capitalize(phenotypes{p}), res.PeakCai, res.BaselineCai, res.Amplitude, res.DecayTime50);
end

%% save figure
saveas(gcf, 'results/calcium_transient_comparison.png');

function str = capitalize(str)
    str(1) = upper(str(1));
end