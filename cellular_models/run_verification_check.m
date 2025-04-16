clc; clear; close all;

%% Simulation Settings
phenotype = 'ventricular';
stimFlag = 0;
tDrugApplication = 10000;
t_fixed = 0:0.001:800;

Y0 = [-0.070, 0.32, 0.0002, 0, 0, 1, 1, 1, 0, 1, 0, ...
      0.75, 0.75, 0, 0.1, 1, 0, 9.2, 0, 0.75, 0.3, 0.9, 0.1];

%% Baseline parameters selected from published values
baseline_params = struct('g_CaL', 8.635702e-5, ...
                         'kNaCa', 6514.47574, ...
                         'VmaxUp', 0.82205, ...
                         'g_irel_max', 55.808061);

%% Solvers and tolerances for verification
solvers = {'ode15s', 'ode23t'};
tolerances = [1e-4, 1e-5, 1e-6];

results = struct();

%% Run convergence check simulations
for s = 1:length(solvers)
    solver = solvers{s};
    for tol = 1:length(tolerances)
        % Explicitly format tolerance string correctly
        tol_str = sprintf('%.0e', tolerances(tol));
        tol_str = strrep(tol_str, '-', 'm');  % replace minus sign to 'm'

        options = odeset('MaxStep', 1e-3, 'InitialStep', 2e-5, ...
                         'RelTol', tolerances(tol), 'AbsTol', tolerances(tol));

        fprintf('Running %s with tolerance %g...\n', solver, tolerances(tol));

        solver_func = str2func(solver);
        [~, Y] = solver_func(@Paci2020, t_fixed, Y0, options, stimFlag, ...
                             tDrugApplication, phenotype, 1,1,1,1,1,1,1, baseline_params);

        results.(solver).(['tol_', tol_str]) = Y(:,1); % Vm only
    end
end

save('results/verification_results.mat', 'results', 't_fixed');

%% Visualization for verification (Vm at final seconds)
figure;
hold on;
colors = lines(length(solvers)*length(tolerances));
legend_entries = {};

counter = 1;
for s = 1:length(solvers)
    solver = solvers{s};
    for tol = 1:length(tolerances)
        % Make sure you use the same explicit formatting again here
        tol_str = sprintf('%.0e', tolerances(tol));
        tol_str = strrep(tol_str, '-', 'm');

        Vm = results.(solver).(['tol_', tol_str]) * 1e3;
        plot(t_fixed, Vm, 'Color', colors(counter,:), 'LineWidth', 1.2);
        legend_entries{counter} = sprintf('%s, tol=%g', solver, tolerances(tol));
        counter = counter + 1;
    end
end

xlabel('Time (s)');
ylabel('Vm (mV)');
title('Verification: Convergence and Solver Comparison');
legend(legend_entries, 'Location', 'best');
grid on;
xlim([795 800]);
saveas(gcf, 'results/Verification_Convergence_Check.png');

fprintf('Verification completed successfully.\n');
