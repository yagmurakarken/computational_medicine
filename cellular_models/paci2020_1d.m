%% paci2020_1d_simulation.m - explicit 1D cable simulation using Paci2020 model
clear; clc;

dt = 0.01e-3; % timestep (10 µs)
t_end = 0.8; % at least 500 ms
N_cells = 10; % for quicker test
stim_period = 1; % stimulus every 1 s
stim_duration = 1e-3; % 1 ms
stim_amplitude = -80e-12; % correct amplitude

stim_cells = 1:2; % cells receiving stimulation

%% Initialize state variables
Y = zeros(23, N_cells);
[Y(:, :), ~] = paci_initial_conditions(N_cells);

%% Diffusion parameters
D = 0.001; % lowered diffusion (cm^2/ms)
cell_length = 0.01; % cm


%% Time loop
n_steps = round(t_end / dt);
time = 0;
V = zeros(N_cells, n_steps); % voltage over time

for step = 1:n_steps
    time = time + dt;

    % Solve cell models
    for idx = 1:N_cells
        isStim = double(mod(time, stim_period) <= stim_duration && ismember(idx, stim_cells));
        [dY_temp, ~] = Paci2020(time, Y(:, idx), isStim, 0, 'ventricular', ...
                                1, 1, 1, 1, 1, 1, 1, struct());

        % Rush-Larsen for gating variables
        gating_vars = [5:17, 19:23];
        Y(gating_vars, idx) = Y(gating_vars, idx) + dt * dY_temp(gating_vars);

        % Euler for remaining variables
        non_gating_vars = setdiff(1:23, gating_vars);
        Y(non_gating_vars, idx) = Y(non_gating_vars, idx) + dt * dY_temp(non_gating_vars);
    end

    % Diffusion: explicit scheme for membrane voltage
    Vm = Y(1, :);
    d2Vm_dx2 = zeros(1, N_cells);
    d2Vm_dx2(2:end-1) = (Vm(3:end) - 2*Vm(2:end-1) + Vm(1:end-2)) / cell_length^2;
    d2Vm_dx2(1) = (Vm(2) - Vm(1)) / cell_length^2;           % left boundary
    d2Vm_dx2(end) = (Vm(end-1) - Vm(end)) / cell_length^2;   % right boundary

    % Voltage update
    Y(1, :) = Y(1, :) + dt * D * d2Vm_dx2;
    Y(1,:) = real(Y(1,:)); % remove imaginary parts
    Y(1,:) = max(min(Y(1,:), 0.05), -0.1); % clamp voltages between -100 mV and +50 mV


    % Store voltage
    V(:, step) = Y(1, :);

    % Progress update
    if mod(step, round(n_steps/10)) == 0
        fprintf('Simulation progress: %.1f%%\n', (step/n_steps)*100);
    end
end

%% Plot voltage traces for first 5 cells
figure; hold on;
time_vec = (1:n_steps)*dt*1000; % ms
for i = 1:5
    plot(time_vec, V(i, :)*1000, 'DisplayName', ['Cell ' num2str(i)]);
end
xlabel('Time (ms)'); ylabel('Membrane potential (mV)');
title('Voltage Traces for First 5 Cells');
legend;

%% Activation time analysis
activation_threshold = 0; % mV
activation_times = nan(1, N_cells);

for cell = 1:N_cells
    idx = find(V(cell, :) >= activation_threshold/1000, 1, 'first');
    if ~isempty(idx)
        activation_times(cell) = idx * dt;
    end
end
disp('Activation times (ms):');
disp(activation_times * 1000);

% Visualize activation time vs. position
cell_positions = (0:N_cells-1) * cell_length;
figure;
plot(cell_positions, activation_times*1000, 'o-');
xlabel('Cell position (cm)');
ylabel('Activation time (ms)');
title('Activation Time vs Cell Position');

% Calculate conduction velocity
valid = ~isnan(activation_times);
if sum(valid) >= 2
    p = polyfit(cell_positions(valid), activation_times(valid), 1);
    CV = 1 / p(1); % cm/s
    fprintf('Conduction velocity: %.3f cm/s\n', CV);
    fprintf('Conduction velocity: %.3f m/s\n', CV / 100);
else
    warning('Not enough data to calculate CV.');
end

%% APD90 calculation
resting_potential = -80e-3;
peak_potential = max(V(1,:));
threshold_90 = peak_potential - 0.9*(peak_potential - resting_potential);

cell_to_check = 1;
crossing_up = find(V(cell_to_check,:) >= 0, 1, 'first');
crossing_down = find(V(cell_to_check, crossing_up:end) <= threshold_90, 1, 'first') + crossing_up - 1;

if ~isempty(crossing_up) && ~isempty(crossing_down)
    APD90 = (crossing_down - crossing_up)*dt * 1000;
    fprintf('APD90: %.2f ms\n', APD90);
else
    disp('APD90 could not be calculated: insufficient repolarization.');
end

%% Helper function
function [Y0, dati0] = paci_initial_conditions(N_cells)
    Y0_single = [-0.0750; 0.3; 0.0002; 0; 0; 1; 1; 1; 0; 1; 0; 0.75; 0.75; 0; 0; 0; 0; 10; 0; 1; 0; 0; 1];
    Y0 = repmat(Y0_single, 1, N_cells);
    dati0 = [];
end

