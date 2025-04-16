clc; clear;
dt = 0.01e-3; % 0.01 ms timestep (same as tissue)
t_end = 10;    % 5 seconds for pacing steady state
t = 0:dt:t_end;
Y = zeros(23, length(t));

% Use provided steady-state initial conditions for paced scenario
Y(:,1) = [-0.0728; 0.0914; 1.4387e-05; 0; 0.00011; 0.9493; 0.99998; 0.99914; 0.0382; 0.4249; ...
          0.0354; 0.7188; 0.6906; 0.0464; 0.1458; 0.8226; 0.00621; 7.4758; 0.00334; 0.1153; ...
          0.03332; 7.2334e-05; 0.9086];

stim_period = 2;      % 1 Hz pacing
stim_duration = 1e-3; % 1 ms stimulus
stim_amplitude = -80e-12;

for step = 1:length(t)-1
    current_time = t(step);
    stimFlag = double(mod(current_time, stim_period) <= stim_duration);

    [dY, ~] = Paci2020(current_time, Y(:,step), stimFlag, 0, 'ventricular', 1,1,1,1,1,1,1, struct());

    gating_vars = [5:17, 19:23];
    non_gating_vars = setdiff(1:23, gating_vars);

    % Rush-Larsen (gating)
    Y(gating_vars,step+1) = Y(gating_vars,step) + dt*dY(gating_vars);

    % Euler (others)
    Y(non_gating_vars,step+1) = Y(non_gating_vars,step) + dt*dY(non_gating_vars);
end

Vm = Y(1,:)*1000; % convert to mV
time = t*1000; % ms

plot(time, Vm);
xlabel('Time (ms)'); ylabel('Voltage (mV)');
title('Single-cell Paci2020 (explicit Rush-Larsen + Euler)');

%% Explicit calculation of APD90 (single cell) - Improved
Vm_mV = Vm; % mV
time_ms = time; % ms

peak_threshold = 0;
[peaks, locs] = findpeaks(Vm_mV, time_ms, 'MinPeakHeight', peak_threshold);

if numel(peaks) < 1
    error('No valid peaks detected.');
end

last_peak_time = locs(end);
AP_start_idx = find(time_ms >= last_peak_time, 1, 'first');
analysis_end_idx = find(time_ms >= last_peak_time + 3000, 1, 'first'); % extended 3 s window

% Find diastolic potential (minimum after peak)
diastolic_potential = min(Vm_mV(AP_start_idx:analysis_end_idx));

% Compute 90% repolarization voltage
V_90 = peaks(end) - 0.9*(peaks(end) - diastolic_potential);

% Find APD90 time index (after peak)
post_peak_Vm = Vm_mV(AP_start_idx:analysis_end_idx);
post_peak_time = time_ms(AP_start_idx:analysis_end_idx);
cross_idx = find(post_peak_Vm <= V_90, 1, 'first');

if ~isempty(cross_idx)
    APD90_single_cell = post_peak_time(cross_idx) - last_peak_time;
    fprintf('Explicit Single-Cell APD90: %.2f ms\n', APD90_single_cell);
else
    disp('Could not determine APD90: insufficient repolarization within analysis window (consider increasing simulation duration).');
end
