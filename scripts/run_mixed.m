%% run_mixed.m 
% Mixed (80% ventricular, 10% atrial, 10% nodal) phenotype simulation

clc; clear; close all;

%% load individual phenotype results
ventricular_data = load('results/ventricular_results.mat');
atrial_data      = load('results/atrial_results.mat');
nodal_data       = load('results/nodal_results.mat');

t = ventricular_data.t; % assuming identical 't'

%% Extract state variables 
Yc_ventricular = ventricular_data.Yc;
Yc_atrial      = atrial_data.Yc;
Yc_nodal       = nodal_data.Yc;

%% Extract ionic and calcium currents 
currents = {'INa', 'ICaL', 'INaCa', 'If', 'IKr', 'IKs', 'IK1', 'Iup', 'Irel'};
for i = 1:length(currents)
    current_name = currents{i};
    ventricular_currents.(current_name) = ventricular_data.(current_name);
    atrial_currents.(current_name)      = atrial_data.(current_name);
    nodal_currents.(current_name)       = nodal_data.(current_name);
end

%% dimension check
if ~isequal(size(Yc_ventricular), size(Yc_atrial), size(Yc_nodal))
    error('dimension mismatch between phenotype simulations.');
end

%% Phenotype mixing proportions
ventricular_fraction = 0.80;
atrial_fraction      = 0.10;
nodal_fraction       = 0.10;

%% Combine state variables 
Yc_mixed = ventricular_fraction * Yc_ventricular + ...
           atrial_fraction * Yc_atrial + ...
           nodal_fraction * Yc_nodal;

%% Combine ionic and calcium currents
for i = 1:length(currents)
    current_name = currents{i};
    mixed_currents.(current_name) = ventricular_fraction * ventricular_currents.(current_name) + ...
                                    atrial_fraction * atrial_currents.(current_name) + ...
                                    nodal_fraction * nodal_currents.(current_name);
end

%% Save mixed results
save('results/mixed_results.mat', 't', 'Yc_mixed', 'mixed_currents');

%% Quick verification
figure;
plot(t, Yc_mixed(:,1)*1e3, 'LineWidth', 1.2); % Vm in mV
xlabel('Time (s)');
ylabel('Vm (mV)');
title('Mixed hPSC-CM AP (80% ventricular, 10% atrial, 10% nodal)');
grid on;
xlim([790 800]);

%% Biomarker Calculation for mixed phenotype
analysis_window = (t >= 780 & t <= 800);
t_analysis  = t(analysis_window);
Vm_analysis = Yc_mixed(analysis_window,1)*1e3;

peak_indices = find(diff(sign(diff(Vm_analysis))) < 0) + 1;
peak_threshold = 0;
valid_peaks = Vm_analysis(peak_indices) > peak_threshold;
pks  = Vm_analysis(peak_indices(valid_peaks));
locs = t_analysis(peak_indices(valid_peaks));

spont_intervals = diff(locs);
mean_spont_interval = mean(spont_intervals);

DMP = zeros(length(locs)-1, 1);
for i = 1:length(locs)-1
    segment = (t_analysis >= locs(i)) & (t_analysis <= locs(i+1));
    DMP(i) = min(Vm_analysis(segment));
end
mean_DMP = mean(DMP);
mean_peak_voltage = mean(pks);

APD90 = zeros(length(locs)-1, 1);
for i = 1:length(locs)-1
    AP_segment_time = t_analysis((t_analysis >= locs(i)) & (t_analysis <= locs(i)+mean_spont_interval));
    AP_segment_Vm = Vm_analysis((t_analysis >= locs(i)) & (t_analysis <= locs(i)+mean_spont_interval));

    peak_voltage = max(AP_segment_Vm);
    dmp_voltage  = min(AP_segment_Vm);
    repol_90_voltage = peak_voltage - 0.9*(peak_voltage - dmp_voltage);

    post_peak_indices = find(AP_segment_time > locs(i));
    repol_cross_idx = find(AP_segment_Vm(post_peak_indices) <= repol_90_voltage, 1, 'first');

    if ~isempty(repol_cross_idx)
        APD90(i) = (AP_segment_time(post_peak_indices(repol_cross_idx)) - locs(i))*1000;
    else
        APD90(i) = NaN;
    end
end
mean_APD90 = nanmean(APD90);

%% Display computed mixed biomarkers
fprintf('Mixed Biomarkers (80%% ventricular, 10%% atrial, 10%% nodal):\n');
fprintf('Diastolic Membrane Potential (DMP): %.2f mV\n', mean_DMP);
fprintf('Peak Voltage: %.2f mV\n', mean_peak_voltage);
fprintf('APD90: %.2f ms\n', mean_APD90);
fprintf('Spontaneous Beating Interval: %.2f s\n', mean_spont_interval);

disp('Mixed simulation completed and successfully.');